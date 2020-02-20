% This routine is about solving equation (75)
function [fncon,termd] = frconim(m,sl,rcon,theta,bendsti,...
    bendratio,eps_ch,consta)
% m - number of points on vesicle
% sl - length of vesicle
% rcon - concentration species
% theta - opening angle
% bendsti - bending stiffness
% bendratio - bending ratio
% eps_ch - small parameter in the flux in equation (65)
% consta - ???

% fncon - ???
% termd - ???
 
% compute the old forcing
termd = fluxj(sl,theta,rcon,bendsti,bendratio,consta,m,eps_ch);
% now taking derivative of term, twice
rcons = fd1(termd,m);       
rconss = fd1(rcons,m);
betarcon = 1;
fcon(1,1:m) = betarcon*rconss(1,1:m)/sl^2;
fcon(1,m+1) = fcon(1,1);

% now we must subtract off the stiffest part in Fourier space.
temp(1,1:m) = fcon(1,1:m);
tempt(1,1:m) = rcon(1,1:m);

temp(1,m+1) = temp(1,1);
tempt(1,m+1) = tempt(1,1);
temp = fft(temp,m);
tempt = fft(tempt,m);

% this is a fourth order derivative
N=pi*2*[0 1:m/2 m/2-1:-1:1];
rlen=(N/sl).^4*consta*eps_ch;

temp4 = temp(1,1:m)+rlen(1,1:m).*tempt(1,1:m);

temp4 = real(ifft(temp4,m));

fncon(1,1:m) = temp4(1,1:m);

%fncon(1,m+1) = fncon(1,1);

end
