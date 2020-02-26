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
 
%setting up RHS of eq (75) 
%computing the variational derivative in eq (13) only for
% termd is the a/e(f'(u)-e^2u_ss) + b'(u)/2 * k^2 terms in eq (13) 
termd = fluxj(sl,theta,rcon,bendsti,bendratio,consta,m,eps_ch);
% now taking derivative of the variational derivative twice
rcons = fd1(termd,m);       
rconss = fd1(rcons,m);
%betarcon is the 1/pe term in eq (75)?
betarcon = 1;
%fcon is u_t in eq (75)
fcon(1:m) = betarcon*rconss(1:m)/sl^2;
fcon(m+1) = fcon(1);

% now we must subtract off the stiffest part in Fourier space. (???)

%temp is u_t in eq (75)
temp(1:m) = fcon(1:m);
%tempt is u, input parameter
tempt(1:m) = rcon(1:m);

temp = [temp temp(1)];
tempt = [tempt tempt(1)];

%temp is now the fourier transform of u_t
% this is a fourth order derivative
temp = fft(temp,m);
%tempt is now the fourier transform of u
tempt = fft(tempt,m);
%N are the fourier modes
N=pi*2*[0 1:m/2 m/2-1:-1:1];
%rlen 
rlen=(N/sl).^4*consta*eps_ch;

temp4 = temp(1,1:m)+rlen(1,1:m).*tempt(1,1:m);

temp4 = real(ifft(temp4,m));

fncon(1,1:m) = temp4(1,1:m);

%fncon(1,m+1) = fncon(1,1);

end
