% This routine is about solving equation (75)
function fncon = frconim(m,sl,rcon,theta,bendsti,...
    bendratio,eps_ch,consta)
% m - number of points on vesicle
% sl - length of vesicle
% rcon - concentration species
% theta - opening angle
% bendsti - bending stiffness
% bendratio - bending ratio
% eps_ch - small parameter in the flux in equation (65)
% consta - constant a in equation (65) and other places

% fncon - The non-linear N_2 defined in equation (67)
 
% setting up RHS of eq (68) 
% termd is the a/eps(f'(u)-eps^2u_ss) + b'(u)/2 * kappa^2 terms in eq (13) 
% ie. termd is the variatiaonl derivative of the energy with
% respect to the lipid concentration
termd = fluxj(sl,theta,rcon,bendsti,bendratio,consta,m,eps_ch);
% now taking derivative of the variational derivative twice
rcons = fd1(termd,m);       
rconss = fd1(rcons,m);
% betarcon is the 1/Pe term in eq (67)
betarcon = 1;
% fcon is the part of N_2 in eq (67) but is still missing the
% a\eps/s_al^4 * u_alalalal term, al is shorthand for alpha here
fcon = betarcon*rconss(1:m)/sl^2;
%fcon(m+1) = fcon(1);

fconHat = fft(fcon);
rconHat = fft(rcon,m);
% N are the fourier modes
N = pi*2*[0 1:m/2 m/2-1:-1:1];
% rlen is the term that needs to multiply the fourier coefficients of
% the lipid species to result in the fourth derivative scaled by the a
% and \eps constants in equation (67)
rlen = (N/sl).^4*consta*eps_ch;
% N2Hat is the Fourier coefficeints of the right hand side of
% equation (67) 
N2Hat = fconHat + rlen.*rconHat;
% Move back to physical (real) space
fncon = real(ifft(N2Hat));

end
