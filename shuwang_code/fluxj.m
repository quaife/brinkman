% THIS ROUTINE COMPUTES THE FORCING ON THE RHS flux OF THE Cahn_Hillard
% Equation for rcon periodic on [0,1] and for closed curves
function term = fluxj(sl,theta,rcon,bendsti,bendratio,consta,m,eps_ch)

% compute the 1st term from free energy
term1(1,1:m)=1/2*(rcon(1,1:m).*(ones(1,m)-rcon(1,1:m)).^2 - ...
                  rcon(1,1:m).^2.*(ones(1,m)-rcon(1,1:m)));
term1(m+1) = term1(1);
% compute the 2nd term
rcons = fd1(rcon,m);
rconss = fd1(rcons,m);

term2(1,1:m) = -eps_ch^2*rconss(1,1:m)/sl^2;
term2(1,m+1) = term2(1,1);
%compute the 3rd term
b0 = bendsti;
b1 = bendsti*bendratio;
rbn(1,1:m) = b0*(ones(1,m)-rcon(1,1:m))+b1*rcon(1,1:m);
rbndu(1,1:m) = (b1-b0)*ones(1,m);
rbn(1,m+1) = rbn(1,1);

dkap = acurv(sl,theta,m);

term3(1:m) = rbndu(1,1:m)/2.*dkap(1,1:m).^2;
term3(1,m+1) = term3(1,1);
% compute the 4th term compute partial kappa partial rcon here

% now put all terms together
term(1,1:m) = consta/eps_ch*(term1(1,1:m)+term2(1,1:m))+term3(1,1:m);
term(1,m+1) = term(1,1);

end
