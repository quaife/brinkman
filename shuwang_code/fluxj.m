% THIS ROUTINE COMPUTES THE FORCING ON THE RHS flux OF THE Cahn_Hillard
% Equation for rcon periodic on [0,1] and for closed curves
function term = fluxj(sl,theta,rcon,bendsti,bendratio,consta,m,eps_ch)
% sl - length of vesicle
% theta - opening angle
% rcon - concentration
% bendsti - bending stiffness
% bendratio - bending ratio
% consta - term 'a' in the paper. It is the line tension scaled by a
% characteristic bending stiffness
% m - number of points on vesicle
% eps_ch - small parameter in equation (65)

%term - ???

% compute the 1st term from free energy
% term1 is f'(u) where f(u) is defined immediately after equation (11).
% This shows up in equation (13)
% It is the derivative of the double-well potential
term1 = 1/2*(rcon(1,1:m).*(ones(1,m)-rcon(1,1:m)).^2 - ...
             rcon(1,1:m).^2.*(ones(1,m)-rcon(1,1:m)));

% compute the 2nd term
rcons = fd1(rcon,m); % derivative of the concentration
rconss = fd1(rcons,m); % second derivative of the concentration
% This is eps^2 * u_{ss} as defined in equation (13)
term2 = -eps_ch^2*rconss(1,1:m)/sl^2;

%compute the 3rd term
b0 = bendsti;
b1 = bendsti*bendratio;
% rbn is not even used after this.
rbn = b0*(ones(1,m) - rcon(1,1:m)) + b1*rcon(1,1:m);
rbndu = (b1 - b0)*ones(1,m); 
% difference between the bending ratios
% Taking a difference seems really strange
rbn = [rbn rbn(1)];

% compute the curvature
dkap = acurv(sl,theta,m);

term3 = 0.5*rbndu.*dkap(1,1:m).^2;
% compute the 4th term compute partial kappa partial rcon here

% now put all terms together. Note that term3 is not multiplied by the
% a/eps term
term(1,1:m) = consta/eps_ch * (term1 + term2) + term3(1,1:m);
term = [term term(1)];

end
