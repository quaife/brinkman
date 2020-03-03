% THIS ROUTINE COMPUTES THE FORCING ON THE RHS flux OF THE Cahn_Hillard
% Equation for rcon periodic on [0,1] and for closed curves
function term = fluxj(sl,theta,rcon,bendsti,bendratio,consta,m,eps_ch)
% sl: length of vesicle
% theta: opening angle
% rcon: concentration
% bendsti: bending stiffness
% bendratio: bending ratio
% consta: term 'a' in the paper. It is the line tension scaled by a
% characteristic bending stiffness
% m: number of points on vesicle
% eps_ch: small parameter in equation (65)

% term: a/eps(f'(u)-eps^2u_ss) + b'(u)/2 * kappa^2 terms in eq (13) 


% dfu is the derivative of the double-well potential
dfu = 1/2*(rcon(1,1:m).*(ones(1,m)-rcon(1,1:m)).^2 - ...
             rcon(1,1:m).^2.*(ones(1,m)-rcon(1,1:m)));
% derivative of the concentration
rcons = fd1(rcon,m); 
% second derivative of the concentration
rconss = fd1(rcons,m);
% term2 is eps^2 * u_ss as defined in equation (13)
term2 = -eps_ch^2*rconss(1,1:m)/sl^2;

% computing the b'(u)/2 * kappa^2 term in eq (13)
b0 = bendsti;
b1 = bendsti*bendratio;
% rbn is b(u).
rbn = b0*(ones(1,m) - rcon(1,1:m)) + b1*rcon(1,1:m);
% rbndu is b'(u) 
rbndu = (b1 - b0)*ones(1,m);
rbn = [rbn rbn(1)];
% compute the curvature
dkap = acurv(sl,theta,m);
% term3 is b'(u)/2 * kappa^2 in eq (13)
term3 = 0.5*rbndu.*dkap(1,1:m).^2;

% term is a/eps(f'(u)-eps^2u_ss) + b'(u)/2 * kappa^2 terms in eq (13) 
term(1,1:m) = consta/eps_ch * (dfu + term2) + term3(1,1:m);
term = [term term(1)];

end
