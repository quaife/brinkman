function [un,vn] = uset_un(sl,theta,rcon,bendsti,bendratio,m)
% sl - length of vesicle
% theta - opening angle
% rcon - lipid concentration species
% bendsti - bending stiffness
% bendratio - bending ratio
% m - number of points on vesicle

% un - First two terms in Equation (14), but with spotaneous curvature
% set to 0
% vn - Last two terms in Equation (13), but with spontaneous curvature
% set to 0

% Calculate curvature first
dkap = acurv(sl,theta,m);
b0 = bendsti;
b1 = bendsti*bendratio;
% %disp('Plotting RCON')
% plot(dkap)
% %pause
% bending coefficient which depends on the lipid concentration that is
% stored in rcon. This is the variable b(u) in equation (10)
rbn = b0*(ones(1,m) - rcon(1:m)) + b1*rcon(1:m);
% plot(rbn)
% %pause
% Note that dkap is not the derivative of the curvature, but is the
% curvature. The variable name is most likely left over from Shuwang's
% Fortran code ie. d = double
cv0 = rbn.*dkap(1:m); % bending times curvature
cv1 = fd1(cv0,m); % derivative of bending times curvature
cv2 = fd1(cv1,m); % second derivative of bending times curvature
% plot(cv2)
% %pause
% derivative of the bending coefficient
bs = fd1(rbn,m);

%plot(dkap)
%plot(cv2(1:m)/sl^2)
%%pause

% un variable is equation (14) with spotaneous curvature set to zero.
% The sl^2 term is needed since fd1 computes derivatives w.r.t. [0,2*pi]
% rather than w.r.t. arclength. Therefore, each s derivative requires
% dividing by the total length of the vesicle (this also relies on
% equispaced points in arclength on the vesicle)
un = -(cv2(1:m)/sl^2 + rbn/2.*dkap(1,1:m).^3);
% %disp('non')
% plot(un)
% vn variable is the second term in equation (13) (differs by a negative
% sign). The last term drops since spontaneous curvature is 0. The first
% term is not in this routine since we are only calculating variations
% due to changes in the vesicle shape and not the lipid species (yet).
vn(1,1:m) = -bs(1,1:m).*dkap(1,1:m).^2/sl/2;
% %figure(2)
% plot(vn)
% %pause
% MIGHT COME FROM THE NEGATIVE SIGN IN EQUATION (23) WHICH HAS A
% NEGATIVE IN FRONT OF THE VARIATION W.R.T. $u$

un = [un un(1)];
vn = [vn vn(1)];

end
