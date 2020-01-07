% this code forms the normal and tan velocity given the angle and
% arclength and returns the normal and tan vel. from the ave. vel.
% integral
function [un,vn] = uset_un(sl,theta,rcon,bendsti,bendratio,m)

% Calculate curvature first
dkap = acurv(sl,theta,m);

b0 = bendsti;
b1 = bendsti*bendratio;

% bending coefficient
rbn(1,1:m) = b0*(ones(1,m) - rcon(1,1:m)) + b1*rcon(1,1:m);

% Note that dkap is not the derivative of the curvature, but is the
% curvature. The variable name is most likely left over from Shuwang's
% Fortran code ie. d = double
cv0(1,1:m) = rbn(1,1:m).*(dkap(1,1:m)); 
cv1 = fd1(cv0,m);
cv2 = fd1(cv1,m);

% derivative of the bending coefficient
bs = fd1(rbn,m);

% compute normal velocity-unconstraint from bending
un(1,1:m) = -(cv2(1,1:m)/sl^2+...
    rbn(1,1:m)/2.*dkap(1,1:m).*(dkap(1,1:m).*dkap(1,1:m)));
vn(1,1:m) = -bs(1,1:m).*dkap(1,1:m).^2/sl/2;

un(1,m+1) = un(1,1);
vn(1,m+1) = vn(1,1);

end
