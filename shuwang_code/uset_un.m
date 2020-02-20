% this code forms the normal and tan velocity given the angle and
% arclength and returns the normal and tan vel. from the ave. vel.
% integral
function [un,vn] = uset_un(sl,theta,rcon,bendsti,bendratio,m)

% Calculate curvature first
dkap = acurv(sl,theta,m);

b0 = bendsti;
b1 = bendsti*bendratio;

% bending coefficient which depends on the lipid concentration that is
% stored in rcon. This is the variable b(u) in equation (10)
rbn(1,1:m) = b0*(ones(1,m) - rcon(1,1:m)) + b1*rcon(1,1:m);

% Note that dkap is not the derivative of the curvature, but is the
% curvature. The variable name is most likely left over from Shuwang's
% Fortran code ie. d = double
cv0(1,1:m) = rbn(1,1:m).*(dkap(1,1:m)); 
cv1 = fd1(cv0,m);
cv2 = fd1(cv1,m);

% derivative of the bending coefficient
bs = fd1(rbn,m);

% compute normal velocity-unconstrained from bending
% un variable is equation (14) with spotaneous curvature set to zero.
% The sl^2 term is needed since fd1 computes derivatives w.r.t. [0,2*pi]
% rather than w.r.t. arclength. Therefore, each s derivative requires
% dividing by the total length of the vesicle (this also relies on
% equispaced points in arclength on the vesicle)
un(1,1:m) = -(cv2(1,1:m)/sl^2+...
    rbn(1,1:m)/2.*dkap(1,1:m).*(dkap(1,1:m).*dkap(1,1:m)));
% vn variable is the second term in equation (13) (differs by a negative
% sign). The last term drops since spontaneous curvature is 0. The first
% term is not in this routine since we are only calculating variations
% due to changes in the vesicle shape and not the lipid species (yet).
vn(1,1:m) = -bs(1,1:m).*dkap(1,1:m).^2/sl/2;

un(1,m+1) = un(1,1);
vn(1,m+1) = vn(1,1);

end
