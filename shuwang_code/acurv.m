% This routine computes the curvature from the tangent angle
% formulation

function dkap = acurv(sl,theta,m)

theta(1,m+1) = theta(1,1);
dkap = fd1x(theta(1,1:m+1),m)/sl;
dkap(1,m+1) = dkap(1,1);

end
