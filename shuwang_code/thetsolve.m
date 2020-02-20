function theta = thetsolve(sxn,syn,sl,t0,ngrid)
% sxn - first x derivative of vesicle
% syn - first y derivative of vesicle
% sl -  length of vesicle
% t0 - initial theta value
% ngrid - number of points on the vesicle
% theta - opening angle between tangent vector and horizontal

ssxn = fd1(sxn,ngrid);
ssyn = fd1(syn,ngrid);

dkap = (sxn.*ssyn-syn.*ssxn)/(sl^3);
% compute curvature
dkap(ngrid+1) = dkap(1);
stheta = sl*dkap - 2*pi;
% derivative of the opening angle
theta = fin(stheta,ngrid);
% integrate to find the angle
theta = theta + t0 + 2*pi*(0:ngrid)/ngrid;
% add the appropriate constant and a linear function to make it the
% acutal opening angle
theta(1) = t0;

end
