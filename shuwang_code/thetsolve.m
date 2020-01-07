function theta = thetsolve(sxn,syn,sl,t0,ngrid)

ssxn = fd1(sxn,ngrid);
ssyn = fd1(syn,ngrid);

dkap = (sxn.*ssyn-syn.*ssxn)/(sl^3);
dkap(ngrid+1) = dkap(1);
stheta = sl*dkap-2*pi;
theta = fin(stheta,ngrid);
theta = theta+t0+pi*2*(0:ngrid)/ngrid;
theta(1) = t0;

end
