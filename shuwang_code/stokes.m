function  vel=stokes(dkap,m,sl,theta,rhs,A,c,c1)
 
tol = 1e-10;  maxit = numel(rhs); 
rhscolume=rhs';
afun(rhscolume);
% based on the variable name, it seems that this is solving for a
% velocity
[vell,flag,relres,iter,resvec] = gmres(@afun,rhscolume,m/2,tol,maxit);
%,@mfun);
vel = vell';
%  flag
%  relres
%  iter
%  %resvec
%  pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vell = afun(pss)
% matvec corresponding to equation (40)

slam = pss';

slams = fd1(slam,m);

forc1 = +slam(1,1:m).*dkap(1,1:m).*sin(theta(1,1:m)) - ...
         slams(1,1:m).*cos(theta(1,1:m))/sl;
forc2 = -slam(1,1:m).*dkap(1,1:m).*cos(theta(1,1:m)) - ...
         slams(1,1:m).*sin(theta(1,1:m))/sl;  

k = A*[forc1 forc2]';
% norm(k)
% pause
forc001_l = integral3(forc1,m);
forc002_l = integral3(forc2,m);

sigma1(1,1:m)=(k(1:m)     )'*c - forc001_l(1:m)*c1;
sigma2(1,1:m)=(k(m+1:2*m) )'*c - forc002_l(1:m)*c1;

uun(1,1:m) = sigma1(1,1:m).*sin(theta(1,1:m)) - ...
             sigma2(1,1:m).*cos(theta(1,1:m));
utn(1,1:m) = sigma1(1,1:m).*cos(theta(1,1:m)) + ...
             sigma2(1,1:m).*sin(theta(1,1:m));
utn(m+1) = utn(1);
utns = fd1(utn,m);

vel(1,1:m) = utns(1,1:m)/sl + uun(1,1:m).*dkap(1,1:m);
vell=vel';
% seems that the output is (u \cdot s)_s + kappa * (u \cdot n)

end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preconditioner (equation (63) and (64)) for equation (40) which is the
% Schur complement for a Lagrange multiplier needed to enforce local
% inextensibility
function zz = mfun(rr)

r = rr';
ssigma(1,1:m) = r(1,1:m);
ssigma(1,m+1) = ssigma(1,1);

ssigma = fft(ssigma,m);
coe = [1 1:1:m/2 m/2-1:-1:1];
ssigma(1,1:m) = ssigma(1,1:m)./coe;
ssigma(1,1:m) = real(ifft(ssigma(1,1:m),m));
z(1,1:m) = ssigma(1,1:m);
zz = z';
end

end
 
