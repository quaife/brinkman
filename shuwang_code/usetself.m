% this code forms the normal and tan velocity given the angle and
% arclength and returns the normal and tan vel. from the ave. vel.
% integral
function [un,ut,rlambdalnew,xo,yo,forc1,forc2,xcc,ycc,...
    area] = usetself(x0,y0,sl,theta,rcon)

global kmatrix velocity bendsti bendratio uinside uoutside m
                 
% compute the curvature
dkap = acurv(sl,theta,m);

[xo,yo] = recon(m,x0,y0,sl,theta);     
%right hand side
sum2 = sum(sin(theta(1,1:m)).*xo(1,1:m) - ...
           cos(theta(1,1:m)).*yo(1,1:m))/2*sl/m;
sumxx = sum(sin(theta(1,1:m)).*(xo(1,1:m).^2))/2*sl/m;
sumyy = sum(-cos(theta(1,1:m)).*(yo(1,1:m).^2))/2*sl/m;
xcc = sumxx/sum2;
ycc = sumyy/sum2;
area = sum2;

% first stokes equation with unconstrained force
ulam = uinside/uoutside;

[forc,fforc] = uset_un(sl,theta,rcon,bendsti,bendratio,m);

% turn the force densty at tangential normal direction into x,y
% coordinate
forc1(1,1:m) = -(forc(1,1:m).*sin(theta(1,1:m)) + ...
                fforc(1,1:m).*cos(theta(1,1:m)));
forc2(1,1:m) = -(-forc(1,1:m).*cos(theta(1,1:m)) + ...
                fforc(1,1:m).*sin(theta(1,1:m)));
% first stokes equation with unconstrained force
% sigma1 and sigma2 are the x y velocity of the original vesicle

tau = [forc1 forc2]' ;
so.x = xo(1,1:m)'+ 1i*yo(1,1:m)'; 

A3 = selfmatrix(so,kmatrix); 

A = 2*A3;

c = sl/pi/m/(1+ulam)/uoutside/2;
c1 = sl/(1+ulam)/uoutside;

forc001_l = integral3(forc1,m);
forc002_l = integral3(forc2,m);

k = A*tau;

sigma1(1,1:m) = k(1:m)'*c - forc001_l(1:m)*c1 + ...
    yo(1,1:m)*velocity;
sigma2(1,1:m) = k(m+1:2*m)'*c - forc002_l(1:m)*c1;

% Now get the right hand side of the linear system which is the local
% incompressiblity condition 
uun(1,1:m) = sigma1(1,1:m).*sin(theta(1,1:m)) - ...
             sigma2(1,1:m).*cos(theta(1,1:m));
utn(1,1:m) = sigma1(1,1:m).*cos(theta(1,1:m)) + ...
             sigma2(1,1:m).*sin(theta(1,1:m));

utn(m+1) = utn(1);
utns = fd1(utn,m);

rhs(1,1:m) = -(utns(1,1:m)/sl + uun(1,1:m).*dkap(1,1:m));

% now solve the linear system
slam = stokes(dkap,m,sl,theta,rhs,A,c,c1);

%suppose we have a box, and the box has cordinates as x_box y_box.
slams = fd1(slam,m);

forcs1(1,1:m) = slam(1,1:m).*dkap(1,1:m).*sin(theta(1,1:m)) - ...
                slams(1,1:m).*cos(theta(1,1:m))/sl;
forcs2(1,1:m) = -slam(1,1:m).*dkap(1,1:m).*cos(theta(1,1:m)) - ...
                slams(1,1:m).*sin(theta(1,1:m))/sl;


forc1(1,1:m) = forc1(1,1:m)+forcs1(1,1:m);  
forc2(1,1:m) = forc2(1,1:m)+forcs2(1,1:m);  

tau=[forc1 forc2]' ;

k=A*tau;

forc001_l = integral3(forc1,m);
forc002_l = integral3(forc2,m);

un(1,1:m) = k(1:m)'*c - forc001_l(1:m)*c1 + yo(1,1:m)*velocity;
ut(1,1:m) = k(m+1:2*m)'*c - forc002_l(1:m)*c1;

rlambdalnew = slam;

end
