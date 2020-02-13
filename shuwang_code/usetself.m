% this code forms the normal and tan velocity given the angle and
% arclength and returns the normal and tan vel. from the ave. vel.
% integral
function [un,ut,rlambdalnew,xo,yo,forc1,forc2,xcc,ycc,...
    area] = usetself(x0,y0,sl,theta,rcon)
% INPUT VARIABLES
% x0 - x-coordinate of the single reference point
% y0 - y-coordinate of the single reference point
% sl - length of the interface
% theta - angle between the tangent vector and horizontal
% rcon - lipid species concentration

% OUTPUT VARIABLES
% un - normal velocity
% ut - tangential velocity
% rlambdalnew - Lagrange multiplier for inextensibility which is given
% as T_s + \kappa V = 0
% xo - x-coordinates coming from the length, angle, and reference point
% yo - y-coordinates coming from the length, angle, and reference point
% forc1 - 
% forc2 - 
% xcc - 
% ycc - 
% area - area inside the interface



global kmatrix velocity bendsti bendratio uinside uoutside m
                 
% compute the initial curvature using the tangent angle
dkap = acurv(sl,theta,m);
% Use the arclength dkap (calculated above), tan angle theta (input)
% to recunstruct the curve to find the initial x and y
[xo,yo] = recon(m,x0,y0,sl,theta);     

%Green's identity on (1/2)Integral(x dot n) to find the area inside 
%the vesicle (Integral about the surface of dx)
sum2 = sum(sin(theta(1,1:m)).*xo(1,1:m) - ...
           cos(theta(1,1:m)).*yo(1,1:m))/2*sl/m;
%Not sure what sumxx and sumyy are yet
sumxx = sum(sin(theta(1,1:m)).*(xo(1,1:m).^2))/2*sl/m;
sumyy = sum(-cos(theta(1,1:m)).*(yo(1,1:m).^2))/2*sl/m;

%Some sort of area ratio?
xcc = sumxx/sum2;
ycc = sumyy/sum2;

%This comes from Green's identity Integral dx about the surface = area
area = sum2;

% viscosity contrast
ulam = uinside/uoutside;

%find the normal and tangent velocities from the average velocity
%integral possibly - unconstrained from bending
%I think calling uset_un has something to do with calculating the 
%variational derivatives in eq (13) & (14) unconstrained from bending in 
%the RHS of eqn(33)
[forc,fforc] = uset_un(sl,theta,rcon,bendsti,bendratio,m);
%now that we have the variational derivatives, we compute the rest of the
%jump condition in eq (33) where n = (sin,-cos), s = (cos, sin)
forc1(1,1:m) = -(forc(1,1:m).*sin(theta(1,1:m)) + ...
                fforc(1,1:m).*cos(theta(1,1:m)));
forc2(1,1:m) = -(-forc(1,1:m).*cos(theta(1,1:m)) + ...
                fforc(1,1:m).*sin(theta(1,1:m)));
%Put forc1 and forc2 into a single array
tau = [forc1 forc2]' ;
%Define new variable so s.t. so'so^T = (x')^2+(y')^2 = ||r'||^2
so.x = xo(1,1:m)' + 1i*yo(1,1:m)'; 
%Construct Stokes matrix
A3 = selfmatrix(so,kmatrix); 
%Where does this 2 come from??
A = 2*A3;
%k is the velocity on the interface corresponding to v^u in eq (32)
k = A*tau; %[[P^u n]]_sigma = tau, then k = stokesMatrix*tau

% -----------------------------------------------------------------
c = sl/pi/m/(1+ulam)/uoutside/2;
c1 = sl/(1+ulam)/uoutside;

forc001_l = integral3(forc1,m);
forc002_l = integral3(forc2,m);

sigma1(1,1:m) = k(1:m)'*c - forc001_l(1:m)*c1 + ...
    yo(1,1:m)*velocity;
sigma2(1,1:m) = k(m+1:2*m)'*c - forc002_l(1:m)*c1;
% -----------------------------------------------------------------

%Calculate v dot n in eq (40)
uun(1,1:m) = sigma1(1,1:m).*sin(theta(1,1:m)) - ...
             sigma2(1,1:m).*cos(theta(1,1:m));
%Calculate v dot s in eq (40)
utn(1,1:m) = sigma1(1,1:m).*cos(theta(1,1:m)) + ...
             sigma2(1,1:m).*sin(theta(1,1:m));
utn(m+1) = utn(1);
%Calculate d/ds(v dot s) in eq (40)
utns = fd1(utn,m);

%Calculate the right hand side of equation (40)
rhs(1,1:m) = -(utns(1,1:m)/sl + uun(1,1:m).*dkap(1,1:m));

% The velocity components in eq (40) are nonlocal linear functions of 
% lambdaTilde. 
% Call stokes which solves the linear system for LambdaTilde in 39 using 
% GMRES. Each iteration of GMRES requires a solution of stokes equation 
slam = stokes(dkap,m,sl,theta,rhs,A,c,c1); %slam is LambdaTilde in eq (39)?
%calculate the fourier derivative
slams = fd1(slam,m);
%We can now calculate the traction jump in first part of equation (39)...?
forcs1(1,1:m) = slam(1,1:m).*dkap(1,1:m).*sin(theta(1,1:m)) - ...
                slams(1,1:m).*cos(theta(1,1:m))/sl;
forcs2(1,1:m) = -slam(1,1:m).*dkap(1,1:m).*cos(theta(1,1:m)) - ...
                slams(1,1:m).*sin(theta(1,1:m))/sl;

%Are we adding the jump conditions in eq (39) and (33) here?
forc1(1,1:m) = forc1(1,1:m) + forcs1(1,1:m);  
forc2(1,1:m) = forc2(1,1:m) + forcs2(1,1:m);  


tau=[forc1 forc2]' ;
%Calculating u tilde in eq (40)?
k=A*tau;

%.....
forc001_l = integral3(forc1,m);
forc002_l = integral3(forc2,m);

%Are we updating u in eq (31) now with the u tilde?
% compute normal and tangential velocities
un(1,1:m) = k(1:m)'*c - forc001_l(1:m)*c1 + yo(1,1:m)*velocity;
ut(1,1:m) = k(m+1:2*m)'*c - forc002_l(1:m)*c1;

rlambdalnew = slam;

end
