function [x,y,theta,rcon,sl] = initiall(...
    shortax,ngrid,concentration,symmetry,smallper)
% x - x coordinate of vesicle
% y - y coordinate of vesicle
% theta - opening angle of vesicle
% rcon - lipid species concentration
% sl - length of vesicle

% shortax - short axis
% ngrid - number of grid points
% concentration - concentration value?
% symmetry - concentration profile?
% smallper - small parameter for introducing parameters in the lipid
% species concentration

tol = 1e-12;
N = (0:1:ngrid);
x0 = cos(2*pi*N/ngrid);
y0 = shortax*sin(2*pi*N/ngrid);
% initial distribution of the points, but not equispaced in arclength

dx = fd1(x0,ngrid);
dy = fd1(y0,ngrid);
% derivatives of the initial vesicle shape functions

slo = sqrt(dx.^2+dy.^2);
slo(ngrid+1) = slo(1);
% arclength
alpha = arc(slo,ngrid);
% find parameter values in [0 1] so that the parameterization gives
% points equispaced in arclength

x = cos(pi*2*alpha);
y = shortax*sin(pi*2*alpha);
% x and y coordinates of the vesicle, but equispaced in arclength
sxn = fd1(x,ngrid);
syn = fd1(y,ngrid);
% compute derivatives of the vesicle shape
% sxn.^2 + syn.^2 should be nearly constant
sl = sum(sqrt(sxn(1:ngrid).^2+syn(1:ngrid).^2));
sl = sl/ngrid;
% compute total length of the vesicle
if abs(sxn(1)) > tol
  t0 = atan(abs(syn(1)/sxn(1)));
else
  t0 = pi/2;
end
% find intitial opening angle of tangent vector
theta = thetsolve(sxn,syn,sl,t0,ngrid);
% compute the opening angle between the tangent vector and the
% horizontal axis
if symmetry == -1
  rcon = concentration + smallper*10*(rand(1,ngrid+1) - 0.5);
  % uniform concentration with a little random noise
elseif symmetry==0
  rcon = concentration + 3*smallper*cos(alpha*pi*2) + ...
    0.5*smallper*cos(3*alpha*pi*2) + ...
    0.5*smallper*cos(4*alpha*pi*2);
  % uniform concentration with a few small Fourier modes
elseif symmetry==1
  rcon = concentration + 5*smallper*cos(2*alpha);
  % uniform concentration with one small even Fourier mode
elseif symmetry==2
  rcon = concentration+ 5*smallper*sin(alpha*pi*2); 
  % uniform concentration with one small odd Fourier mode
end

end
