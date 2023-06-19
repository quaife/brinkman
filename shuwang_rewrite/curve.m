classdef curve
% This class implements that basic calculus on the curve.
% The basic data structure is a matrix X in which the columns 
% represent periodic C^{\infty} closed curves with N points, 
% X(1:n,j) is the x-coordinate of the j_th curve and X(n+1:N,j) 
% is the y-coordinate of the j_th curve; here n=N/2
% X coordinates do not have to be periodic, but the curvature,
% normals, etc that they compute will be garbage.

properties
IK;
end

methods

function o = curve(N)
    o.IK = o.modes(N);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y] = getXY(o,X)
% [x,y] = getXY(X) gets the [x,y] component of curves X
N = size(X,1)/2;
x = X(1:N,:);
y = X(N+1:end,:);

end % getXY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function V = setXY(o,x,y)
% V = setXY(x,y) sets the [x,y] component of vector V on the curve

N = size(x,1);
V=zeros(2*N,size(x,2));
V(1:N,:) = x;
V(N+1:end,:) = y;

end % setXY


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Dx,Dy] = getDXY(o,X)
% [Dx,Dy]=getDXY(X), compute the derivatives of each component of X 
% these are the derivatives with respect the parameterization 
% NOT arclength. Divide by length to find the derivative with respect to 
% arclength

x = X(1:end/2,:);
y = X(end/2+1:end,:);
N = size(x,1);
nv = size(x,2);
Dx = o.diffFT(x);
Dy = o.diffFT(y);

end % getDXY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [jacobian,tangent,curvature] = diffProp(o,X)
% [jacobian,tangent,curvature] = diffProp(X) returns differential
% properties of the curve each column of the matrix X. Each column of 
% X uld be a closed curve defined in plane. The tangent is the 
% normalized tangent.
%
% EXAMPLE:
%    n = 128; nv = 3;
%    X = boundary(n,'nv',nv,'curly');
%    oc = curve;
%    [k t s] = oc.diffProp(X);

N = size(X,1)/2;
nv = size(X,2);

% get the x y components
[Dx,Dy] = o.getDXY(X);

jacobian = sqrt(Dx.^2 + Dy.^2); 

if nargout>1  % if user requires tangent
  tangent = o.setXY( Dx./jacobian, Dy./jacobian);
end

if nargout>2  % if user requires curvature
  %IK = o.modes(N);
  [DDx,DDy] = o.getDXY([Dx;Dy]);
  DDx = DDx(:); DDy = DDy(:);
  curvature = (Dx.*DDy - Dy.*DDx)./(jacobian.^3);
end
% curvature of the curve

end % diffProp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [reducedArea,area,length] = geomProp(o,X)
% [reducedArea area length] = geomProp(X) calculate the length, area 
% and the reduced volume of domains inclose by columns of X. 
% Reduced volume is defined as 4*pi*A/L. 
% EXAMPLE:
%   X = boundary(64,'nv',3,'curly');
%   c = curve(X);
%   [rv A L] = c.geomProp(X);

[x,y] = o.getXY(X);
N = size(x,1);
[Dx,Dy] = o.getDXY(X);
length = sum(sqrt(Dx.^2 + Dy.^2))/N;
area = 0.5*sum(x.*Dy - y.*Dx)/N;

reducedArea = 4*pi*area./length.^2;

end % geomProp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [alpha,X] = initConfig(o,N,equispaced,varargin)       
% [alpha,X] = initConfig(n,varargin) returns N coordinates of boundary
% points. alpha are the values in [0,2*pi] that give points that are
% equispacced given a canonical parameterization, and X are the
% coordinates in the format [x;y].
% Available OPTIONS are
%
%   'curly'       - returns a somewhat wiggly vesicle.
%   'ellipse'     - elliptical shape
%   if neither of the above three are requested, it will create a
%   vesicle with the requsted reduced area, which is 0.65 if one
%   particarular value is not requested
%   'angle'       - inclination angle of the vesicle(s) form the
%                   vertical position, 
%   'scale'       - multiplies the size of the output boundary
%   'choke'       - returns a choked domain.
%   'longchoke'   - returns a long coked domain
%   'doublechoke' - returns a domain with two chokes
%   'contracting' - returns a domain for a couette apparatus.
%   'tube'        - returns a domain that is an ellongated ellipse


options = varargin;
if(any(strcmp(options,'angle')))
  theta = options{find(strcmp(options,'angle'))+1};
else
  theta = 0;
end
% rotate the vesicle by a certain angle

if(any(strcmp(options,'center')))
  cen = options{find(strcmp(options,'center'))+1};
else
  cen = [(0:0);zeros(1,1)];
end
% pick the center of the vesicles

if(any(strcmp(options,'reducedArea')))
  ra = options{find(strcmp(options,'reducedArea'))+1};
else 
  ra = 0.65; 
end
% desired reduced area

if(any(strcmp(options,'shortax')))
  shortax = options{find(strcmp(options,'shortax'))+1};
else
  shortax = 0.5;
end

if(any(strcmp(options,'scale')))
  scale = options{find(strcmp(options,'scale'))+1};
else
  scale = 1;
end
% scale the size of the vesicle (reduced area is invariant under
% scale)

if(any(strcmp(options,'folds')))
  folds = options{find(strcmp(options,'folds'))+1};
else
  folds = 1;
end
% number of folds for a star shape

if(any(strcmp(options,'parameter')))
  alpha = options{find(strcmp(options,'parameter'))+1};
else
  alpha = (0:N-1)'/N;
end
% Discretization in parameter space

if(any(strcmp(options,'geometry')))
  geometry = options{find(strcmp(options,'geometry'))+1};
else
  geometry = 'ellipse';
end
% Discretization in parameter space

if strcmp(geometry,'curly')
%if any(strcmp(options,'curly'))
  a = 1; b = 3*a; c = 0.85; 
  r = 0.5*sqrt( (a*cos(2*pi*alpha)).^2 + (b*sin(2*pi*alpha)).^2) + ...
      .07*cos(12*2*pi*alpha);
  x = scale*c*r.*cos(2*pi*alpha);
  y = scale*r.*sin(2*pi*alpha);
  X0 = [x;y];
  % radius of curly vesicle

elseif strcmp(geometry,'star')
%elseif any(strcmp(options,'star'))
  radius = 1 + 0.2*cos(folds*2*pi*alpha);
  X0 = [radius.*cos(2*pi*alpha);radius.*sin(2*pi*alpha)];
  % a star that comes very close to intersecting itself at the origin

elseif strcmp(geometry,'choke')
  a = 3; b = 1; c = 0.6; order = 8;
  % Parameterize t so that geometry is closer to equi-spaced in
  % arclength
  t = alpha;
  % parameters for the boundary
  r = (cos(2*pi*t).^order + sin(2*pi*t).^order).^(-1/order);
  x = a*r.*cos(2*pi*t); y = b*r.*sin(2*pi*t);
  ind = abs(x) < 1;
  y(ind) = y(ind).*(1-c*cos(pi*x(ind)))/(1+c);
  X0 = 3.5*[x;y];
  % choked domain.  a and b control the length and height.  c
  % controls the width of the gap, and order controls the
  % regularity

elseif strcmp(geometry,'longchoke')
  a = 30; b = 7/2; c = 0.67; order = 8;
  len = 25;
  % Parameterize t so that geometry is closer to equi-spaced in
  % arclength
  t = alpha;
  % parameters for the boundary
  r = (cos(2*pi*t).^order + sin(2*pi*t).^order).^(-1/order);
  x = a*r.*cos(2*pi*t); y = b*r.*sin(2*pi*t);
  ind = abs(x) < len;
  scalRat = 2*c/(1+c)*(0.5-0.5*cos(pi*x(ind(1:end/2))/len)).^10 + ...   
            (1-c)/(1+c);
  y(ind) = y(ind).*[scalRat;scalRat];
  X0 = [x;y];
  % choked domain.  a and b control the length and height.  c
  % controls the width of the gap, and order controls the
  % regularity
%   
elseif strcmp(geometry,'doublechoke')
%elseif any(strcmp(options,'doublechoke'))
  a = 10; b = 3; c = 0.6; order = 8;
  shift = pi/2 + 0.1;
  % parameters for the boundary
  t = (0:N-1)'*2*pi/N;
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*cos(t); y = b*r.*sin(t);
  ind = abs(x-shift) < pi/2;
  y(ind) = y(ind).*(1-c*cos(2*(x(ind)-shift)))/(1+c);
  ind = abs(x+shift) < pi/2;
  y(ind) = y(ind).*(1-c*cos(2*(x(ind)+shift)))/(1+c);
  X0 = [x;y];
  % choked domain.  a and b control the length and height.  c
  % controls the width of the gap, and order controls the
  % regularity.  shift controls the distance between the two
  % regions where the domain is restricted

elseif strcmp(geometry,'contracting')
%elseif any(strcmp(options,'contracting'))
  w = 0.5; % width of the opening
  ell1 = 3.0; % length before contracting region
  ell2 = 12; % length (in x) of contracting region
  ell3 = 3.0; % length after contracting region
  angd = 14; % angle in degrees of contracting region
  ang = angd*pi/180;

  W = 2*ell2*tan(ang) + w;  % total width of the geometry

  % Set up 10 marker points
  z1 = 0 - 1i*(W/2);
  z2 = z1 + ell1;
  z3 = z2 + (ell2 + 1i*ell2*tan(ang));
  z4 = z3 - 1i*ell2*tan(ang);
  z5 = z4 + ell3;
  z6 = z5 + 1i*W;
  z7 = z6 - ell3;
  z8 = z7 - 1i*ell2*tan(ang);
  z9 = z8 - (ell2 - 1i*ell2*tan(ang));
  z10 = z9 - ell1;

  dz = 0.1; % length spacing between successive nodes

  % setup roughly equispaced points along the geometry
  len1 = abs(z1 - z2);
  zz1 = linspace(z1,z2,round(len1/dz)); zz1 = zz1(1:end-1);
  len2 = abs(z2 - z3);
  zz2 = linspace(z2,z3,round(len2/dz)); zz2 = zz2(1:end-1);
  len3 = abs(z3 - z4);
  zz3 = linspace(z3,z4,round(len3/dz)); zz3 = zz3(1:end-1);
  len4 = abs(z4 - z5);
  zz4 = linspace(z4,z5,round(len4/dz)); zz4 = zz4(1:end-1);
  len5 = abs(z5 - z6);
  zz5 = linspace(z5,z6,round(len5/dz)); zz5 = zz5(1:end-1);
  len6 = abs(z6 - z7);
  zz6 = linspace(z6,z7,round(len6/dz)); zz6 = zz6(1:end-1);
  len7 = abs(z7 - z8);
  zz7 = linspace(z7,z8,round(len7/dz)); zz7 = zz7(1:end-1);
  len8 = abs(z8 - z9);
  zz8 = linspace(z8,z9,round(len8/dz)); zz8 = zz8(1:end-1);
  len9 = abs(z9 - z10);
  zz9 = linspace(z9,z10,round(len9/dz)); zz9 = zz9(1:end-1);
  len10 = abs(z10 - z1);
  zz10 = linspace(z10,z1,round(len10/dz)); zz10 = zz10(1:end-1);

  z = [zz1,zz2,zz3,zz4,zz5,zz6,zz7,zz8,zz9,zz10];
  winsize = 10; % Gaussian window size
  zPad = [z(end-winsize+1:end),z,z(1:winsize)];
  % pad so that periodicity is incorporated into smoothing
  z = smoothdata(zPad,'gaussian',winsize);
  % smooth padded data
  z = z(winsize+1:end - winsize);
  % remove the padding
  nz = numel(z);
  nzmodes = [(0:nz/2) (-nz/2+1:-1)];
  zft = fft(z)/nz;
  z = zeros(N,1);
  for k = 1:nz
     z = z + zft(k)*exp(2*pi*1i*nzmodes(k)*alpha); 
  end
  %z = interpft(z,N);
  % interpolate with FFT to the requested number of discretization
  % points
  X0 = [real(z);imag(z)];
  % store as (x,y) rather than x + 1i*y

elseif strcmp(geometry,'tube')
%elseif any(strcmp(options,'tube'))
  a = 3; b = 1; order = 8;
  % parameters for the boundary
%  Nsides = ceil(0.5*b/(2*a+2*b)*N);
%  Ntop = (N-4*Nsides)/2;
%  t1 = linspace(0,0.2*pi,Nsides+1); t1 = t1(1:end-1)';
%  t2 = linspace(0.2*pi,pi-0.2*pi,Ntop+1); t2 = t2(1:end-1)';
%  t3 = linspace(pi-0.2*pi,pi+0.2*pi,2*Nsides+1); t3 = t3(1:end-1)';
%  t4 = linspace(pi+0.2*pi,2*pi-0.2*pi,Ntop+1); t4 = t4(1:end-1)';
%  t5 = linspace(2*pi-0.2*pi,2*pi,Nsides+1); t5 = t5(1:end-1)';
%  t = [t1;t2;t3;t4;t5];
  t = alpha;
  % Parameterize t so that geometry is closer to 
  % equispaced in arclength
  r = (cos(2*pi*t).^order + sin(2*pi*t).^order).^(-1/order);
  x = a*r.*cos(2*pi*t); y = b*r.*sin(2*pi*t);
  X0 = 3.5*[x;y];
  % rounded off cylinder.  a and b control the length and height 
  % and order controls the regularity
else
  X0 = [cos(2*pi*alpha);shortax*sin(2*pi*alpha)]; 
  %X0 = o.ellipse(N,ra);
  % build a vesicle of reduced area ra with N points
end 
% end of building reference vesicles.  Only need to rotate
% and shift them as desired

% X0 has a reference shape that needs to be scaled and rotated
X = zeros(2*N,1);
X(1:N) = scale*(cos(theta) * X0(1:N) - sin(theta) * X0(N+1:2*N)) + cen(1,1);
X(N+1:2*N) = scale*(sin(theta) * X0(1:N) +  cos(theta) * X0(N+1:2*N)) + cen(2,1);
% Rotate and scale vesicle
if ~equispaced
  sa = o.diffProp(X); % compute arclength of the vesicle shape
  alpha = o.arc(sa);
  % find parameter values that give discretization points that are
  % nearly equispaced in arclength\
 
  [alpha,X] = o.initConfig(N,true,'angle',theta,...
       'reducedArea',ra,'shortax',shortax,'scale',scale,...
       'folds',folds,'parameter',alpha, 'center', cen,...
       'geometry',geometry);
  % [alpha, X] = o.initConfig(N,true,options,'parameter',alpha);
   
end
% plot(X(1:N),X(N+1:2*N))
% pause
end % initConfig

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alpha = arc(o,sa)
% alpha = arc(sa) computes the parameter values of alpha that can be
% used to distribute points evenly in arclength

N = numel(sa);
length = sum(sa)/N; 
% total length of vesicle. This is spectrally accurate since it is the
% trapezoid rule applied to a periodic function
tol = 1e-13; % tolerance
%IK = o.modes(N);

intsa = o.intFT(sa,o.IK);

f = (0:N-1)'*length/N - intsa;
% we are interested in fiding the N roots of f. First one is at 0
df = -sa;
% derivative of f so that we can apply Newton's method
fh = fft(f)/N;

dfh = fft(df)/N;
% Fourier coeffients of f and its derivative

alpha = zeros(N,1);

for j = 2:N
  alpha(j) = alpha(j-1);
  for k = 1:100 % apply Newton's method
    falpha = real(sum(fh.*exp(o.IK*alpha(j)))) + ...
        ((j-1)*length/N - alpha(j)*length);
    % f at alpha
    dfalpha = real(sum(dfh.*exp(o.IK*alpha(j))));
    % derivative of f at alpha
    alpha(j) = alpha(j) - falpha/dfalpha;
    % apply Newton iteration
    if abs(abs(falpha) < tol)
      break
    end
  end
end

end % arc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X0 = ellipse(o,N,ra)
% X0 = o.ellipse(N,ra) finds the ellipse (a*cos(theta),sin(theta)) so
% that the reduced area is ra.  Uses N points.  Parameter a is found 
% by using bisection method

t = (0:N-1)'*2*pi/N;
a = (1 - sqrt(1-ra^2))/ra;
% initial guess using approximation length = sqrt(2)*pi*sqrt(a^2+1)

X0 = [cos(t);a*sin(t)];

[raNew,~,~] = o.geomProp(X0);

cond = abs(raNew - ra)/ra < 1e-4;
maxiter = 10;
iter = 0;
while (~cond && iter < maxiter);
  iter = iter + 1;
  if (raNew > ra)
    a = 0.9 * a;
  else
    a = 1.05 * a;
  end
  % update the major axis
  X0 = [cos(t);a*sin(t)];
  % Compute new possible configuration
  [raNew,~,~] = o.geomProp(X0);
  % compute new reduced area
  cond = abs(raNew - ra) < 1e-2;
  % check if the residual is small enough
end
% iteration quits if reduced area is achieved within 1% or 
% maxiter iterations have been performed

end % ellipse


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rcon = initialConcentration(o,N,alpha,concentration,symmetry)
% rcon = initialConcentration(N,alpha,concentration,symmetry) sets up
% the initial lipid concentration defined on the vesicle.
if concentration == 0
  smallper = 0;
else
  smallper = 5e-2;
end
if symmetry == -1
  rng(10245621);
  rcon = smallper*10*(rand(N,1) - 0.5);
  rcon = rcon - mean(rcon) + concentration;
  % uniform concentration with a little random noise
elseif symmetry == 0
  rcon = concentration + 3*smallper*cos(2*pi*alpha) + ...
    smallper*0.5*cos(3*2*pi*alpha) + ...
    smallper*0.5*cos(4*2*pi*alpha);
  % uniform concentration with a little noise at specific even Fourier
  % modes
elseif symmetry==1
  rcon = concentration + 5*smallper*cos(2*2*pi*alpha);
  % uniform concentration with one small even Fourier mode
elseif symmetry==2
  rcon = concentration+ 5*smallper*sin(4*pi*alpha);
  % uniform concentration with one small odd Fourier mode
elseif symmetry == 3
  rcon = ones(N,1);
  rcon(1:ceil(concentration*N)) = 0;
elseif symmetry == 4
  rcon = ones(N,1);
  rcon(1:ceil(concentration*N/2)) = 0;
  rcon(end-ceil(concentration*N/2):end) = 0;
elseif symmetry == 5
  rcon = zeros(N,1);
  rcon(1:ceil(concentration*N/2)) = 1;
  rcon(end-ceil(concentration*N/2):end) = 1;
elseif symmetry == 6
  rcon = ones(N,1);  
  rcon(ceil(concentration*N):end) = 0;
elseif symmetry == 7
  rcon = ones(N,1);  
  rcon(1:ceil(concentration*N)) = 0;
end


end % initialConcentration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [L,theta,cur] = computeOpeningAngle(o,N,X)
% [L,theta,cur] = computeOpeningAngle(N,X) finds the length, opening
% angle between the tangent vector and the x-axis, and the curvature. N
% is the number of discretization points and X = [x;y] are the
% discretization points.  This routine requires the points to be
% equispaced in arclength

[Dx,Dy] = o.getDXY(X);
% first derivative of the shape
%if abs(Dx(1)) > 1e-8
%  t0 = atan(Dy(1)/Dx(1));
%else
%  t0 = pi/2;
%end
t0 = angle(Dx(1) + 1i*Dy(1));
% initial opening angle
[DDx,DDy] = o.getDXY([Dx;Dy]);

% second derivative of the shape;
L = sum(sqrt(Dx.^2 + Dy.^2))/N;

% arclength which should be nearly constant since we constructed
% discretizations points that are equally spaced in arclength

cur = (Dx.*DDy - Dy.*DDx)/L^3;
% compute curvature

theta = o.intFT(L*cur-2*pi,o.IK);
% integrate the curature to find the opening angle

theta = theta + t0 + (0:N-1)'*2*pi/N;
% add in the initial opening angle

end % computeOpeningAngle


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = recon(o,N,x0,y0,L,theta)
% X = recon(N,x0,y0,L,theta) computes the (x,y) coordinates of a
% geometry with length L, opening angle theta, and its first point at
% (x0,y0)

% IK = o.modes(N);
x = o.intFT(L*cos(theta),o.IK);
y = o.intFT(L*sin(theta),o.IK);
avx = sum(L*cos(theta))/N;
avy = sum(L*sin(theta))/N;
x = x - (0:N-1)'*avx/N + x0;
y = y - (0:N-1)'*avy/N + y0;

X = o.setXY(x,y);

end % recon


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function df = diffFT(o,f)
% df = diffFT(f,IK) Computes the first derivative of a periodic function
% f using fourier transform. IK is used to speed up the code.  It is the
% index of the fourier modes so that fft and ifft can be used
df = real(ifft(o.IK.*fft(f)));

end % diffFT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function intf = intFT(o,f,IK)
% intf = intFT(f,IK) computes the indefinite integral of a periodic
% function f using fourier transform

N = numel(f);
fh = fft(f);

zeroMode = -sum(fh(2:end)./IK(2:end));
% zero mode of the integral of f
intfh = [zeroMode;fh(2:end)./IK(2:end)];
% non-zero modes of the integral of f

intf = real(fh(1)/N * (0:N-1)'/N + ifft(intfh));
% sum of the linear term and the contribution from the periodic part

end % intFT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IK = modes(o,N)
% IK = modes(N) builds the order of the fourier modes required for using
% fft and ifft to do spectral differentiation
IK = 2*pi*1i*[(0:N/2-1) -N/2 (-N/2+1:-1)]';
% diagonal term for Fourier differentiation

end % modes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dcur = fdcur(o,ves,unt)
% compute fourier derivative of opening tangent angle using "trick" by
% subtracting a function that goes from 0 to 2*pi over the range of
% alpha values
% This routine is equation (61) in the Sohn et al 2010 JCP paper

N = ves.N;
%IK = o.modes(N);

dthetads = o.rmLinearPart(N, ves.theta);
dcur = intFT(o,[dthetads.*unt],o.IK);

end %fdcur

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dx = rmLinearPart(o,N,theta)
% N = ves.N;
% theta = ves.theta;
%IK = o.modes(N);

% function with the same increase over a period of 2*pi as the angle
% theta

% need to filter theta a bit
IK = o.IK;
%o.IK(abs(o.IK) > 2*pi*N/4) = 0;

dx = o.diffFT([theta - 2*pi*(0:1:N-1)'/N]) + 2*pi;
%semilogy(abs(fftshift(fft(theta - 2*pi*(0:1:N-1)'/N))))
%semilogy(abs(fftshift(fft(dx))))
%pause

o.IK = IK;

end %rmLinearPart

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fntheta = fthetaim(o,ves,un,ut,fdotn)

N = ves.N;
L = ves.L;  
theta = ves.theta;
%IK = o.modes(N);
% compute right hand side of the ODE for \theta in equation (8) in Sohn
% et al JCP 2010
ftheta = o.forceTheta(ves,un,ut,fdotn);
% now we must subtract off the stiffest part in Fourier space.
% To the power of 3 term comes from the third-derivative of the function
% that L is being applied to in equation (53)
rlen = -ves.bendsti*(abs(o.IK)/L).^3/4 - ves.SPc*ves.bendsti*(o.IK/L).^4; 
% last term in eq(54)
var1 = fft(ftheta); % fft of first 2 terms in (54)
var2 = fft(theta-(2*pi*(0:N-1))'/N);
fntheta = real(ifft(var1 - rlen.*var2));
%Krasney filter applied to smooth out spurious high frequency terms
fntheta = o.kfilter(fntheta,N);

end %fthetaim

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ftheta = forceTheta(o,ves,un,ut,fdotn)
%forcing term for the ODE for theta 
N = ves.N;
L = ves.L;
%IK = o.modes(N);
%compute first derivative of the opening angle using Fourier
%differentiation. This is actually the curvature
dthetads = o.rmLinearPart(N,ves.theta)/L;
%compute the derivative of the normal velocity
dunds = o.diffFT(un + ves.SPc*fdotn)/L;    
% Check T_s = -V * curvature (local inextensibility condition) (see
% equation (20) in Sohn et al JCP 2010)
ftheta = -dunds + dthetads.*ut;
    
end % forceTheta

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filtered = kfilter(o,ftheta,N)
%fourier filter with krasny filtering with 1 input vector 
%USING FAST FOURIER TRANSFORM
tol = 1e-8;
b = fft(ftheta);
b(abs(b)/N < tol) = 0;
b(N/2:N/2+2) = 0;
c = real(b.*exp(-10*([0:2:N N-2:-2:2]'/N).^25));
d = imag(b.*exp(-10*([1:2:N-1 0 N-1:-2:3]'/N).^25));
d(N/2-1) = 0; d(N/2+3) = 0;
filtered = real(ifft(c+1i*d));

end %kfilter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function N2 = frconim(o,ves,epsch,consta)
%This routine is about solving equation (75) for N2Hat.
N = ves.N;
L = ves.L;
rcon = ves.rcon;
%IK = o.modes(N);
%setting up RHS of eq (68) 
%termd is the a/eps(f'(u)-eps^2u_ss) + b'(u)/2 * kappa^2 terms in eq (13) 
%ie. termd is the variational derivative of the energy with
%respect to the lipid concentration
termd = o.fluxj(ves,epsch,consta);
%now taking derivative of the variational derivative twice
rcons = o.diffFT(termd)/L;
rconss = o.diffFT(rcons)/L;
%invPe is the 1/Pe term in eq (67)
invPe = 1;
%fcon is the part of N_2 in eq (67) but is still missing the
fcon = invPe*rconss;
fconHat = fft(fcon);
rconHat = fft(rcon);
%rlen is the term that needs to multiply the fourier coefficients of
%the lipid species to result in the fourth derivative scaled by the a
%and \eps constants in equation (67)
rlen = invPe*(o.IK/L).^4*consta*epsch;
%N2 is equation (67) in the physical space.
N2Hat = fconHat + rlen.*rconHat;
N2Hat(1) = 0;
N2 = real(ifft(N2Hat));

end %frconim

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RHSflux = fluxj(o,ves,epsch,consta)
%This routine computes the forcing on the RHS flux OF THE Cahn_Hillard
%Equation for rcon periodic on [0,1] and for closed curves. RHSflux are the
%a/eps(f'(u)-eps^2u_ss) + b'(u)/2 * kappa^2 terms in eq (13) 

N = ves.N;
L = ves.L;
rcon = ves.rcon;
%IK = o.modes(N);
%dfu is the derivative of the double-well potential
dfu = 0.5*(rcon.*(ones(N,1)-rcon).^2 - rcon.^2.*(ones(N,1)-rcon));
%derivative of the concentration 
rcons = o.diffFT(rcon)/L;
%second derivative of the concentration
rconss = o.diffFT(rcons)/L;
%term2 is eps^2 * u_ss as defined in equation (13)
term2 = -epsch^2*rconss;
%computing the b'(u)/2 * kappa^2 term in eq (13)
b0 = ves.bendsti;
b1 = ves.bendsti*ves.bendratio;
%rbn is b(u).
%rbn = b0*(ones(N,1)-rcon) + b1*rcon;
rbn = (b1 - b0)/2*tanh(3*(rcon - 0.5)) + (b0 + b1)/2;
%rbndu is b'(u) 
%rbndu = (b1 - b0)*ones(N,1);
rbndu = (3/2)*(b1 - b0)*sech(3*(rcon-0.5)).^2; 
%compute the curvature
dkap = o.acurv(N,ves.theta,L);
dkap = ves.cur; % TODO: COMMENT OUT LATER
%term3 is b'(u)/2 * kappa^2 in eq (13)
term3 = 0.5*rbndu.*dkap.^2;
%RHSflux is a/eps(f'(u)-eps^2u_ss) + b'(u)/2 * kappa^2 terms in eq (13) 
RHSflux = consta/epsch * (dfu + term2) + term3;

end % fluxj

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function force = CHforce(o,ves,epsch,consta)
% This routine computes the forcing on the RHS flux OF THE Cahn_Hillard,
% but only the phase contribution. In particular, it computes
% a/eps(f'(u)-eps^2u_ss) as defined in eq (13) 

N = ves.N;
L = ves.L;
rcon = ves.rcon;
% dfu is the derivative of the double-well potential
fu = 0.25*rcon.^2.*(1-rcon).^2;
%dfu = 0.5*(rcon.*(ones(N,1)-rcon).^2 - rcon.^2.*(ones(N,1)-rcon));
% derivative of the concentration 
rcons = o.diffFT(rcon)/L;
% second derivative of the concentration
%rconss = o.diffFT(rcons)/L;
% phase_gradient is eps^2 * u_ss as defined in equation (13)
phase_gradient = -epsch^2/2*rcons.^2;
%phase_gradient = -epsch^2*rconss;
%force = consta/epsch * (dfu + phase_gradient);
force = consta/epsch * (fu + phase_gradient);


end % CHforce

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dkap = acurv(o,N,theta,L)
% This routine computes the curvature from the tangent angle formulation

dkap = o.rmLinearPart(N,theta)/L;

end %acurv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CoM = centerOfMass(o,X)
% This function returns the center of mass of the vesicle
    
  N = length(X)/2;  
  [~,Area,L] = o.geomProp(X);
  [~,tangent,~] = o.diffProp(X);
  normal = [tangent(N+1:end);-tangent(1:N)];
%   plot(X(1:N), X(N+1:end))
%   hold on
%   quiver(X(1:N), X(N+1:end),normal(1:N),normal(N+1:end))
%   axis equal
%   pause
  CoM = [0.5*L*sum(X(1:N).^2.*normal(1:N))/(N*Area),0.5*L*...
        sum(X(N+1:end).^2.*normal(N+1:end))/(N*Area)];

end %centerOfMass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % methods

end % classdef
