classdef curve
% This class implements that basic calculus on the curve.
% The basic data structure is a matrix X in which the columns 
% represent periodic C^{\infty} closed curves with N points, 
% X(1:n,j) is the x-coordinate of the j_th curve and X(n+1:N,j) 
% is the y-coordinate of the j_th curve; here n=N/2
% X coordinates do not have to be periodic, but the curvature,
% normals, etc that they compute will be garbage.

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y] = getXY(o,X)
% [x,y] = getXY(X) gets the [x,y] component of curves X
N = size(X,1)/2;
x = X(1:N,:);
y = X(N+1:end,:);

end % getXY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function V = setXY(o,x,y)
% V = setXY(x,y) sets the [x,y] component of vector V on the curve
N = size(x,1);
V=zeros(2*N,size(x,2));
V(1:N,:) = x;
V(N+1:end,:) = y;

end % setXY


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Dx,Dy] = getDXY(o,X)
% [Dx,Dy]=getDXY(X), compute the derivatives of each component of X 
% these are the derivatives with respect the parameterization 
% not arclength
x = X(1:end/2,:);
y = X(end/2+1:end,:);
N = size(x,1);
nv = size(x,2);
IK = o.modes(N);
Dx = o.diffFT(x,IK);
Dy = o.diffFT(y,IK);

end % getDXY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [jacobian,tangent,curvature] = diffProp(o,X)
% [jacobian,tangent,curvature] = diffProp(X) returns differential
% properties of the curve each column of the matrix X. Each column of 
% X should be a closed curve defined in plane. The tangent is the 
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
  IK = o.modes(N);
  [DDx,DDy] = oc.getDXY([Dx,Dy])
  curvature = (Dx.*DDy - Dy.*DDx)./(jacobian.^3);
end
% curvature of the curve

end % diffProp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
length = sum(sqrt(Dx.^2 + Dy.^2))*2*pi/N;
area = sum(x.*Dy - y.*Dx)*pi/N;

reducedArea = 4*pi*area./length.^2;

end % geomProp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

options = varargin;

if(any(strcmp(options,'angle')))
  theta = options{find(strcmp(options,'angle'))+1};
else
  theta = 0;
end
% rotate the vesicle by a certain angle

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

if(any(strcmp(options,'parameter')));
  alpha = options{find(strcmp(options,'parameter'))+1};
else
  alpha = (0:N-1)'*2*pi/N;
end
% Discretization in parameter space

if any(strcmp(options,'curly'))
  a = 1; b = 3*a; c = 0.85; 
  r = 0.5*sqrt( (a*cos(alpha)).^2 + (b*sin(alpha)).^2) + ...
      .07*cos(12*alpha);
  x = scale*c*r.*cos(alpha);
  y = scale*r.*sin(alpha);
  X0 = [x;y];
  % radius of curly vesicle

elseif any(strcmp(options,'star'))
  radius = 1 + 0.2*cos(folds*alpha);
  X0 = [radius.*cos(alpha);radius.*sin(alpha)];
  % a star that comes very close to intersecting itself
  % at the origin

elseif any(strcmp(options,'ellipse'))
  X0 = [cos(alpha);shortax*sin(alpha)];

else
  X0 = o.ellipse(N,ra);
  % build a vesicle of reduced area ra with N points
end 
% end of building reference vesicles.  Only need to rotate
% and shift them as desired

% X0 has a reference shape that needs to be scaled and rotated
X = zeros(2*N,1);
X(1:N) = scale*(cos(theta) * X0(1:N) - sin(theta) * X0(N+1:2*N));
X(N+1:2*N) = scale*(sin(theta) * X0(1:N) +  cos(theta) * X0(N+1:2*N));
% Rotate and scale vesicle

if ~equispaced
%  figure(1); clf
%  plot(X(1:end/2),X(end/2+1:end),'b-o')
%  axis equal
%  figure(2); clf;
%  plot(X(1:end/2),X(end/2+1:end),'b--')
%  axis equal
%  hold on
  sa = o.diffProp(X); % compute arclength of the vesicle shape
  alpha = o.arc(sa);
  % find parameter values that give discretization points that are
  % nearly equispaced in arclength
  [alpha,X] = o.initConfig(N,'true',options{1},'angle',theta,...
       'reducedArea',ra,'shortax',shortax,'scale',scale,...
       'folds',folds,'parameter',alpha);
end
%figure(2); hold on;
%plot(X(1:end/2),X(end/2+1:end),'ro')
%axis equal;

end % initConfig

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alpha = arc(o,sa)
% alpha = arc(sa) computes the parameter values of alpha that can be
% used to distribute points evenly in arclength

N = numel(sa);
length = sum(sa)*2*pi/N; 
% total length of vesicle. This is spectrally accurate since it is the
% trapezoid rule applied to a periodic function
tol = 1e-13; % tolerance
IK = o.modes(N);

intsa = o.intFT(sa,IK);

f = (0:N-1)'*2*pi/N*(length/2/pi) - intsa;
% we are interested in fiding the N roots of f. First one is at 0
df = -sa;
% derivative of f so that we can apply Newton's method
fh = fft(f)/N;
dfh = fft(df)/N;
% Fourier coeffients of f and its derivative

IK = o.modes(N); 
% Fourier modes needed for doing fourier interpolation

alpha = zeros(N,1);

for j = 2:N
  alpha(j) = alpha(j-1);
  for k = 1:100 % apply Newton's method
    falpha = real(sum(fh.*exp(IK*alpha(j)))) + ...
        ((j-1)*length/N - alpha(j)*length/2/pi);
    % f at alpha
    dfalpha = real(sum(dfh.*exp(IK*alpha(j))));
    % derivative of f at alpha
    alpha(j) = alpha(j) - falpha/dfalpha;
    % apply Newton iteration
    if abs(abs(falpha) < tol)
      break
    end
  end
end


end % arc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X0 = ellipse(o,N,ra)
% X0 = o.ellipse(N,ra) finds the ellipse (a*cos(theta),sin(theta)) so
% that the reduced area is ra.  Uses N points.  Parameter a is found 
% by using bisection method

t = (0:N-1)'*2*pi/N;
a = (1 - sqrt(1-ra^2))/ra;
% initial guess using approximation length = sqrt(2)*pi*sqrt(a^2+1)

X0 = [a*cos(t);sin(t)];

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rcon = initialConcentration(o,N,alpha,concentration,symmetry)
% rcon = initialConcentration(N,alpha,concentration,symmetry) sets up
% the initial lipid concentration defined on the vesicle.

if concentration == 0
  smallper = 0;
else
  smallper = 5e-2;
end

if symmetry == -1
  rcon = concentration + smallper*10*(rand(N,1) - 0.5);
  % uniform concentration with a little random noise
elseif symmetry == 0
  rcon = concentration + 3*smallper*cos(alpha) + ...
    0.5*smallper*cos(3*alpha) + ...
    0.5*smallper*cos(4*alpha);
  % uniform concentration with a little noise at specific even Fourier
  % modes
elseif symmetry==1
  rcon = concentration + 5*smallper*cos(2*alpha);
  % uniform concentration with one small even Fourier mode
elseif symmetry==2
  rcon = concentration+ 5*smallper*sin(alpha);
  % uniform concentration with one small odd Fourier mode
end

end % initialConcentration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [L,theta,cur] = computeOpeningAngle(o,N,X)
% [L,theta,cur] = computeOpeningAngle(N,X) finds the length, opening
% angle between the tangent vector and the x-axis, and the curvature. N
% is the number of discretization points and X = [x;y] are the
% discretization points.  This routine requires the points to be
% equispaced in arclength

[Dx,Dy] = o.getDXY(X);
% first derivative of the shape
if abs(Dx(1)) > 1e-12
  t0 = atan(Dy(1)/Dx(1));
else
  t0 = pi/2;
end
% initial opening angle
[DDx,DDy] = o.getDXY([Dx;Dy]);
% second derivative of the shape

L = sum(sqrt(Dx.^2 + Dy.^2))*2*pi/N;
% arclength which should be nearly constant since we constructed
% discretizations points that are equally spaced in arclength

cur = (Dx.*DDy - Dy.*DDx)/(L/2/pi)^3;

% compute curvature
IK = o.modes(N);
theta = L/2/pi*o.intFT(cur,IK);
% integrate the curature to find the opening angle

theta = theta + t0;
% add in the initial opening angle

end % computeOpeningAngle


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = recon(o,N,x0,y0,L,theta)
% X = recon(N,x0,y0,L,theta) computes the (x,y) coordinates of a
% geometry with length L, opening angle theta, and its first point at
% (x0,y0)

IK = o.modes(N);
x = x0 + o.intFT(L*cos(theta),IK)/2/pi;
y = y0 + o.intFT(L*sin(theta),IK)/2/pi;
X = o.setXY(x,y);

end % recon


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function df = diffFT(o,f,IK)
% df = diffFT(f,IK) Computes the first derivative of a periodic function
% f using fourier transform. IK is used to speed up the code.  It is the
% index of the fourier modes so that fft and ifft can be used

df = real(ifft(IK.*fft(f)));

end % diffFT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function intf = intFT(o,f,IK)
% intf = intFT(f,IK) computes the indefinite integral of a periodic
% function f using fourier transform

N = numel(f);
fh = fft(f);

zeroMode = -sum(fh(2:end)./IK(2:end));
% zero mode of the integral of f
intfh = [zeroMode;fh(2:end)./IK(2:end)];
% non-zero modes of the integral of f

intf = real(fh(1)/N * (0:N-1)'*2*pi/N + ifft(intfh));
% sum of the linear term and the contribution from the periodic part

end % intFT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IK = modes(o,N)
% IK = modes(N) builds the order of the fourier modes required for using
% fft and ifft to do spectral differentiation

IK = 1i*[(0:N/2-1) -N/2 (-N/2+1:-1)]';
% diagonal term for Fourier differentiation

end % modes

end % methods

end % classdef
