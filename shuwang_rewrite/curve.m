classdef curve
% This class implements that basic calculus on the curve.
% The basic data structure is a matrix X in which the columns 
% represent periodic C^{\infty} closed curves with N points, 
% X(1:n,j) is the x-coordinate of the j_th curve and X(n+1:N,j) 
% is the y-coordinate of the j_th curve; here n=N/2
% X coordinates do not have to be periodic, but the curvature,
% normals, etc that they compute will be garbage.

methods

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
% not arclength
x = X(1:end/2,:);
y = X(end/2+1:end,:);
N = size(x,1);
nv = size(x,2);
IK = o.modes(N);
Dx = o.diffFT(x,IK);
Dy = o.diffFT(y,IK);

end % getDXY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
length = sum(sqrt(Dx.^2 + Dy.^2))*2*pi/N;
area = sum(x.*Dy - y.*Dx)*pi/N;

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
  alpha = (0:N-1)'/N;
end
% Discretization in parameter space

if any(strcmp(options,'curly'))
  a = 1; b = 3*a; c = 0.85; 
  r = 0.5*sqrt( (a*cos(2*pi*alpha)).^2 + (b*sin(2*pi*alpha)).^2) + ...
      .07*cos(12*2*pi*alpha);
  x = scale*c*r.*cos(2*pi*alpha);
  y = scale*r.*sin(2*pi*alpha);
  X0 = [x;y];
  % radius of curly vesicle

elseif any(strcmp(options,'star'))
  radius = 1 + 0.2*cos(folds*2*pi*alpha);
  X0 = [radius.*cos(2*pi*alpha);radius.*sin(2*pi*alpha)];
  % a star that comes very close to intersecting itself at the origin

elseif any(strcmp(options,'ellipse'))
  X0 = [cos(2*pi*alpha);shortax*sin(2*pi*alpha)];
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
%  figure(1); clf;
%  plot(X(1:end/2),X(end/2+1:end),'b--')
%  axis equal;
%  hold on;
  sa = o.diffProp(X); % compute arclength of the vesicle shape
  alpha = o.arc(sa);
  % find parameter values that give discretization points that are
  % nearly equispaced in arclength
  [alpha,X] = o.initConfig(N,'true',options{1},'angle',theta,...
       'reducedArea',ra,'shortax',shortax,'scale',scale,...
       'folds',folds,'parameter',alpha);
end
%clf
%figure(1); hold on;
%plot(X(1:end/2),X(end/2+1:end),'ro')
%axis equal;
%pause

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
IK = o.modes(N);

intsa = o.intFT(sa,IK);

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
    falpha = real(sum(fh.*exp(IK*alpha(j)))) + ...
        ((j-1)*length/N - alpha(j)*length);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
  rcon = concentration + smallper*10*(rand(N,1) - 0.5);
  % uniform concentration with a little random noise
elseif symmetry == 0
  rcon = concentration + 3*smallper*cos(2*pi*alpha) + ...
    0.5*smallper*cos(3*2*pi*alpha) + ...
    0.5*smallper*cos(4*2*pi*alpha);
  % uniform concentration with a little noise at specific even Fourier
  % modes
elseif symmetry==1
  rcon = concentration + 5*smallper*cos(2*2*pi*alpha);
  % uniform concentration with one small even Fourier mode
elseif symmetry==2
  rcon = concentration+ 5*smallper*sin(2*pi*alpha);
  % uniform concentration with one small odd Fourier mode
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
if abs(Dx(1)) > 1e-12
  t0 = atan(Dy(1)/Dx(1));
else
  t0 = pi/2;
end
% initial opening angle
[DDx,DDy] = o.getDXY([Dx;Dy]);
% second derivative of the shape

L = sum(sqrt(Dx.^2 + Dy.^2))/N;
% arclength which should be nearly constant since we constructed
% discretizations points that are equally spaced in arclength

cur = (Dx.*DDy - Dy.*DDx)/L^3;
% compute curvature

IK = o.modes(N);
theta = L*o.intFT(cur,IK);
% integrate the curature to find the opening angle

theta = theta + t0;
% add in the initial opening angle

end % computeOpeningAngle


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = recon(o,N,x0,y0,L,theta)
% X = recon(N,x0,y0,L,theta) computes the (x,y) coordinates of a
% geometry with length L, opening angle theta, and its first point at
% (x0,y0)

IK = o.modes(N);
x = x0 + o.intFT(L*cos(theta),IK);%/2/pi;
y = y0 + o.intFT(L*sin(theta),IK);%/2/pi;
X = o.setXY(x,y);

end % recon


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function df = diffFT(o,f,IK)
% df = diffFT(f,IK) Computes the first derivative of a periodic function
% f using fourier transform. IK is used to speed up the code.  It is the
% index of the fourier modes so that fft and ifft can be used

df = real(ifft(IK.*fft(f)));

end % diffFT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function intf = intFT(o,f,IK)
% intf = intFT(f,IK) computes the indefinite integral of a periodic
% function f using fourier transform

N = numel(f);
fh = fft(f,N);

zeroMode = -sum(fh(2:end)./IK(2:end));
% zero mode of the integral of f
intfh = [zeroMode;fh(2:end)./IK(2:end)];
% non-zero modes of the integral of f

intf = real(fh(1)/N * (0:N-1)'/N + ifft(intfh,N));
% sum of the linear term and the contribution from the periodic part

% disp('Here')
% plot(intf)
% pause
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
%compute fourier derivative of opening tangent angle using "trick" by
%subtracting a function that goes from 0 to 2*pi over the range of
%alpha values
% This routine is equation (61) in the Sohn et al 2010 JCP paper
N = ves.N;
IK = o.modes(N);

dthetads = o.rmLinearPart(ves);
dcur = intFT(o,[dthetads.*unt],IK);

end %fdcur

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dx = rmLinearPart(o,ves)
N = ves.N;
theta = ves.theta;
IK = o.modes(N);

% function with the same increase over a period of 2*pi as the angle
% theta
dx = o.diffFT([theta - 2*pi*(0:1:N-1)'/N],IK) + 2*pi;

end %rmLinearPart

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fntheta = fthetaim(o,ves,unt,utn)
N = ves.N;
L = ves.L;  
theta = ves.theta;
IK = o.modes(N);
%compute right hand side of the ODE for \theta in equation (8) in Sohn
%et al JCP 2010
ftheta = o.forceTheta(ves,unt,utn);   
%now we must subtract off the stiffest part in Fourier space.
%To the power of 3 term comes from the third-derivative of the function
%that L is being applied to in equation (53)
rlen = ves.bendsti*(2*pi*[0 1:N/2 N/2-1:-1:1]'/L).^3/4;
%fntheta - ???
var1 = fft(ftheta,N);
var2 = fft(theta-[2*pi*(0:N-1)]'/N,N);
fntheta = real(ifft(var1+rlen.*var2));

%Krasney filter applied to smooth out spurious high frequency terms
fntheta = o.kfilter(fntheta,N);
% plot(fntheta)
% pause
end %fthetaim

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ftheta = forceTheta(o,ves,unt,utn)

N = ves.N;
L = ves.L;
IK = o.modes(N);
%compute first derivative of the opening angle using Fourier
%differentiation. This is actually the curvature
dthetads = o.rmLinearPart(ves);
%compute the derivative of the normal velocity
dunds = o.diffFT(unt,IK);    
%Check T_s = -V * curvature (local inextensibility condition) (see
%equation (20) in Sohn et al JCP 2010)
ftheta = -dunds/L + dthetads.*utn/L;    
    
end %forceTheta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function filtered = kfilter(o,ftheta,N)
%fourier filter with krasny filtering with 1 input vector 
%USING FAST FOURIER TRANSFORM
tol = 1e-10;
b = fft(ftheta,N);
b(abs(b)/N < tol) = 0;
b(N/2:N/2+2) = 0;

c = real(b.*exp(-10*([0:2:N N-2:-2:2]'/N).^25));
d = imag(b.*exp(-10*([1:2:N-1 0 N-1:-2:3]'/N).^25));
d(N/2-1) = 0; d(N/2+3) = 0;

filtered = real(ifft(c+1i*d,N));

end %kfilter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function N2Hat = frconim(o,ves,epsch,consta)
%This routine is about solving equation (75) for N2Hat.
N = ves.N;
L = ves.L;
rcon = ves.rcon;
IK = o.modes(N);
%setting up RHS of eq (68) 
%termd is the a/eps(f'(u)-eps^2u_ss) + b'(u)/2 * kappa^2 terms in eq (13) 
%ie. termd is the variatiaonl derivative of the energy with
%respect to the lipid concentration
termd = o.fluxj(ves,epsch,consta);
%now taking derivative of the variational derivative twice
rcons = o.diffFT(termd,IK);
rconss = o.diffFT(rcons,IK);
%invPe is the 1/Pe term in eq (67)
invPe = 1;
%fcon is the part of N_2 in eq (67) but is still missing the
%a\eps/s_al^4 * u_alalalal term, (al is shorthand for alpha here)
fcon = invPe*rconss/L^2;
fconHat = fft(fcon);
rconHat = fft(rcon,N);
%FM are the fourier modes
FM = 2*pi*[0 1:N/2 N/2-1:-1:1]';
%rlen is the term that needs to multiply the fourier coefficients of
%the lipid species to result in the fourth derivative scaled by the a
%and \eps constants in equation (67)
rlen = (FM/L).^4*consta*epsch;
%N2Hat is the Fourier coefficeints of the right hand side of
%equation (67) in the physical space.
N2Hat = real(ifft(fconHat + rlen.*rconHat));

end %frconim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RHSflux = fluxj(o,ves,epsch,consta)
%This routine computes the forcing on the RHS flux OF THE Cahn_Hillard
%Equation for rcon periodic on [0,1] and for closed curves. RHSflux are the
%a/eps(f'(u)-eps^2u_ss) + b'(u)/2 * kappa^2 terms in eq (13) 

N = ves.N;
L = ves.L;
rcon = ves.rcon;
IK = o.modes(N);

%dfu is the derivative of the double-well potential
dfu = 0.5*(rcon.*(ones(N,1)-rcon).^2 - rcon.^2.*(ones(N,1)-rcon));
%derivative of the concentration 
rcons = o.diffFT(rcon,IK);
%second derivative of the concentration
rconss = o.diffFT(rcons,IK);
%term2 is eps^2 * u_ss as defined in equation (13)
term2 = -epsch^2*rconss/L^2;
%computing the b'(u)/2 * kappa^2 term in eq (13)
b0 = ves.bendsti;
b1 = ves.bendsti*ves.bendratio;
%rbn is b(u).
rbn = b0*(ones(N,1)-rcon) + b1*rcon;
%rbndu is b'(u) 
rbndu = (b1 - b0)*ones(N,1);
%compute the curvature
dkap = o.acurv(ves);
%term3 is b'(u)/2 * kappa^2 in eq (13)
term3 = 0.5*rbndu.*dkap.^2;
%RHSflux is a/eps(f'(u)-eps^2u_ss) + b'(u)/2 * kappa^2 terms in eq (13) 
RHSflux = consta/epsch * (dfu + term2) + term3;

end %fluxj
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dkap = acurv(o,ves)
%This routine computes the curvature from the tangent angle
%formulation

dkap = o.rmLinearPart(ves)/ves.L;

end %acurv

end % methods

end % classdef
