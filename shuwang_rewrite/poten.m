classdef poten
%This class defines single layers for stokes kernels on 2D periodic curves.
%It also defines the matricies that map a density function defined on the
%boundary of a curve to the layer potential evaluated on the curve. This
%class also has the main routine that evaluates layer-potentials using 
%near-singular integration.

properties
N; % points per curve
end % properties

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function po = poten(N)
% o = poten(N) is a constructor the initializes the class
po.N = N;
end % constructor: poten

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function kmatrix = oddEvenMatrix(o)
% kmatrix = oddEvenMatrix builds a matrix that is 0 if the row and
% column have the same parity (both even or odd) and is 1 otherwise.
% This matrix is used to do odd-even integration with kernels that have
% a weak singularity at the diagonal terms such as the Stokes
% single-layer potential

kmatrix = ones(2,2*o.N);
kmatrix(1,1:2:end) = 0;
kmatrix(2,2:2:end) = 0;
kmatrix = repmat(kmatrix,o.N,1);

end % oddEvenMatrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function selfmatrix = StokesMatrixLogless(o,ves)
% selfmatrix = StokesMatrixLogless builds the Stokes matrix without the
% log singularity. Contains the rightmost kernel in equation (43).

X = ves.X;
% define points
pts = X(1:o.N) + 1i*X(o.N+1:end);
% r is the target minus source for the single-layer potential with ones
% on the main diagonal. The x-coordinate is in the real part and the y-
% coordinate is in the imaginary part
r = repmat(pts,[1 o.N]) - repmat(transpose(pts),[o.N 1]) + 1*eye(o.N);
% rt is (alpha - alpha')/2 in equation (43) where alpha and alpha' go
% from 0 to 2*pi. We are only going from 0 to pi to account for the
% divided by 2 in the argument of sin in equation (43). We put ones on
% the diagonal, the limit which is something smooth (see equation (44))
rt = repmat((1:o.N)'*pi/o.N , [1 o.N]) - ...
     repmat((1:o.N)*pi/o.N, [o.N 1]) + 1*eye(o.N);
% Compute the denominator of equation (43). Putting ones in the diagonal
% of r and denom results in taking log of something non-zero.
denom43 = 2*abs(sin(rt));
irr = 1./(conj(r).*r);    % 1/r^2
d1 = real(r); % x-coordinate of r
d2 = imag(r); % y-coordinate of r
Ilogr = log(abs(r)./denom43);  % log(1/r) diag block
salpha = ves.L/(2*pi); % limiting value for log part of kernel
Ilogr(1:o.N+1:end) = log(salpha); % correct the diagonal terms
oc = curve(o.N); %shorthand for curve class
% Compute the derivatives of x and y
[dx,dy] = oc.getDXY(ves.X);
% Divide by length to get the arclength derivative
dx = dx/ves.L; dy = dy/ves.L;
A11 = d1.^2.*irr; % (1,1) block ie. rx^2/(rx^2 + ry^2)
A12 = d1.*d2.*irr; % (1,2) block ie. rx*ry/(rx^2 + ry^2)
A22 = d2.^2.*irr; % (2,2) block ie. ry^2/(rx^2 + ry^2)
% Limiting values of the above 3 blocks
A11(1:o.N+1:end) = dx.^2;
A12(1:o.N+1:end) = dx.*dy;
A22(1:o.N+1:end) = dy.^2;
% Compute the selfmatrix
selfmatrix = [-Ilogr + A11 A12; A12 -Ilogr + A22];

end % StokesMatrixLogless

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Symm_sigma = IntegrateLogKernel(o,sigma)
%IntegrateLogKernel returns the symms operator applied to sigma, where
%sigma is a density function.

%take the fft of sigma
sigmah = fft(sigma); 
%define the Fourier coefficients
if length(sigma) == 2*o.N
    coeff = [[(0:o.N/2-1)';(-o.N/2:-1)'];[(0:o.N/2-1)';(-o.N/2:-1)']];
else
    coeff = [(0:o.N/2-1)';(-o.N/2:-1)'];
end
Symm_sigmah = -sigmah./abs(coeff);
Symm_sigmah(1) = 0;
Symm_sigma = real(ifft(Symm_sigmah));

end %IntegrateLogKernel


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vel = StokesSLP(o,geom,StokesMat,tau)
% compute the velocity due to the single-layer potential of the density
% function tau

N = geom.N;
vel = zeros(2*N,1);
L = geom.L;

k = StokesMat*tau;

LogKernel1 = o.IntegrateLogKernel(tau(1:N));
LogKernel2 = o.IntegrateLogKernel(tau(N+1:end));

vel(1:N) = L/(4*pi*N)*k(1:N) - L/(8*pi)*LogKernel1;
vel(N+1:end) = L/(4*pi*N)*k(N+1:end) - L/(8*pi)*LogKernel2;

end % StokesSLP


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = StokesDLP(o,geom)

oc = curve(geom.N);
% Geometry positions
[x,y] = oc.getXY(geom.X);

% number of points on geometry
N = geom.N;

% tangent vector
tx = cos(geom.theta);
ty = sin(geom.theta);

% curvature
cur = geom.cur';

% target points
xtar = x(:,ones(N,1))';
ytar = y(:,ones(N,1))';

% source points
xsou = x(:,ones(N,1));
ysou = y(:,ones(N,1));

% tangent at sources
txsou = tx';
tysou = ty';

% rx and ry terms
diffx = xtar - xsou;
diffy = ytar - ysou;
% 1 over the distance to the power of 4
rho4 = (diffx.^2 + diffy.^2).^(-2);
% set diagonal terms to 0. Will fix with limiting value involving the
% curvature a bit later
rho4(1:N+1:N.^2) = 0;

kernel = diffx.*(tysou(ones(N,1),:)) - ...
        diffy.*(txsou(ones(N,1),:));
kernel = kernel.*rho4;

% (1,1) component
D11 = kernel.*diffx.^2;
% diagonal limiting term
D11(1:N+1:N.^2) = 0.5*cur.*txsou.^2;

% (1,2) component
D12 = kernel.*diffx.*diffy;
% diagonal limiting term
D12(1:N+1:N.^2) = 0.5*cur.*txsou.*tysou;

% (2,2) component
D22 = kernel.*diffy.^2;
% diagonal limiting term
D22(1:N+1:N.^2) = 0.5*cur.*tysou.^2;

% build matrix with four blocks
D = [D11 D12; D12 D22];
% scale with the arclength spacing and divide by pi
D = 1/pi*D*geom.L/N;


end % StokesDLP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [N0] = StokesN0mat(o, ves)
% N0 = stokesN0matrix(vesicle) generates the the integral operator with
% kernel normal(x) \otimes normal(y) which removes the rank one
% defficiency of the double-layer potential.  Need this operator for
% solid walls

% Use N0 if solving (-1/2 + DLP)\eta = f where f has no flux through
% the boundary.  By solving (-1/2 + DLP + N0)\eta = f, we guarantee
% that \eta also has no flux through the boundary.  This is not
% required, but it means we're enforcing one addition condition on eta
% which removes the rank one kernel.  DLP is the double-layer potential
% for stokes equation

N = ves.N;
X = ves.X;
oc = curve(N);

%get the tangent
[~,tangent,~] = oc.diffProp(X);
%compute the normal
normal = [tangent(N+1:end);-tangent(1:N)];
normal = normal(:,ones(2*N,1));
jac = oc.diffProp(X);
jacobian = [jac;jac];
jacobian = jacobian(:,ones(2*N,1));
N0 = normal.*normal'.*jacobian'*2*pi/N;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START OF ROUTINES THAT EVALUATE STOKES LAYER POTENTIALS AT TARGET
% POINTS THAT DIFFER FROM THE SOURCCE POINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vel = StokesSLPtar(o,ves,trac,Xtar)

oc = curve(ves.N);

[xsou,ysou] = oc.getXY(ves.X);
[fx,fy] = oc.getXY(trac);

[xtar,ytar] = oc.getXY(Xtar);
[itar,jtar] = size(xtar);
xtar = xtar(:); ytar = ytar(:);

velx = zeros(size(xtar));
vely = zeros(size(xtar));


for k = 1:numel(xtar)
  rx = xtar(k) - xsou;
  ry = ytar(k) - ysou;
  rho2 = rx.^2 + ry.^2;
  rdotf = rx.*fx + ry.*fy;

  velx(k) = sum(-0.5*log(rho2).*fx + rdotf./rho2.*rx)/(4*pi);
  vely(k) = sum(-0.5*log(rho2).*fy + rdotf./rho2.*ry)/(4*pi);
end
% reshape so that it is the original size that came into the input
velx = reshape(velx,itar,jtar);
vely = reshape(vely,itar,jtar);

velx = velx * ves.L/ves.N;
vely = vely * ves.L/ves.N;
vel = [velx;vely];

end % StokesSLPtar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vel = StokesDLPtar(o,geom,eta,Xtar)
% evaluate the double-layer potential due to the geometry stored in geom
% using the density function eta, and evaluating at the target points in
% Xtar

oc = curve(geom.N);

[xsou,ysou] = oc.getXY(geom.X);
[~,tangent,~] = oc.diffProp(geom.X);
normal = [tangent(geom.N+1:end);-tangent(1:geom.N)];
[nx,ny] = oc.getXY(normal);

[xtar,ytar] = oc.getXY(Xtar);
[itar,jtar] = size(xtar);
[xeta,yeta] = oc.getXY(eta);
xtar = xtar(:); ytar = ytar(:);
xeta = xeta(:); yeta = yeta(:); 

velx = zeros(size(xtar));
vely = zeros(size(xtar));

for k = 1:numel(xtar)
  rx = xtar(k) - xsou;
  ry = ytar(k) - ysou;
  rho2 = rx.^2 + ry.^2;
  rdotn = rx.*nx + ry.*ny;
  rdoteta = rx.*xeta + ry.*yeta;
  kernel = (rdotn./rho2).*(rdoteta./rho2);
  velx(k) = -1/pi*sum(kernel.*rx);
  vely(k) = -1/pi*sum(kernel.*ry);
end
% reshape so that it is the original size that came into the input
velx = reshape(velx,itar,jtar);
vely = reshape(vely,itar,jtar);

vel = [velx;vely]*geom.L/geom.N;

end % StokesDLPtar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF ROUTINES THAT EVALUATE STOKES LAYER POTENTIALS AT TARGET POINTS
% THAT DIFFER FROM THE SOURCCE POINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LP = nearSingInt(o,vesicleSou,f,selfMat,...
    NearStruct,kernel,kernelDirect,vesicleTar,tEqualS,idebug)
% 
% vesicleSou : Source points
% f : density
% selfMat: function selfMat(f) : self-interactions
% NearStruct : output of getZone()
% kernel : fmm kernel evaluations between source and target points
% kernelDirect: direct kernel evaluations between source and target points
% tEqualsS: source and target points coincide
% idebug: debugging flag to plot near singular potential for all points in the
% near zone
%
% LP =
% nearSingInt(vesicle,f,selfMat,NearStruct,kernel,kernelDirect,vesicleTar,tEqualS,idebug)
% computes a layer potential due to f at all points in vesicleTar.X.  If
% tEqualS==true, then the vesicleTar == vesicleSou and the self-vesicle
% interaction is skipped.  selfMat is the diagonal of the potential
% needed to compute the layer potential of each vesicle indepenedent of
% all others.  kernel and kernelDirect are two (possibly the same)
% routines that compute the layer potential.  kernelDirect always uses
% the direct method whereas kernel may use an FMM-accelerated method.
% NearStruct is a structure containing the variables
% zone,dist,nearest,icp,argnear which are required by near-singular
% integration (they keep everything sorted and precomputed) Everything
% is in the 2*N x nv format Can pass a final argument if desired so that
% plots of the near-singular integration algorithm are displayed

if (tEqualS && size(vesicleSou.X,2) == 1)
  LP = zeros(size(vesicleSou.X));
  return
end
% only a single vesicle, so velocity on all other vesicles will always
% be zero

if (nargin == 9)
  idebug = false;
end

dist = NearStruct.dist;
zone = NearStruct.zone;
nearest = NearStruct.nearest;
icp = NearStruct.icp;
argnear = NearStruct.argnear;

Xsou = vesicleSou.X; % source positions
Nsou = size(Xsou,1)/2; % number of source points
nvSou = size(Xsou,2); % number of source 'vesicles'
Xtar = vesicleTar.X; % target positions
Ntar = size(Xtar,1)/2; % number of target points
nvTar = size(Xtar,2); % number of target 'vesicles'

h = vesicleSou.L/Nsou; % arclength term

Nup = Nsou*ceil(sqrt(Nsou));
vself = selfMat(f);

% upsample to N^(3/2).  
Xup = [interpft(Xsou(1:Nsou,:),Nup);...
       interpft(Xsou(Nsou+1:2*Nsou,:),Nup)];
fup = [interpft(f(1:Nsou,:),Nup);...
       interpft(f(Nsou+1:2*Nsou,:),Nup)];

% Compute velocity due to each vesicle independent of others.  This is
% needed when using near-singular integration since we require a value
% for the layer-potential on the vesicle of sources 

% allocate space for storing velocity at intermediate points needed
% by near-singular integration

vesicleUp = capsules(Xup,[]);
% Build an object with the upsampled vesicle

LP = o.lagrangeInterp;
interpOrder = size(LP,1);
% lagrange interpolation order
p = ceil((interpOrder+1)/2);
% want half of the lagrange interpolation points to the left of the
% closest point and the other half to the right

if tEqualS % sources == targets
  if nvSou > 1
    for k = 1:nvSou
      K = [(1:k-1) (k+1:nvSou)];
      farField(:,k) = kernelDirect(vesicleUp,fup,Xtar(:,k));
    end
    % This is a huge savings if we are using a direct method rather than
    % the fmm to evaluate the layer potential. The speedup is more than
    % N^{1/2}, where N is the resolution of the vesicles that we are
    % computing with
  else
    farField = zeros(2*Ntar,nvTar);
  end

else % sources ~= targets
  farField = kernel(vesicleUp,fup,Xtar);
  % evaluate layer potential due to all 'vesicles' at all points in
  % Xtar;
end
% Use upsampled trapezoid rule to compute layer potential

nearField = zeros(2*Ntar,nvTar);

beta = 1.1;
% small buffer to make sure Lagrange interpolation points are
% not in the near zone
for k1 = 1:nvSou
  if tEqualS % sources == targets
    K = [(1:k1-1) (k1+1:nvTar)];
    % skip diagonal vesicle
  else % sources ~= targets
    K = (1:nvTar);
    % consider all vesicles
  end
  for k2 = K
    J = find(zone{k1}(:,k2) == 1);
    % set of points on vesicle k2 close to vesicle k1
    if (numel(J) ~= 0)
      indcp = icp{k1}(J,k2);
      % closest point on vesicle k1 to each point on vesicle k2 
      % that is close to vesicle k1
      for j = 1:numel(J)
        pn = mod((indcp(j)-p+1:indcp(j)-p+interpOrder)' - 1,Nsou) + 1;
        % index of points to the left and right of the closest point
        v = filter(1,[1 -full(argnear{k1}(J(j),k2))],...
          LP*vself(pn,k1));
        vel(J(j),k2,k1) = v(end);  
        % x-component of the velocity at the closest point
        v = filter(1,[1 -full(argnear{k1}(J(j),k2))],...
          LP*vself(pn+Nsou,k1));
        vel(J(j)+Ntar,k2,k1) = v(end);
        % y-component of the velocity at the closest point
      end
%     compute values of velocity at required intermediate points
%     using local interpolant
      
      if ((numel(J) + numel(fup)) >= 512 && numel(J) > 32)
        potTar = kernel(vesicleUp,fup,...
           [Xtar(J,k2);Xtar(J+Ntar,k2)]);
      else
        potTar = kernelDirect(vesicleUp,fup,...
           [Xtar(J,k2);Xtar(J+Ntar,k2)]);
      end
      % Need to subtract off contribution due to vesicle k1 since its
      % layer potential will be evaulted using Lagrange interpolant of
      % nearby points
      nearField(J,k2) =  nearField(J,k2) - ...
          potTar(1:numel(J));
      nearField(J+Ntar,k2) =  nearField(J+Ntar,k2) - ...
          potTar(numel(J)+1:end);
      
      XLag = zeros(2*numel(J),interpOrder - 1);
      % initialize space for initial tracer locations
      for i = 1:numel(J)
        nx = (Xtar(J(i),k2) - nearest{k1}(J(i),k2))/...
            dist{k1}(J(i),k2);
        ny = (Xtar(J(i)+Ntar,k2) - nearest{k1}(J(i)+Ntar,k2))/...
            dist{k1}(J(i),k2);
        XLag(i,:) = nearest{k1}(J(i),k2) + ...
            beta*h*nx*(1:interpOrder-1);
        XLag(i+numel(J),:) = nearest{k1}(J(i)+Ntar,k2) + ...
            beta*h*ny*(1:interpOrder-1);
        % Lagrange interpolation points coming off of vesicle k1 All
        % points are behind Xtar(J(i),k2) and are sufficiently far from
        % vesicle k1 so that the Nup-trapezoid rule gives sufficient
        % accuracy
      end

      if (numel(XLag)/2 > 100)
        lagrangePts = kernel(vesicleUp,fup,XLag);
      else
        lagrangePts = kernelDirect(vesicleUp,fup,XLag);
      end
      % evaluate velocity at the lagrange interpolation points
      
      for i = 1:numel(J)
        Px = LP*[vel(J(i),k2,k1) lagrangePts(i,:)]';
        Py = LP*[vel(J(i)+Ntar,k2,k1) lagrangePts(i+numel(J),:)]';
        % Build polynomial interpolant along the one-dimensional
        % points coming out of the vesicle
        dscaled = full(dist{k1}(J(i),k2)/(beta*h*(interpOrder-1)));
        % Point where interpolant needs to be evaluated

        v = filter(1,[1 -dscaled],Px);
        nearField(J(i),k2) = nearField(J(i),k2) + ...
            v(end);
        v = filter(1,[1 -dscaled],Py);
        nearField(J(i)+Ntar,k2) = nearField(J(i)+Ntar,k2) + ...
            v(end);
        % Assign higher-order results coming from Lagrange 
        % integration to velocity at near point.  Filter is faster
        % than polyval

        if idebug
          figure(2); clf; hold on;
          plot(Xsou(1:Nsou,:),Xsou(Nsou+1:end,:),'r.','markersize',10)
          plot(Xtar(1:Ntar,:),Xtar(Ntar+1:end,:),'k.','markersize',10)
          plot(Xtar(J,k2),Xtar(Ntar+J,k2),'b.','markersize',10)
          plot(XLag(1:numel(J),:),XLag(numel(J)+1:end,:),'kx','markersize',10)
          plot(XLag(i,:),XLag(numel(J)+i,:),'gx','markersize',10)
          axis equal

          figure(1); clf; hold on
          plot((0:interpOrder-1)*beta*h,...
              real([vel(J(i),k2,k1) lagrangePts(i,:)]),'g-o')
          plot((0:interpOrder-1)*beta*h,...
              real([vel(J(i)+Ntar,k2,k1) lagrangePts(i+numel(J),:)]),'r--o')
          pause
        end
        % DEBUG: PASS IN idebug=true INTO THIS ROUTINE AND THEN YOU CAN
        % SEE THE INTERPOLATION POINTS AND CHECK THE SMOOTHNESS OF THE
        % INTERPOLANT

      end % i
    end % numel(J) ~= 0
    % Evaluate layer potential at Lagrange interpolation
    % points if there are any
  end % k2
end % k1
% farField

LP = farField + nearField;

% **********************************************************************
% Add kernel due to far points and near points.  Far points were
% upsampled if source==vesicle so need to truncate here.  We are 
% only using Ntar target points.  Note that it is only the sources 
% that were upsampled

end % nearSingInt





end % methods


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static)

function LP = lagrangeInterp
% interpMap = lagrangeInterp builds the Lagrange interpolation
% matrix that takes seven function values equally distributed
% in [0,1] and returns the seven polynomial coefficients

interpMat = zeros(7);
LP(1,1) = 6.48e1;
LP(1,2) = -3.888e2;
LP(1,3) = 9.72e2;
LP(1,4) = -1.296e3;
LP(1,5) = 9.72e2;
LP(1,6) = -3.888e2;
LP(1,7) = 6.48e1;

LP(2,1) = -2.268e2;
LP(2,2) = 1.296e3;
LP(2,3) = -3.078e3;
LP(2,4) = 3.888e3;
LP(2,5) = -2.754e3;
LP(2,6) = 1.0368e3;
LP(2,7) = -1.62e2;

LP(3,1) = 3.15e2;
LP(3,2) = -1.674e3;
LP(3,3) = 3.699e3;
LP(3,4) = -4.356e3;
LP(3,5) = 2.889e3;
LP(3,6) = -1.026e3;
LP(3,7) = 1.53e2;

LP(4,1) = -2.205e2;
LP(4,2) = 1.044e3;
LP(4,3) = -2.0745e3;
LP(4,4) = 2.232e3;
LP(4,5) = -1.3815e3;
LP(4,6) = 4.68e2;
LP(4,7) = -6.75e1;

LP(5,1) = 8.12e1;
LP(5,2) = -3.132e2;
LP(5,3) = 5.265e2;
LP(5,4) = -5.08e2;
LP(5,5) = 2.97e2;
LP(5,6) = -9.72e1;
LP(5,7) = 1.37e1;

LP(6,1) = -1.47e1;
LP(6,2) = 3.6e1;
LP(6,3) = -4.5e1;
LP(6,4) = 4.0e1;
LP(6,5) = -2.25e1;
LP(6,6) = 7.2e0;
LP(6,7) = -1e0;

LP(7,1) = 1e0;
% rest of the coefficients are zero

end % lagrangeInterp

end % methods (Static)


end % classdef

