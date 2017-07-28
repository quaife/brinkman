classdef poten 
% this class defines single and double layers for various kernels 
% (stokes, laplace) on 2D periodic curves.  Also defines the 
% integrals required for pressure and stress.
% Defines the matricies that map a density function defined on the
% boundary of a curve to the layer potential evaluated on the 
% curve, and defines the operators that takes a density function 
% defined on the boundary of a curve and returns the layer 
% potential at arbitrary target points.
% This class also has the main routine that evaluates layer
% potentials using near-singular integration.
    
properties
  qw; 
  % quadrature weights for logarithmic singularity
  qp; 
  % quadrature points for logarithmic singularity (Alpert's rule)
  interpMat;  
  % interpolation matrix used for near-singular integration
  % This matrix replaces the need to use polyfit
  G
  % single-layer potential due to the inner geometries
  invG
  % inverse of the matrix G required in the preconditioner
  fmm
  % use the fmm or not
end % properties

methods 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = poten(geom,fmm,bieSolve,computeEuler)
% o = poten(N,fmm): constructor; N is the number of points per curve.

o.interpMat = o.lagrangeInterp;
% load in the interpolation matrix which is precomputed
% with 7 interpolation points.  This is required for the near-singular
% interation scheme
accuracyOrder = 8;
o.qw = o.quadratureS(accuracyOrder,geom.N);
o.qp = o.qw(:,2:end);
o.qw = o.qw(:,1);
if bieSolve || computeEuler
  o.G = o.stokesSLmatrix(geom);
else
  o.G = [];
end

if bieSolve
  for k = 1:geom.nv
    o.invG(:,:,k) = pinv(o.G(:,:,k));
  end
else
  o.invG = zeros(size(o.G));
end

o.fmm = fmm;
end % poten: constructor



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function invGf = matVecInvBD(o,f,innerGeom,DLPpreco)
% invGf = matVecInvBD(f,innerGeom,DLPpreco) is used as the
% block-diagonal preconditioner due to the inner circles and outer
% geometry.  Does inverse of each term using a circle who has the same
% circumference as the geometry, and the outer boundary is
% preconditioned with the precomputed matrix DLPpreco

invGf = f;
% want the part corresponding to the outer walls to be the identity

sigmah = zeros(2*innerGeom.N,1);
modes = [(0:innerGeom.N/2-1) (-innerGeom.N/2:-1)]';
for k = 1:innerGeom.nv
  istart = (k-1)*2*innerGeom.N+1;
  iend = istart + 2*innerGeom.N - 1;
  invGF(istart:iend) = o.invG(:,:,k)*f(istart:iend);
end % k = exclusions
istart = iend + 1;
iend = numel(f);
invGf(istart:iend) = DLPpreco * f(istart:iend);
% use precomputed matrix DLPpreco to precondition the outer boundary
% density function

invGf = real(invGf);

end % matVecInvBD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Gf = matVecMultiply(o,f,innerGeom,outerGeom,...
    NearI2I,NearI2O,NearO2I)
% Gf = matVecMultiply(f,innerGeom,outerGeom,NearI2O,NearO2I) is the
% main confined flow matrix-vector-multiplication routine.  f is the
% density function, innerGeom and outerGeom are objects corresponding
% to the inner and outer boundaries, respectively, and NearI2O and
% NearO2I are near-singular integration structures required to do
% inner to outer (I2O) and outer to inner (O2I) interactions
global iteration
iteration = iteration + 1;
fprintf('Iteration number is %d\n',iteration);

Ninner = innerGeom.N;
nv = innerGeom.nv;
Nouter = outerGeom.N;
Gfinner = zeros(2*Ninner,nv);
Gfouter = zeros(2*Nouter,1);
innerEta = zeros(2*Ninner,nv);
outerEta = zeros(2*Nouter,1);
% allocate space for density function and layer potentials

for k = 1:nv
  istart = (k-1)*2*Ninner + 1;
  iend = istart + 2*Ninner - 1;
  innerEta(:,k) = f(istart:iend);
end
istart = nv*2*Ninner + 1;
iend = istart + 2*Nouter - 1;
outerEta = f(istart:iend);
% unstack f so that it is one x-coordinate and one y-coordinate per
% column

Gfinner = Gfinner + o.exactStokesSLdiag(innerGeom,innerEta);
% diagonal term from exclusions

Gfouter = Gfouter - 0.5*outerEta;
% jump term
Gfouter = Gfouter + o.exactStokesDLdiag(outerGeom,outerEta);
% double-layer potential
Gfouter = Gfouter + o.exactStokesN0diag(outerGeom,outerEta);
% rank one modification to remove null space

if ~o.fmm
  stokesSLP = o.nearSingInt(...
      innerGeom,innerEta,@o.exactStokesSLdiag,...
      NearI2I,@o.exactStokesSL,innerGeom,1,'inner');
  % velocity on inner posts due to all non-self inner posts
  stokesSLPtar = o.nearSingInt(...
      innerGeom,innerEta,@o.exactStokesSLdiag,...
      NearI2O,@o.exactStokesSL,outerGeom,0,'inner');
  % velocity on outer geometry due to all inner posts
  stokesDLPtar = o.nearSingInt(...
      outerGeom,outerEta,@o.exactStokesDLdiag,...
      NearO2I,@o.exactStokesDL,innerGeom,0,'outer');
  % velocity on inner posts due to outer geometry
  % velocity on inner posts due to outer geometry
else
  stokesSLP = o.nearSingInt(...
      innerGeom,innerEta,@o.exactStokesSLdiag,...
      NearI2I,@o.exactStokesSLfmm,innerGeom,1,'inner');
  % velocity on inner posts due to all non-self inner posts
  stokesSLPtar = o.nearSingInt(...
      innerGeom,innerEta,@o.exactStokesSLdiag,...
      NearI2O,@o.exactStokesSLfmm,outerGeom,0,'inner');
  % velocity on outer geometry due to all inner posts
  stokesDLPtar = o.nearSingInt(...
      outerGeom,outerEta,@o.exactStokesDLdiag,...
      NearO2I,@o.exactStokesDLfmm,innerGeom,0,'outer');
  % velocity on inner posts due to outer geometry
end

Gfinner = Gfinner + stokesSLP;
% add in contribution from all other exclusions
Gfouter = Gfouter + stokesSLPtar;
% add in contribution from all exclusions to outer boundary

Gfinner = Gfinner + stokesDLPtar;
% add in contribution from the outer boundary to all the exclusions

theta = (0:Ninner-1)'*2*pi/Ninner;
for k = 1:nv
  Gfinner(:,k) = Gfinner(:,k) + (2*pi/Ninner)^2*...
      ([cos(theta);sin(theta)]'*innerEta(:,k))*[cos(theta);sin(theta)];
end
% rank one modification to remove null space.  Each exclusion results
% in a one dimensional null space

Gf = [Gfinner(:);Gfouter(:)];

end % matVecMultiply



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Gf = exactStokesSLdiag(o,geom,f)
% Gf = exactStokesSLdiag(geom,f) computes the single-layer potential due
% to a bunch of circles in the object geom.  f is the density function.

Gf = zeros(2*geom.N,geom.nv);

for k = 1:geom.nv
  Gf(:,k) = o.G(:,:,k) * f(:,k);
end % k = exclusions

end % exactStokesSLdiag



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stokesDLP = exactStokesDLdiag(o,geom,f)
% stokesDLP = exactStokesDLdiag(geom,f) computes the double-layer
% potential when the source and target points coincide.  f is the
% density function, the trapezoid rule is used, and the correct
% limiting term is used for the diagonal term.  The result will have
% spectral accuracy as the kernel is smooth

oc = curve;
[x,y] = oc.getXY(geom.X);
[nx,ny] = oc.getXY(geom.normal);
[tx,ty] = oc.getXY(geom.xt);
[fx,fy] = oc.getXY(f);
N = size(x,1);
sa = geom.sa;
fDOTt = fx.*tx + fy.*ty;

stokesDLP = zeros(2*N,1);

for j = 1:N
  rho2 = (x(j)-x).^2 + (y(j)-y).^2;

  coeff = ((x(j) - x).*nx + (y(j) - y).*ny).*...
          ((x(j) - x).*fx + (y(j) - y).*fy).*...   
          sa./rho2.^2/pi*2*pi/N;
  coeff(j) = 0;
  % Set diagonal term to one to avoid dividing by zero

  stokesDLP(j) = sum(coeff.*(x(j) - x));
  stokesDLP(j+N) = sum(coeff.*(y(j) - y));
  stokesDLP(j) = stokesDLP(j) - ...
      1/2/pi*geom.cur(j)*fDOTt(j)*tx(j)*sa(j)*2*pi/N;
  stokesDLP(j+N) = stokesDLP(j+N) - ...
      1/2/pi*geom.cur(j)*fDOTt(j)*ty(j)*sa(j)*2*pi/N;
end

end % exactStokesDLdiag


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stokesN0 = exactStokesN0diag(o,geom,f)
% stokesN0 = exactStokesN0diag(geom,f) computes the rank one
% modification to the double-layer potential to remove the null
% space.  f is the density function.

oc = curve;
[nx,ny] = oc.getXY(geom.normal);
[fx,fy] = oc.getXY(f);
N = size(nx,1);
sa = geom.sa;
INTnDOTf = sum((nx.*fx + ny.*fy).*sa)*2*pi/N;

stokesN0 = zeros(2*N,1);

for j = 1:N
  stokesN0(j) = INTnDOTf*nx(j);
  stokesN0(j+N) = INTnDOTf*ny(j);
end

end % exactStokesN0diag


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function velAndDef = interpolateLayerPot(o,t,zIn,...
    eulerX,eulerY,u,v,u_x,u_y,v_x,v_y,T,...
    xmthresh,xpthresh,ymthresh,ypthresh,defGradient)
% vel = interpolateLayerPot(t,zIn,eulerX,eulerY,u,v,T) interpolates a
% given velocity field (u,v) defined at the points (eulerX,eulerY)
% at the set of points defined in zIn.  t is the current time and T
% is the time horizion which is used to print a progress bar

%message = ['ode45 ' num2str(t/T*100,'%04.1f') ' %% completed '];
%nmess = numel(message);
%fprintf(repmat('\b',1,nmess));
%fprintf(message);
%% print how far along the simulation has prcoeeded

x = zIn(1); y = zIn(2);
z1 = zIn(3); z2 = zIn(4); z3 = zIn(5); z4 = zIn(6);

if x <= xmthresh
  % tracer location is before the point where simulation is stopped
  velx = 0;
  vely = 0;
  F11 = 0; F12 = 0; F21 = 0; F22 = 0;
elseif x >= xpthresh
  % tracer location is past the point where simulation is stopped
  velx = 0;
  vely = 0;
  F11 = 0; F12 = 0; F21 = 0; F22 = 0;
elseif y <= ymthresh
  % tracer location is outside of the solid wall
  velx = 0;
  vely = 0;
  F11 = 0; F12 = 0; F21 = 0; F22 = 0;
elseif y >= ypthresh
  % tracer location is outside of the solid wall
  velx = 0;
  vely = 0;
  F11 = 0; F12 = 0; F21 = 0; F22 = 0;
else
  % tracer location is fine, so we do an interpolation
  [ny,nx] = size(eulerX);

%  [~,imaxX] = min(abs(eulerX(1,:) - (x + 0.1)));
%  [~,iminX] = min(abs(eulerX(1,:) - (x - 0.1)));
%  [~,imaxY] = min(abs(eulerY(:,1) - (y + 0.1)));
%  [~,iminY] = min(abs(eulerY(:,1) - (y - 0.1)));
% This is too slow

  xmin = eulerX(1,1);
  dx = eulerX(1,2) - xmin;
  ymin = eulerY(1,1);
  dy = eulerY(2,1) - ymin;
  a = ceil(((x + 0.1)-xmin)/dx);
  b = ceil(((x - 0.1)-xmin)/dx);
  c = ceil(((y + 0.1)-ymin)/dy);
  d = ceil(((y - 0.1)-ymin)/dy);
  imaxX = ceil((a + b)/2) + 5;
  iminX = ceil((a + b)/2) - 5;
  imaxY = ceil((c + d)/2) + 5;
  iminY = ceil((c + d)/2) - 5;
  % Take a very minimal window to help with interpolation that is done
  % next

  imaxX = max(10,imaxX);
  iminX = min(nx-10,iminX);
  imaxY = max(10,imaxY);
  iminY = min(ny-10,iminY);
  % find a small window that contains the interpolation points

%  tic
  velx = interp2FAST(eulerX(iminY:imaxY,iminX:imaxX),...
                     eulerY(iminY:imaxY,iminX:imaxX),...
                     u(iminY:imaxY,iminX:imaxX),x,y);
  vely = interp2FAST(eulerX(iminY:imaxY,iminX:imaxX),...
                     eulerY(iminY:imaxY,iminX:imaxX),...
                     v(iminY:imaxY,iminX:imaxX),x,y);
%  toc
%  tic
%  velx = interp2(eulerX(iminY:imaxY,iminX:imaxX),...
%                     eulerY(iminY:imaxY,iminX:imaxX),...
%                     u(iminY:imaxY,iminX:imaxX),x,y,'cubic');
%  vely = interp2(eulerX(iminY:imaxY,iminX:imaxX),...
%                     eulerY(iminY:imaxY,iminX:imaxX),...
%                     v(iminY:imaxY,iminX:imaxX),x,y,'cubic');
%  toc
  % interpolate the velocity field

  if defGradient
    velx_x = interp2FAST(eulerX(iminY:imaxY,iminX:imaxX),...
                       eulerY(iminY:imaxY,iminX:imaxX),...
                       u_x(iminY:imaxY,iminX:imaxX),x,y);
    velx_y = interp2FAST(eulerX(iminY:imaxY,iminX:imaxX),...
                       eulerY(iminY:imaxY,iminX:imaxX),...
                       u_y(iminY:imaxY,iminX:imaxX),x,y);
    vely_x = interp2FAST(eulerX(iminY:imaxY,iminX:imaxX),...
                       eulerY(iminY:imaxY,iminX:imaxX),...
                       v_x(iminY:imaxY,iminX:imaxX),x,y);
    vely_y = interp2FAST(eulerX(iminY:imaxY,iminX:imaxX),...
                       eulerY(iminY:imaxY,iminX:imaxX),...
                       v_y(iminY:imaxY,iminX:imaxX),x,y);
  else
    velx_x = 0; velx_y = 0;
    vely_x = 0; vely_y = 0;
  end
  % interpolate the gradient of the velocity field

  F11 = velx_x*z1 + velx_y*z3;
  F12 = velx_x*z2 + velx_y*z4;
  F21 = vely_x*z1 + vely_y*z3;
  F22 = vely_x*z2 + vely_y*z4;

end

velAndDef = [velx;vely;F11;F12;F21;F22];
% stack the output appropriately

if velx~=velx
  fprintf('\n Problem with Interpolant\n');
end


end % interpolateLayerPot



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vel = layerEvalVelocity(o,Xtar,...
    innerGeom,outerGeom,sigmaInner,sigmaOuter)
% vel = layerEvalVelocity(t,Xtar,xmThresh,xpThresh,innerGeom,outerGeom,
% sigmaInner,sigmaOuter) computes the velocity due to the inner and
% outer boundaries at a set of fixed Eulerian points.  t is the
% current time (not needed), Xtar is the set of target points,
% xmThresh and xpThresh are upper and lower bounds on the y
% componenet.  Target points whose y component lies outside the window
% [xmThresh,xpThresh] is automatically assigned a velocity of 0

targetPnts = bodies(Xtar);
% Build an object for the target points

[~,NearI2T] = innerGeom.getZone(targetPnts,2);
[~,NearO2T] = outerGeom.getZone(targetPnts,2);
% near singular integration structures due to the inner exclusions and
% the outer geometry

if ~o.fmm
  vel1 = o.nearSingInt(innerGeom,sigmaInner,@o.exactStokesSLdiag,...
      NearI2T,@o.exactStokesSL,targetPnts,0,'inner');
  % velocity due to the exclusions
  vel2 = o.nearSingInt(outerGeom,sigmaOuter,@o.exactStokesDLdiag,...
      NearO2T,@o.exactStokesDL,targetPnts,0,'outer');
  % velocity due to the outer wall
else
  vel1 = o.nearSingInt(innerGeom,sigmaInner,@o.exactStokesSLdiag,...
      NearI2T,@o.exactStokesSLfmm,targetPnts,0,'inner');
  % velocity due to the exclusions
  vel2 = o.nearSingInt(outerGeom,sigmaOuter,@o.exactStokesDLdiag,...
      NearO2T,@o.exactStokesDL,targetPnts,0,'outer');
  % velocity due to the outer wall
end

vel = vel1 + vel2;
% add velocity due to the inner and outer boundaries

end % layerEvalVelocity


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LP = nearSingInt(o,souPts,f,diagLP,...
    NearStruct,kernel,tarPts,tEqualS,side,trash)
% LP = nearSingInt(souPts,f,diagLP,...
% NearStruct,kernel,tarPts,tequalS,side,trash) computes a layer
% potential due to f at all points in tarPts.X.  If tEqualS==true,
% then the tarPts == souPts and the self-interaction is skipped.
% diagLP is the diagonal of the potential needed to compute the layer
% potential of each source curve indepenedent of all others.
% NearStruct contains all the variables required to do near-singular
% integration (they keep everything sorted and precomputed) Everything
% is in the 2*N x nv format Can pass a final argument if desired so
% that plots of the near-singular integration algorithm are displayed

dist = NearStruct.dist;
zone = NearStruct.zone;
nearest = NearStruct.nearest;
icp = NearStruct.icp;
argnear = NearStruct.argnear;

Xsou = souPts.X; % source positions
Nsou = size(Xsou,1)/2; % number of source points
nvSou = size(Xsou,2); % number of source curves
Xtar = tarPts.X; % target positions
Ntar = size(Xtar,1)/2; % number of target points
nvTar = size(Xtar,2); % number of target curves

h = souPts.length/Nsou; % arclength term

%Nup = Nsou*ceil(sqrt(Nsou));
Nup = Nsou*ceil(Nsou^(1/4));
% think that upsample by N^{1/4} gives good enough results for
% circles40 upsample to N^(1+\alpha) where a good value of \alpha is up
% in the air.  Nup has to be a multiple of N

vself = diagLP(souPts,f);
% Compute velocity due to each curve independent of others.
% This is needed when using near-singular integration since
% we require a value for the layer-potential on the curve of 
% sources 
if strcmp(side,'outer')
  vself = vself - 0.5*f;
  % if dealing with outer boundary, need to account for the jump in
  % the double-layer potential
end

Xup = [interpft(Xsou(1:Nsou,:),Nup);interpft(Xsou(Nsou+1:2*Nsou,:),Nup)];
fup = [interpft(f(1:Nsou,:),Nup);interpft(f(Nsou+1:2*Nsou,:),Nup)];
% upsample positions, traction jump

souUpPts = bodies(Xup);
% Build an object with the upsampled curve

interpOrder = size(o.interpMat,1);
% lagrange interpolation order
p = ceil((interpOrder+1)/2);
% want half of the lagrange interpolation points to the left of the
% closest point and the other half to the right
vel = zeros(2*Ntar,nvTar,nvSou);
% allocate space for storing velocity at intermediate points
% needed by near-singular integration

for k1 = 1:nvSou
  if tEqualS % sources == targets
    K = [(1:k1-1) (k1+1:nvTar)];
    % skip diagonal curve
  else % sources ~= targets
    K = (1:nvTar);
    % consider all curves
  end
  for k2 = K
    J = find(zone{k1}(:,k2) == 1);
    % set of points on curve k2 close to curve k1
    indcp = icp{k1}(J,k2);
    % closest point on curve k1 to each point on curve k2 
    % that is close to curve k1
    for j = 1:numel(J)
      pn = mod((indcp(j)-p+1:indcp(j)-p+interpOrder)' - 1,Nsou) + 1;
      % index of points to the left and right of the closest point
      v = filter(1,[1 -full(argnear{k1}(J(j),k2))],...
        o.interpMat*vself(pn,k1));
      vel(J(j),k2,k1) = v(end);
      % x-component of the velocity at the closest point
      v = filter(1,[1 -full(argnear{k1}(J(j),k2))],...
        o.interpMat*vself(pn+Nsou,k1));
      vel(J(j)+Ntar,k2,k1) = v(end);
      % y-component of the velocity at the closest point
    end
  end
end
% compute values of velocity at required intermediate points
% using local interpolant

if tEqualS % sources == targets
  farField = kernel(souUpPts,fup);
  % evaluate layer potential at all targets except ignore the
  % diagonal term
else % sources ~= targets
%  tic
  [~,farField] = kernel(souUpPts,fup,Xtar,1:nvSou);
%  toc
%  pause
  % evaluate layer potential due to all curves at all points
  % in Xtar;
end
% Use upsampled trapezoid rule to compute layer potential

nearField = zeros(2*Ntar,nvTar);
% Initialize potential at near points to zero

beta = 1.1;
% small buffer to make sure Lagrange interpolation points are
% not in the near zone

for k1 = 1:nvSou
  if tEqualS % sources == targets
    K = [(1:k1-1) (k1+1:nvTar)];
    % skip diagonal curve
  else % sources ~= targets
    K = (1:nvTar);
    % consider all curves
  end
  for k2 = K
    J = find(zone{k1}(:,k2) == 1);
    if (numel(J) ~= 0)
      [~,potTar] = kernel(souUpPts,fup,...
          [Xtar(J,k2);Xtar(J+Ntar,k2)],k1);
      % Need to subtract off contribution due to curve k1 since
      % its layer potential will be evaulted using Lagrange
      % interpolant of nearby points
      nearField(J,k2) =  - potTar(1:numel(J));
      nearField(J+Ntar,k2) =  - potTar(numel(J)+1:end);

      XLag = zeros(2*numel(J),interpOrder - 1);
      for i = 1:numel(J)
        nx = (Xtar(J(i),k2) - nearest{k1}(J(i),k2))/...
            dist{k1}(J(i),k2);
        ny = (Xtar(J(i)+Ntar,k2) - nearest{k1}(J(i)+Ntar,k2))/...
            dist{k1}(J(i),k2);

        XLag(i,:) = nearest{k1}(J(i),k2) + ...
            beta*h(k1)*nx*(1:interpOrder-1);
        XLag(i+numel(J),:) = nearest{k1}(J(i)+Ntar,k2) + ...
            beta*h(k1)*ny*(1:interpOrder-1);
        % Lagrange interpolation points coming off of curve k1
        % All points are behind Xtar(J(i),k2) and are sufficiently
        % far from curve k1 so that the Nup-trapezoid rule gives
        % sufficient accuracy
      end
      [~,lagrangePts] = kernel(souUpPts,fup,XLag,k1);
      % evaluate velocity at the lagrange interpolation points
      x = XLag(1,:);
      y = XLag(2,:);

      for i = 1:numel(J)
        Px = o.interpMat*[vel(J(i),k2,k1) ...
            lagrangePts(i,:)]';
        Py = o.interpMat*[vel(J(i)+Ntar,k2,k1) ...
            lagrangePts(i+numel(J),:)]';
        % Build polynomial interpolant along the one-dimensional
        % points coming out of the curve 
        dscaled = full(dist{k1}(J(i),k2)/...
            (beta*h(k1)*(interpOrder-1)));
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

        if (nargin == 10)
          figure(2); clf; hold on;
          plot(Xsou(1:Nsou,:),Xsou(Nsou+1:end,:),'r.')
          plot(Xtar(1:Ntar,:),Xtar(Ntar+1:end,:),'k.')
          plot(Xtar(J,k2),Xtar(Ntar+J,k2),'b.')
          plot(XLag(1:numel(J),:),XLag(numel(J)+1:end,:),'kx')
          plot(XLag(i,:),XLag(numel(J)+i,:),'gx')
          axis equal
          axis([Xtar(J(i),k2)-1e-1 Xtar(J(i),k2)+1e-1 ...
                Xtar(J(i)+Ntar,k2)-1e-1 Xtar(J(i)+Ntar,k2)+1e-1]);
%          axis([-2 7 13 17])
%          axis([-1 6 27.5 29])

          figure(1); clf; hold on
          plot((0:interpOrder-1)*beta*h(k1),...
              [vel(J(i),k2,k1) lagrangePts(i,:)],'g-o')
          plot((0:interpOrder-1)*beta*h(k1),...
              [vel(J(i)+Ntar,k2,k1) lagrangePts(i+numel(J),:)],'r--o')
          pause(0.01)
          fprintf('paused')
          pause()
        end
        % DEBUG: PASS IN A DUMMY VARIABLE INTO THIS ROUTINE AND THEN
        % YOU CAN SEE THE INTERPOLATION POINTS AND CHECK THE SMOOTHNESS
        % OF THE INTERPOLANT

      end % i
    end % numel(J) ~= 0
    % Evaluate layer potential at Lagrange interpolation
    % points if there are any
  end % k2
end % k1

if tEqualS % souPts == tarPts
  LP = farField(1:Nup/Ntar:end,:) + nearField;
else % souPts ~= tarPts
  LP = farField + nearField;
end
% Add kernel due to far points and near points.  Far points were
% upsampled if source==target so need to truncate here.  We are 
% only using Ntar target points.  Note that it is only the sources 
% that were upsampled


end % nearSingInt



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stokesSLP,stokesSLPtar] = exactStokesSL(o,geom,f,Xtar,K1)
% [stokesSLP,stokesSLPtar] = exactStokesSL(geom,f,Xtar,K1) computes
% the single-layer potential due to f around all geoms except itself.
% Also can pass a set of target points Xtar and a collection of geoms
% K1 and the single-layer potential due to geoms in K1 will be
% evaluated at Xtar.  Everything but Xtar is in the 2*N x nv format
% Xtar is in the 2*Ntar x ncol format

X = geom.X; % Vesicle positions
sa = geom.sa; % Arclength term

if nargin == 5
  Ntar = size(Xtar,1)/2;
  ncol = size(Xtar,2);
  stokesSLPtar = zeros(2*Ntar,ncol);
else
  K1 = [];
  stokesSLPtar = [];
  ncol = 0;
  % if nargin ~= 5, the user does not need  the velocity at arbitrary
  % points
end

den = f.*[sa;sa]*2*pi/geom.N;
% multiply by arclength term

for k2 = 1:ncol % loop over columns of target points
  for j = 1:Ntar % loop over rows of target points
    dis2 = (Xtar(j,k2) - X(1:geom.N,K1)).^2 + ... 
        (Xtar(j+Ntar,k2) - X(geom.N+1:2*geom.N,K1)).^2;
    diffxy = [Xtar(j,k2) - X(1:geom.N,K1) ; ...
        Xtar(j+Ntar,k2) - X(geom.N+1:2*geom.N,K1)];
    % distance squared and difference of source and target location

    logpart = log(dis2)/2;
    % first part of single-layer potential for Stokes
    
    val = logpart.*den(1:geom.N,K1);
    stokesSLPtar(j,k2) = -sum(val(:));
    val = logpart.*den(geom.N+1:2*geom.N,K1);
    stokesSLPtar(j+Ntar,k2) = -sum(val(:));
    % log term in stokes single-layer potential

    rdotf = (diffxy(1:geom.N,:).*den(1:geom.N,K1) + ...
        diffxy(geom.N+1:2*geom.N,:).*...
        den(geom.N+1:2*geom.N,K1))./dis2;
    % second part of single-layer potential for Stokes

    val = rdotf.*diffxy(1:geom.N,:);
    stokesSLPtar(j,k2) = stokesSLPtar(j,k2) + sum(val(:));
    val = rdotf.*diffxy(geom.N+1:2*geom.N,:);
    stokesSLPtar(j+Ntar,k2) = stokesSLPtar(j+Ntar,k2) + sum(val(:));
    % r \otimes r term of the stokes single-layer potential
  end % j

end % k2
% Evaluate single-layer potential at arbitrary target points
stokesSLPtar = 1/(4*pi)*stokesSLPtar;
% 1/4/pi is the coefficient in front of the single-layer potential


stokesSLP = zeros(2*geom.N,geom.nv); % Initialize to zero
if nargin == 3
  for k1 = 1:geom.nv % geom of targets
    K = [(1:k1-1) (k1+1:geom.nv)];
    % Loop over all geoms except k1
    for j=1:geom.N
      dis2 = (X(j,k1) - X(1:geom.N,K)).^2 + ...
          (X(j+geom.N,k1) - X(geom.N+1:2*geom.N,K)).^2;
      diffxy = [X(j,k1) - X(1:geom.N,K) ; ...
          X(j+geom.N,k1) - X(geom.N+1:2*geom.N,K)];
      % distance squared and difference of source and target location

      logpart = log(dis2)/2;
      % first part of single-layer potential for Stokes

      val = logpart.*den(1:geom.N,K);
      stokesSLP(j,k1) = -sum(val(:));
      val = logpart.*den(geom.N+1:2*geom.N,K);
      stokesSLP(j+geom.N,k1) = -sum(val(:));
      % logarithm terms in the single-layer potential

      rdotf = (diffxy(1:geom.N,:).*den(1:geom.N,K) + ...
          diffxy(geom.N+1:2*geom.N,:).*...
          den(geom.N+1:2*geom.N,K))./dis2;
      % second part of single-layer potential for Stokes

      val = rdotf.*diffxy(1:geom.N,:);
      stokesSLP(j,k1) = stokesSLP(j,k1) + sum(val(:));
      val = rdotf.*diffxy(geom.N+1:2*geom.N,:);
      stokesSLP(j+geom.N,k1) = stokesSLP(j+geom.N,k1) + ...
          sum(val(:));
      % r \otimes r term of the single-layer potential

    end % j
  end % k1
  % Evaluate single-layer potential at geoms but oneself
  stokesSLP = 1/(4*pi)*stokesSLP;
  % 1/4/pi is the coefficient in front of the single-layer potential
end % nargin == 3

end % exactStokesSL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stokesDLP,stokesDLPtar] = ...
    exactStokesDL(o,geom,f,Xtar,K1)
% [stokesDLP,stokesDLPtar] = exactStokesDL(geom,f,Xtar,K1) computes
% the double-layer potential due to f around all geoms except itself.
% Also can pass a set of target points Xtar and a collection of geoms
% K1 and the double-layer potential due to geoms in K1 will be
% evaluated at Xtar.  Everything but Xtar is in the 2*N x nv format
% Xtar is in the 2*Ntar x ncol format

nv = geom.nv; % number of geoms
X = geom.X; % Vesicle positions
normal = geom.normal; % Outward normal
sa = geom.sa; % Jacobian

if nargin == 5
  Ntar = size(Xtar,1)/2;
  ncol = size(Xtar,2);
  stokesDLPtar = zeros(2*Ntar,ncol);
else
  K1 = [];
  stokesDLPtar = [];
  ncol = 0;
  Ntar = 0;
  % if nargin ~= 5, the user does not need the velocity at arbitrary
  % points
end

den = f.*[sa;sa]*2*pi/geom.N;

for k2 = 1:ncol % loop over columns of target points
  for j = 1:Ntar % loop over rows of target points
    diffxy = [Xtar(j,k2) - X(1:geom.N,K1) ; ...
        Xtar(j+Ntar,k2) - X(geom.N+1:2*geom.N,K1)];
    dis2 = diffxy(1:geom.N,:).^2 + ...
        diffxy(geom.N+1:2*geom.N,:).^2;
    % difference of source and target location and distance squared

    coeff = (diffxy(1:geom.N,:).*normal(1:geom.N,K1) + ...
      diffxy(geom.N+1:2*geom.N,:).*...
      normal(geom.N+1:2*geom.N,K1))./dis2.^2.* ...
      (diffxy(1:geom.N,:).*den(1:geom.N,K1) + ...
      diffxy(geom.N+1:2*geom.N,:).*...
      den(geom.N+1:2*geom.N,K1));
    % \frac{(r \dot n)(r \dot density)}{\rho^{4}} term

    stokesDLPtar(j,k2) = stokesDLPtar(j,k2) + ...
        sum(sum(coeff.*diffxy(1:geom.N,:)));
    stokesDLPtar(j+Ntar,k2) = stokesDLPtar(j+Ntar,k2) + ...
        sum(sum(coeff.*diffxy(geom.N+1:2*geom.N,:)));
    % r \otimes r term of the single-layer potential
  end % j
end % k2
stokesDLPtar = stokesDLPtar/pi;
% double-layer potential due to geoms indexed over K1 
% evaluated at arbitrary points

if nargin == 3
  stokesDLP = zeros(2*geom.N,nv);
  for k1 = 1:nv
    K = [(1:k1-1) (k1+1:nv)];
    for j=1:geom.N
      diffxy = [X(j,k1) - X(1:geom.N,K) ; ...
          X(j+geom.N,k1) - X(geom.N+1:2*geom.N,K)];
      dis2 = diffxy(1:geom.N,:).^2 + ...
          diffxy(geom.N+1:2*geom.N,:).^2;
      % difference of source and target location and distance squared

      coeff = (diffxy(1:geom.N,:).*normal(1:geom.N,K) + ...
        diffxy(geom.N+1:2*geom.N,:).*...
        normal(geom.N+1:2*geom.N,K))./dis2.^2 .* ...
        (diffxy(1:geom.N,:).*den(1:geom.N,K) + ...
        diffxy(geom.N+1:2*geom.N,:).*...
        den(geom.N+1:2*geom.N,K));
      % \frac{(r \dot n)(r \dot density)}{\rho^{4}} term

      stokesDLP(j,k1) = stokesDLP(j,k1) + ...
          sum(sum(coeff.*diffxy(1:geom.N,:)));
      stokesDLP(j+geom.N,k1) = stokesDLP(j+geom.N,k1) + ...
          sum(sum(coeff.*diffxy(geom.N+1:2*geom.N,:)));
      % double-layer potential for Stokes
    end
  end

  stokesDLP = stokesDLP/pi;
  % 1/pi is the coefficient in front of the double-layer potential
else
  stokesDLP = [];
end
% double-layer potential due to all geoms except oneself

end % exactStokesDL


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stokesSLP = exactStokesSLDirect(o,geom,f)
% stokesSLP = exactStokesSLDirect(geom,f) is a speedup of
% exactStokesSL which can only be used when the source and target
% points coincide.  It computes the single-layer potential using the
% trapezoid rule and assigning a value of 0 to the diagonal term.

N = geom.N; % number of points per geom
nv = geom.nv; % number of geoms
X = geom.X; % geom positions
oc = curve;
[x,y] = oc.getXY(X); % seperate x and y coordinates

[fx,fy] = oc.getXY(f);
stokesSLP = zeros(2*N,nv);

for k = 1:nv
  for j = 1:N
    xSou = x(:,k);
    ySou = y(:,k);
    fxSou = fx(:,k);
    fySou = fy(:,k);

    rho = (x(j,k) - xSou).^2 + (y(j,k) - ySou).^2;
    rho(j) = 1;
    logpart = -1/2*log(rho);
    rdotf = ((x(j,k) - xSou).*fxSou + (y(j,k) - ySou).*fySou)./rho;

    stokesSLP(j,k) = sum(logpart.*fxSou + rdotf.*(x(j,k) - xSou));
    stokesSLP(j+N,k) = sum(logpart.*fySou + rdotf.*(y(j,k) - ySou));
  end
end
stokesSLP = stokesSLP/4/pi;

end % exactStokesSLDirect


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stokesSLP,stokesSLPtar] = ...
    exactStokesSLfmm(o,geom,f,Xtar,K)
% [stokesSLP,stokeSLPtar] = exactStokesSLfmm(geom,f,Xtar,K) uses the FMM
% to compute the single-layer potential due to all geometries except
% itself geometry is a class of object capsules and f is the density
% function NOT scaled by arclength term.  Xtar is a set of points where
% the single-layer potential due to all geometriess in index set K needs
% to be evaulated
global fmms
fmms = fmms + 1;
% count the total number of calls to fmm

oc = curve;
[x,y] = oc.getXY(geom.X); % seperate x and y coordinates

den = f.*[geom.sa;geom.sa]*2*pi/geom.N;

iprec = 4;

if (nargin == 5)
  stokesSLP = [];
else
  [f1,f2] = oc.getXY(den);
  % need to multiply by arclength term.  Seperate it into
  % x and y coordinate

  [u,v] = stokesSLPfmm(f1(:),f2(:),x(:),y(:),...
      x(:),y(:),1,iprec);
  stokesSLP = zeros(2*geom.N,geom.nv); % initialize
  for k = 1:geom.nv
    is = (k-1)*geom.N+1;
    ie = k*geom.N;
    stokesSLP(1:geom.N,k) = u(is:ie);
    stokesSLP(geom.N+1:2*geom.N,k) = v(is:ie);
  end
  % Wrap the output of the FMM into the usual 
  % [[x1;y1] [x2;y2] ...] format

  for k = 1:geom.nv
    [u,v] = stokesSLPfmm(f1(:,k),f2(:,k),...
        x(:,k),y(:,k),x(:,k),y(:,k),1,iprec);
    stokesSLP(:,k) = stokesSLP(:,k) - [u;v];
  end
  % Subtract potential due to each geometry on its own.  Nothing to
  % change here for potential at Xtar
end

if nargin == 3
  stokesSLPtar = [];
else
  [xsou,ysou] = oc.getXY(geom.X(:,K)); 
  % seperate x and y coordinates at geomtries indexed by K
  [Ntar,ncol] = size(Xtar);
  Ntar = Ntar/2;
  [xtar,ytar] = oc.getXY(Xtar);
  % x and y coordinates of the target points
  [f1,f2] = oc.getXY(den(:,K));
  % seperate x and y coordinates at geomtries indexed by K

  [u,v] = stokesSLPfmm(f1(:),f2(:),xsou(:),ysou(:),...
      xtar(:),ytar(:),0,iprec);

  stokesSLPtar = zeros(2*Ntar,ncol); % initialize
  for k = 1:ncol
    is = (k-1)*Ntar+1;
    ie = is + Ntar - 1;
    stokesSLPtar(1:Ntar,k) = u(is:ie);
    stokesSLPtar(Ntar+1:2*Ntar,k) = v(is:ie);
  end
  % Wrap the output of the FMM in the usual format
  % for the target points
end


end % exactStokesSLfmm




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stokesSLP,stokesSLPtar] = ...
    exactStokesSLfmmOld(o,geom,f,Xtar,K)
% [stokesSLP,stokeSLPtar] = exactStokesSLfmm(geom,f,Xtar,K) uses 
% the FMM to compute the single-layer potential due to all geoms
% except itself geom is a class of object bodies and f is the 
% density function NOT scaled by arclength term.  Xtar is a set of 
% points where the single-layer potential due to all geoms in index 
% set K needs to be evaulated
global fmms

fmms = fmms + 1;
% count the total number of calls to fmm

N = geom.N; % number of points per geom
nv = geom.nv; % number of geoms
X = geom.X; % geom positions
oc = curve;
[x,y] = oc.getXY(X); % seperate x and y coordinates

den = f.*[geom.sa;geom.sa]*2*pi/N;

if (nargin == 5)
  stokesSLP = [];
else
  [f1,f2] = oc.getXY(den);
  % need to multiply by arclength term.  Seperate it into
  % x and y coordinate

  [u,v] = stokesSLPfmm(f1(:),f2(:),x(:),y(:));
  stokesSLP = zeros(2*N,nv); % initialize
  for k = 1:nv
    is = (k-1)*N+1;
    ie = k*N;
    stokesSLP(1:N,k) = u(is:ie);
    stokesSLP(N+1:2*N,k) = v(is:ie);
  end
  % Wrap the output of the FMM into the usual 
  % [[x1;y1] [x2;y2] ...] format

  if (N <= 256)
    stokesSLP = stokesSLP - o.exactStokesSLDirect(geom,den);
  else
    for k = 1:nv
      [u,v] = stokesSLPfmm(f1(:,k),f2(:,k),x(:,k),y(:,k));
      stokesSLP(:,k) = stokesSLP(:,k) - [u;v];
    end
  end
  % Subtract potential due to each geom on its own.  Nothing
  % to change here for potential at Xtar
end

if nargin == 3
  stokesSLPtar = [];
else
  [x,y] = oc.getXY(X(:,K)); 
  % seperate x and y coordinates at geoms indexed by K
  [Ntar,ncol] = size(Xtar);
  Ntar = Ntar/2;
  x2 = Xtar(1:Ntar,:);
  x = [x(:);x2(:)];
  y2 = Xtar(Ntar+1:2*Ntar,:);
  y = [y(:);y2(:)];
  % Stack the x and y coordinates of the target points
  [f1,f2] = oc.getXY(den(:,K));
  % seperate x and y coordinates at geoms indexed by K
  f1 = [f1(:);zeros(Ntar*ncol,1)];
  f2 = [f2(:);zeros(Ntar*ncol,1)];
  % pad density function with zeros so that Xtar doesn't
  % affect the single-layer potential


  [u,v] = stokesSLPfmm(f1,f2,x,y);
  stokesSLPtar = zeros(2*Ntar,ncol); % initialize
  for k = 1:ncol
    is = N*numel(K) + (k-1)*Ntar+1;
    ie = is + Ntar - 1;
    stokesSLPtar(1:Ntar,k) = u(is:ie);
    stokesSLPtar(Ntar+1:2*Ntar,k) = v(is:ie);
  end
  % Wrap the output of the FMM in the usual format
  % for the target points
end

end % exactStokesSLfmm




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stokesDLP,stokesDLPtar] = ...
    exactStokesDLfmm(o,geom,f,Xtar,K)
% [stokesDLP,stokeDLPtar] = exactStokesDLfmm(geom,f,Xtar,K) uses the FMM
% to compute the double-layer potential due to all geoms except itself
% geom is a class of object bodies and f is the density function NOT
% scaled by arclength term.  Xtar is a set of points where the
% single-layer potential due to all geoms in index set K needs to be
% evaulated
global fmms

fmms = fmms + 1;
% count the total number of calls to fmm

N = geom.N; % number of points per geom
nv = geom.nv; % number of geoms
X = geom.X; % geom positions
oc = curve;
[x,y] = oc.getXY(X); % seperate x and y coordinates
[nx,ny] = oc.getXY(geom.normal);
% seperate the x and y coordinates of the normal vector

den = f.*[geom.sa;geom.sa]*2*pi/N;

if (nargin == 5)
  stokesDLP = [];
else
  stokesDLP = zeros(2*N,nv); % initialize

  [fx,fy] = oc.getXY(den);
  % need to multiply by arclength term.  Seperate it into
  % x and y coordinate
  dip1 = 0.25/pi*(fy - 1i*fx).*(nx + 1i*ny);
  dip2 = -1i*0.5/pi*(fx.*nx + fy.*ny);

  vel = stokesDLPfmm(dip1(:),dip2(:),x(:),y(:));
  u = -imag(vel);
  v = real(vel);
  for k = 1:nv
    is = (k-1)*N+1;
    ie = k*N;
    stokesDLP(1:N,k) = u(is:ie);
    stokesDLP(N+1:2*N,k) = v(is:ie);
  end
  % Wrap the output of the FMM into the usual 
  % [[x1;y1] [x2;y2] ...] format


  for k = 1:nv
    dip1 = 0.25/pi*(fy(:,k) - 1i*fx(:,k)).*(nx(:,k) + 1i*ny(:,k));
    dip2 = -1i*0.5/pi*(fx(:,k).*nx(:,k) + fy(:,k).*ny(:,k));

    vel = stokesDLPfmm(dip1,dip2,x,y);
    u = -imag(vel);
    v = real(vel);

    stokesDLP(:,k) = stokesDLP(:,k) - [u;v];
  end
  % Subtract potential due to each geom on its own.  Nothing
  % to change here for potential at Xtar
end

if nargin == 3
  stokesDLPtar = [];
else
  [x,y] = oc.getXY(X(:,K)); 
  % seperate x and y coordinates at geoms indexed by K
  [nx,ny] = oc.getXY(geom.normal(:,K));
  [Ntar,ncol] = size(Xtar);
  Ntar = Ntar/2;
  x2 = Xtar(1:Ntar,:);
  x = [x(:);x2(:)];
  y2 = Xtar(Ntar+1:2*Ntar,:);
  y = [y(:);y2(:)];
  % Stack the x and y coordinates of the target points
  [fx,fy] = oc.getXY(den(:,K));
  % seperate x and y coordinates at geoms indexed by K
  fx = [fx(:);zeros(Ntar*ncol,1)];
  fy = [fy(:);zeros(Ntar*ncol,1)];
  % pad density function with zeros so that Xtar doesn't
  % affect the double-layer potential
  nx = [nx(:);zeros(Ntar*ncol,1)];
  ny = [ny(:);zeros(Ntar*ncol,1)];
  % pad the normal vector with zeros so that Xtar doesn't
  % affect the double-layer potential

  dip1 = 0.25/pi*(fy - 1i*fx).*(nx + 1i*ny);
  dip2 = -1i*0.5/pi*(fx.*nx + fy.*ny);

  vel = stokesDLPfmm(dip1(:),dip2(:),x(:),y(:));
  u = -imag(vel);
  v = real(vel);

  stokesDLPtar = zeros(2*Ntar,ncol); % initialize
  for k = 1:ncol
    is = N*numel(K) + (k-1)*Ntar+1;
    ie = is + Ntar - 1;
    stokesDLPtar(1:Ntar,k) = u(is:ie);
    stokesDLPtar(Ntar+1:2*Ntar,k) = v(is:ie);
  end
  % Wrap the output of the FMM in the usual format
  % for the target points
end


end % exactStokesDLfmm



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [laplaceDLP,laplaceDLPtar] = ...
    exactLaplaceDL(o,geom,f,Xtar,K1)
% pot = exactLaplaceDL(geom,f,Xtar,K1) 
% computes the double-layer laplace potential due to f around all
% geoms except itself.  Also can pass a set of target points Xtar
% and a collection of geoms K1 and the double-layer potential due 
% to geoms in K1 will be evaluated at Xtar.
% Everything but Xtar is in the 2*N x nv format
% Xtar is in the 2*Ntar x ncol format

nv = geom.nv; % number of geoms
X = geom.X; % Vesicle positions
normal = geom.normal; % Outward normal
sa = geom.sa; % Jacobian

if nargin == 5
  Ntar = size(Xtar,1)/2;
  ncol = size(Xtar,2);
  laplaceDLPtar = zeros(2*Ntar,ncol);
else
  K1 = [];
  laplaceDLPtar = [];
  ncol = 0;
  Ntar = 0;
  % if nargin ~= 5, user does not need the layer potential at arbitrary
  % points
end

den = f.*[sa;sa]*2*pi/geom.N;
% multiply by arclength term

for k2 = 1:ncol % loop over columns of target points
  for j = 1:Ntar % loop over rows of target points
    dis2 = (X(1:geom.N,K1) - Xtar(j,k2)).^2 + ... 
        (X(geom.N+1:2*geom.N,K1) - Xtar(j+Ntar,k2)).^2;
    diffxy = [X(1:geom.N,K1) - Xtar(j,k2) ; ...
        X(geom.N+1:2*geom.N,K1) - Xtar(j+Ntar,k2)];
    % distance squared and difference of source and target location

    coeff = (diffxy(1:geom.N,:).*normal(1:geom.N,K1) + ...
      diffxy(geom.N+1:2*geom.N,:).*...
      normal(geom.N+1:2*geom.N,K1))./dis2;
    % this is the kernel of the double-layer potential for
    % laplace's equation
    val = coeff.*den(1:geom.N,K1);
    laplaceDLPtar(j,k2) = sum(val(:));
    val = coeff.*den(geom.N+1:2*geom.N,K1);
    laplaceDLPtar(j+Ntar,k2) = sum(val(:));
  end % j

end % k2
% Evaluate double-layer potential at arbitrary target points
laplaceDLPtar = 1/(2*pi)*laplaceDLPtar;
% 1/2/pi is the coefficient in front of the double-layer potential

laplaceDLP = zeros(geom.N,geom.nv); % Initialize to zero
% if we only have one geom, geoms of course can not collide
% Don't need to run this loop in this case
if (nargin == 3 && geom.nv > 1)
  for k1 = 1:geom.nv % geom of targets
    K = [(1:k1-1) (k1+1:geom.nv)];
    % Loop over all geoms except k1

    for j=1:geom.N
      dis2 = (X(1:geom.N,K) - X(j,k1)).^2 + ...
          (X(geom.N+1:2*geom.N,K) - X(j+geom.N,k1)).^2;
      diffxy = [X(1:geom.N,K) - X(j,k1); ...
          X(geom.N+1:2*geom.N,K) - X(j+geom.N,k1)];
      % distance squared and difference of source and target location

      coeff = (diffxy(1:geom.N,:).*normal(1:geom.N,K) + ...
        diffxy(geom.N+1:2*geom.N,:).*...
        normal(geom.N+1:2*geom.N,K))./dis2;
      % this is the kernel of the double-layer potential for
      % laplace's equation
      val = coeff.*den(1:geom.N,K);
      laplaceDLP(j,k1) = sum(val(:));
    end % j
  end % k1
  % Evaluate double-layer potential at geoms but oneself
  laplaceDLP = 1/(2*pi)*laplaceDLP;
  % 1/2/pi is the coefficient in front of the double-layer potential
end % nargin == 3

laplaceDLP = [laplaceDLP;laplaceDLP];

end % exactLaplaceDL



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = stokesDLmatrix(o,geom)
% D = stokesDLmatrix(geom), generate double-layer potential for Stokes
% geom is a data structure defined as in the bodies class D is
% (2N,2N,nv) array where N is the number of points per curve and nv is
% the number of curves in X 


N = geom.N; % Number of points per geom 
nv = geom.nv; % Number of geoms
X = geom.X; % Vesicle positions
oc = curve;
[x,y] = oc.getXY(X);

normal = geom.normal; % Normal vector
tang = geom.xt; % Tangent vector
sa = geom.sa; % Jacobian

D = zeros(2*N,2*N,nv);
for k=1:nv  % Loop over curves
  for j=1:N % Loop over targets
    rho2 = (x(j,k)-x(:,k)).^2 + (y(j,k)-y(:,k)).^2;
    rho2(j) = 1;
    % Set diagonal term to one to avoid dividing by zero

    coeff = ((x(j,k) - x(:,k)).*normal(1:N,k) + ...
        (y(j,k) - y(:,k)).*normal(N+1:2*N,k)).*...
        sa(:,k)./rho2.^2/pi;
    % part of kernel of the double-layer potential

    D(j,:,k) = 2*pi/N*[coeff.*(x(j,k) - x(:,k)).^2; ...
      coeff.*(x(j,k) - x(:,k)).*(y(j,k) - y(:,k))]';
    D(j+N,:,k) = 2*pi/N*[coeff.*(y(j,k) - y(:,k)).*(x(j,k)-x(:,k)); ...
      coeff.*(y(j,k) - y(:,k)).^2]';
    % Build double-layer potential matrix D

    rc = [j j+N];
    D(rc,rc,k) = -2*pi/N*sa(j,k)*geom.cur(j,k)/2/pi*...
      [tang(j,k)^2 tang(j,k)*tang(j+N,k);...
      tang(j+N,k)*tang(j,k) tang(j+N,k)^2];
    % Diagonal term requires the limiting value.  Otherwise, above formula
    % would divide by zero
  end % j
end % k


end % stokesDLmatrix


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N0 = stokesN0matrix(o,geom)
% N0 = stokesN0matrix(geom) generates the the integral 
% operator with kernel normal(x) \otimes normal(y) which removes
% the rank one defficiency of the double-layer potential.  Need
% this operator for solid walls

N = geom.N; % Number of points per geom 
nv = geom.nv; % Number of geoms
X = geom.X; % Vesicle positions
oc = curve;
[x,y] = oc.getXY(X);

normal = geom.normal; % Normal vector
sa = geom.sa; % Jacobian


N0 = zeros(2*N,2*N,nv);
for k = 1:1 % Loop over curves
  % Only want to form N0 for the outer boundary.
  % Rotlets and Stokeslets take care of null space for inner 
  % boundaries
  for j = 1:2*N % Loop over targets
    N0(j,:,k) = normal(j,k)*normal(:,k).*...
        [sa(:,k);sa(:,k)]*2*pi/N;
    % compute the columns of the modification to the double-layer
    % potential.
  end
end
% Use N0 if solving (-1/2 + DLP)\eta = f where f has no flux through
% the boundary.  By solving (-1/2 + DLP + N0)\eta = f, we guarantee
% that \eta also has no flux through the boundary.  This is not
% required, but it means we're enforcing one addition condition on eta
% which removes the rank one kernel.  DLP is the double-layer potential
% for stokes equation

end % stokesN0matrix


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function G = stokesSLmatrix(o,vesicle)
% G = stokesSLmatrix(vesicle) generate single layer potential for
% Stokes vesicle is a data structure defined as in the curve class
% G is (2N,2N,nv) array where N is the number of points per curve
% and nv is the number of curves in X 

oc = curve;
[x,y] = oc.getXY(vesicle.X);
% seperate x and y coordinates of vesicle

G = zeros(2*vesicle.N,2*vesicle.N,vesicle.nv);
% initalize single-layer potential to zero
for k=1:vesicle.nv  % Loop over curves
  for j=1:vesicle.N % Loop over targets
    ind = 1 + mod(j-1 + (0:vesicle.N-1),vesicle.N);
    xin = o.qp*x(ind,k);
    yin = o.qp*y(ind,k);

    rho = (xin-x(j,k)).^2 + (yin-y(j,k)).^2;
    br = 1./rho;

    logpart = -1/2 * o.qw .* log(rho);
    
    G(j,ind,k) = (logpart+...
        (o.qw.*( ( x(j,k) - xin).^2 .* br)))'*o.qp;
    G(j, vesicle.N+ind, k) = ...
        (o.qw.*(x(j,k)-xin).*(y(j,k)-yin) .* br)'*o.qp;
    G(j+vesicle.N, 1:vesicle.N ,k) = G(j,vesicle.N+1:2*vesicle.N,k);
    G(j+vesicle.N,vesicle.N+ind,k) = (logpart+...
        (o.qw.*( ( y(j,k) - yin).^2 .* br)))'*o.qp;
  end % j
  sak = repmat(vesicle.sa(:,k)',2*vesicle.N,2);
  % 2 copies of Jacobian
  
  G(:,:,k) = G(:,:,k).*sak;
  % multiply both components of G by the Jacobian
  % don't need to divide by 4*pi as it is stored in o.qw
end % k

end % stokesSLmatrix



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function qw = quadratureS(o,q,N);
% qw = quadratureS(q) generates the quadrature rules for a function
% with o.N/2 points and a logarithmic singularity at the origin.  q
% controls the accuracy.  All rules are from Alpert 1999.

[v,u,a] = o.getWeights(q,2);
[x,w] = o.regularWeights(a);

m = N/2;
n = m-2*a+1;
h = 1/m;

evalpots1 = v*h;
evalpots3 = 1-x*h;

ys = pi*[evalpots1;evalpots3];
wt = pi*h*[u;w];

wt = [wt; flipud(wt)]; ys = [ys; 2*pi-flipud(ys)];
of = fft1;
A = of.sinterpS(2*m, ys); 
h = pi/m; 
yt = [a*h:h:(m - a)*h]';
%regular points away from singularity
wt = [wt; h*ones(2*length(yt),1)]/4/pi;
%quadrature weights at regular points away from singularity
lyt = 2*length(yt); 
B = sparse(lyt, 2*m); 
pos = 1+[(a:m-a)'; (m+a:2*m-a)'];

for k = 1:lyt
  B(k, pos(k)) = 1;
end
A = [sparse(A); B];
qw = [wt, A];



end % quadratureS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v,u,a] = getWeights(o,q,ker)
% [v,u,a] = getWeights(q,ker) loads quadrature rules for 
% different types of singularities.  All rules come from 
% Bradley Alpert's papers.  We are interested in nodesLogx.dat;

switch ker
case 0
  xp = load('nodesr.dat');
  lenth = [2;4;8];
  par = [2;4;7];
case 1
  xp = load('nodes_sqrtx.dat');
  lenth = [4;8;16];
  par = [3;5;10];
case 2
  xp = load('nodesLog.dat');
  lenth = [3;7;15];
  par = [2;5;10];
end

switch q
case 4
  v = xp(1:lenth(1), 1);
  u = xp(1:lenth(1), 2);
  a = par(1);
case 8
  v = xp(1+lenth(1):lenth(1)+lenth(2),1);
  u = xp(1+lenth(1):lenth(1)+lenth(2),2);
  a = par(2);
case 16
  v = xp(1+lenth(2)+lenth(1):sum(lenth),1);
  u = xp(1+lenth(2)+lenth(1):sum(lenth),2);
  a = par(3);
end


end % getWeights

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,w] = regularWeights(o,a)
% [x,w] = regularWeights(a) gets quadrature rules
% for regular functions (Alpert's quadrature rules).  These are
% needed for integrating away from singularities

par = [2;3;4;5;7;10];
for i = 1:length(par)
  if par(i)==a
    key = i;
  end
end

lenth = [2;3;4;6;8;12];
xp = load('nodesRegular.dat');
if key==1
  starting = 1;
else
  starting = sum(lenth(1:key-1))+1; 
end

x = xp( starting:starting+lenth(key)-1,1);
w = xp(starting:starting+lenth(key)-1,2);

end % regularWeights

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
