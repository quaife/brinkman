classdef capsules < handle
% This class implements standard calculations that need to
% be done to a vesicle, solid wall, or a collection of arbitrary
% target points (such as tracers or pressure/stress targets)
% Given a vesicle, the main tasks that can be performed are
% computing the required derivatives (bending, tension, surface
% divergence), the traction jump, the pressure and stress, 
% and constructing structures required for near-singluar
% integration

properties
N;          % number of points per vesicle
nv;         % number of vesicles
X;          % positions of vesicles
sig;        % tension of vesicles
u;          % velocity field of vesicles
kappa;      % bending modulus
viscCont;   % viscosity contrast
xt;         % tangent unit vector
sa;         % Jacobian
isa;        % inverse of Jacobian
length;     % minimum length over all vesicles
cur;        % curvature
center;     % center of the point required for stokeslets
            % and rotlets
IK;         % index of Fourier modes for fft and ifft
            % that are needed repetatively

end %properties

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = capsules(X,sigma,u,kappa,viscCont)
% capsules(X,sigma,u,kappa,viscCont) sets parameters and options for
% the class; no computation takes place here.  
%
% sigma and u are not needed and typically unknown, so just set them to
% an empty array. (Sorry!)

o.N = size(X,1)/2;              % points per vesicle
o.nv = size(X,2);               % number of vesicles
o.X = X;                        % position of vesicle
oc = curve;
[o.sa,o.xt,o.cur] = oc.diffProp(o.X);
o.isa = 1./o.sa;
% Jacobian, tangent, and curvature
o.sig = sigma;          % Tension of vesicle
o.u = u;                % Velocity of vesicle
o.kappa = kappa;        % Bending modulus
o.viscCont = viscCont;  % Viscosity contrast
o.center = [mean(X(1:o.N,:));mean(X(o.N+1:2*o.N,:))];
% center of vesicle.  Required for center of rotlets and
% stokeslets in confined flows

[~,~,len] = oc.geomProp(X);
o.length = min(len);
% minimum arclength needed for near-singular integration
o.IK = fft1.modes(o.N,o.nv);
% ordering of the fourier modes.  It is faster to compute once here and
% pass it around to the fft differentitation routine

end % capsules

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ben,ten,adh] = computeEnergies(vesicle,W0,d0)

N = vesicle.N;

ben = vesicle.kappa.*sum(vesicle.cur.^2.*vesicle.sa)*2*pi/N;
ten = sum(vesicle.sig.*vesicle.sa)*2*pi/N;

adh = zeros(vesicle.nv,1);
adhPotential = zeros(N,1);
for j = 1:N
  dist2 = (vesicle.X(j,1) - vesicle.X(1:N,2)).^2 + ...
          (vesicle.X(j+N,1) - vesicle.X(N+1:2*N,2)).^2;
  adhPotential(j) = W0*sum((d0./dist2).^2.*((d0./dist2).^2 - 2).*...
    vesicle.sa(:,2))*2*pi/N;
end
adh(1) = sum(adhPotential.*vesicle.sa(:,1))*2*pi/N;
adh(2) = sum(adhPotential.*vesicle.sa(:,1))*2*pi/N;

end % computeEnergies


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = tracJump(o,f,sigma)
% tracJump(f,sigma) computes the traction jump where the derivatives
% are taken with respect to a linear combiation of previous time steps
% which is stored in object o Xm is 2*N x nv and sigma is N x nv
  
f = o.bendingTerm(f) + o.tensionTerm(sigma);

end % tracJump

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ben = bendingTerm(o,f)
% ben = bendingTerm(f) computes the term due to bending
% -kappa*fourth-order derivative

ben = -[curve.arcDeriv(f(1:o.N,:),4,o.isa,o.IK);...
  curve.arcDeriv(f(o.N+1:2*o.N,:),4,o.isa,o.IK)]*diag(o.kappa);

end % bendingTerm
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ten = tensionTerm(o,sig)
% ten = tensionTerm(o,sig) computes the term due to tension (\sigma *
% x_{s})_{s}

ten = [curve.arcDeriv(sig.*o.xt(1:o.N,:),1,o.isa,o.IK);...
    curve.arcDeriv(sig.*o.xt(o.N+1:2*o.N,:),1,o.isa,o.IK)];

end % tensionTerm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = surfaceDiv(o,f)
% divf = surfaceDiv(f) computes the surface divergence of f with respect
% to the vesicle stored in object o.  f has size N x nv

oc = curve; 
[fx,fy] = oc.getXY(f);

[tangx,tangy] = oc.getXY(o.xt);
f = curve.arcDeriv(fx,1,o.isa,o.IK).*tangx + ...
  curve.arcDeriv(fy,1,o.isa,o.IK).*tangy;

end % surfaceDiv

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function adh = adhesionTerm(o,Wconst,d0)
% adh = adhsionTerm(Wconst,d0) computes the term due to adhesion.
% Wconst is the strength of the adhesion force and d0 corresponds to
% the distance where the adhesive force is the strongest

adh = zeros(2*o.N,o.nv);
for k = 1:o.nv
  notk = [(1:k-1) (k+1:o.nv)];
  adhesive_x = o.X(1:o.N,notk);
  adhesive_y = o.X(o.N+1:2*o.N,notk);
  ds = o.sa(:,notk)*2*pi/o.N;
  
  for j = 1:o.N
    dist2 = (o.X(j,k) - adhesive_x).^2 + ...
        (o.X(j+o.N,k) - adhesive_y).^2;
%    Dw = -Wconst * 4*d0^2.*(d0^2 - dist2)./dist2.^(2.5);
    Dw = -Wconst*4*d0^2.*(d0^2 - dist2)./dist2.^3;

    adh(j,k) = sum(sum(-Dw .* (o.X(j,k) - adhesive_x).*ds));
    adh(j+o.N,k) = sum(sum(-Dw .* (o.X(j+o.N,k) - adhesive_y).*ds));
  end

end

end % adhesionTerm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ben,Ten,Div] = computeDerivs(o)
% [Ben,Ten,Div] = computeDerivs computes the matricies that takes a
% periodic function and maps it to the fourth derivative, tension, and
% surface divergence all with respect to arclength.  Everything in this
% routine is matrix free at the expense of having repmat calls

persistent D1 

% persistent variable for matrix that does Fourier differentiation (not with respect to arclength)
if isempty(D1);
  D1 = fft1.fourierDiff(o.N);
end
% Build  Fourier differentiation matrix if it is not constructed yet.
% It can be used throught the whole simulation as it only depends on
% the spatial discretization N

Ben = zeros(2*o.N,2*o.N,o.nv);
Ten = zeros(2*o.N,o.N,o.nv);
Div = zeros(o.N,2*o.N,o.nv);

for k = 1:o.nv
  % compute single arclength derivative matrix
  isa = o.isa(:,k);
  arcDeriv = isa(:,ones(o.N,1)).*D1;
  % This line is equivalent to repmat(o.isa(:,k),1,o.N).*D1 but much
  % faster.

  D4 = arcDeriv*arcDeriv; D4 = D4*D4;
  Ben(:,:,k) = [D4 zeros(o.N); zeros(o.N) D4];

  Ten(:,:,k) = [arcDeriv*diag(o.xt(1:o.N,k));...
               arcDeriv*diag(o.xt(o.N+1:end,k))];

  Div(:,:,k) = [diag(o.xt(1:o.N,k))*arcDeriv ...
                diag(o.xt(o.N+1:end,k))*arcDeriv];
end

Ben = real(Ben);
Ten = real(Ten);
Div = real(Div);
% Imaginary part should be 0 since we are preforming a real operation

end % computeDerivs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sigma,eta,RS,vesVel,iter] = computeSigAndEta(...
    vesicle,tt,walls)
% computes the tension, density function (if flow is confined), and
% velocity of the vesicle for the vesicle configuration given by
% vesicle.X and solid wall configuration given by walls.X

N = vesicle.N; % points per vesicle
nv = vesicle.nv; % number of vesicles
if tt.confined
  Nbd = walls.N; % points per solid wall
  nvbd = walls.nv; % number of solid walls
else
  Nbd = 0;
  nvbd = 0;
end

op = tt.op;
[Ben,Ten,Div] = vesicle.computeDerivs;

% compute self bending, tension, and divegence terms
tt.Galpert = op.stokesSLmatrix(vesicle);
% single-layer potential
tt.D = op.stokesDLmatrix(vesicle);
% double-layer potential

if tt.near
  if tt.confined
    [tt.NearV2V,tt.NearV2W] = vesicle.getZone(walls,3);
    [tt.NearW2W,tt.NearW2V] = walls.getZone(vesicle,3);
  else
    tt.NearV2V = vesicle.getZone([],1);
    tt.NearW2V = [];
    tt.NearV2W = [];
    tt.NearW2W = [];
  end
  % near-singular integration structures
else
  tt.NearV2V = [];
  tt.NearW2V = [];
  tt.NearV2W = [];
  tt.NearW2W = [];
  % empty since we are not doing near-singular integration
end

if tt.confined
  rhs = zeros(2*N*nv + N*nv + 2*Nbd*nvbd + 3*(nvbd-1),1);
else
  rhs = zeros(2*N*nv + N*nv,1);
end
% initalize the right-hand side for [vesicle velocity;tension] stacked
% over each vesicle followed by the density function followed by the
% rotlets and stokeslets

if ~tt.fmm
  kernel = @op.exactStokesSL;
  kernelDirect = @op.exactStokesSL;
else
  kernel = @op.exactStokesSLfmm;
  kernelDirect = @op.exactStokesSL;
end
% kernel for single-layer potential.  Only difference is if the FMM is
% used or not

f = vesicle.tracJump(vesicle.X,zeros(N,nv));
% traction jump due only to the position

if ~tt.near
  Fslp = kernel(vesicle,f);
  if tt.confined
    [~,FSLPwall] = kernel(vesicle,f,[],walls.X,(1:nv));
  end
else
  SLP = @(X) op.exactStokesSLdiag(vesicle,tt.Galpert,X);
  Fslp = op.nearSingInt(vesicle,f,SLP,...
      tt.NearV2V,kernel,kernelDirect,vesicle,true,false);
  
  if tt.confined
    FSLPwall = op.nearSingInt(vesicle,f,SLP,...
      tt.NearV2W,kernel,kernelDirect,walls,false,false);
    % Evaluate single-layer potential due to all vesicles on
    % the solid walls WITH near-singular integration
  end
end
% single-layer potential due to all vesicles except the diagonal one

if ~tt.confined
  Ffar = tt.farField(vesicle.X);
else
  Ffar = zeros(2*N,nv);
  U = tt.farField(walls.X);
  % no slip boundary condition for velocity on solid walls
end
% If doing unbounded flow, add in the background velocity

velBen = Fslp + op.exactStokesSLdiag(vesicle,tt.Galpert,f) + Ffar;


for k = 1:nv
  istart = (k-1)*3*N + 1;
  iend = istart + 2*N - 1;
  rhs(istart:iend) = velBen(:,k);
end

for k = 1:nvbd
  istart = 3*nv*N+1 + (k-1)*2*Nbd;
  iend = istart + 2*Nbd - 1;
  rhs(istart:iend) = U(:,k) - FSLPwall(:,k);
end

alpha = (1+vesicle.viscCont)/2;
for k = 1:nv
  if any(vesicle.viscCont ~= 1)
    [tt.bdiagVes.L(:,:,k),tt.bdiagVes.U(:,:,k)] = lu( ...
      [(alpha(k)*eye(2*N) - tt.D(:,:,k)) ...
          -tt.Galpert(:,:,k)*Ten(:,:,k); ...
      Div(:,:,k) zeros(N)]);
  else
    [tt.bdiagVes.L(:,:,k),tt.bdiagVes.U(:,:,k)] = lu(...
      [eye(2*N) -tt.Galpert(:,:,k)*Ten(:,:,k); ...
      Div(:,:,k) zeros(N)]);
  end
end

% build vesicle part of the block-diagonal preconditioner
% GMRES

if 1
tic
[sigDen,F,R,I] = gmres(@(X) tt.sigDenMatVec(X,vesicle,walls),rhs,...
    [],tt.gmresTol,min(tt.gmresMaxIter,N*nv+2*Nbd*nvbd + 3*(nvbd-1)),...
    @tt.preconditionerBD);
toc
%pause
% solve for tension and density function with block-diagonal 
% preconditioning
iter = I(2);
end

% Solve with LU DECOMPOSITION
if 0
sigDen = tt.preconditionerBD(rhs);
iflag = 0;
R = 0;
I = [1;1];
resvec = [0];
iter = I(2);
end

vesVel = zeros(2*N,nv);
sigma = zeros(N,nv); % tension
eta = zeros(2*Nbd,nvbd); % solid wall density
RS = zeros(3,nvbd); % rotlets and stokeslets

for k = 1:nv
  vesVel(:,k) = sigDen((k-1)*3*N+1:(3*k-1)*N);
  sigma(:,k) = sigDen((3*k-1)*N+1:3*k*N);
end

for k = 1:nvbd
  eta(:,k) = sigDen(3*nv*N+(k-1)*2*Nbd+1:3*nv*N+k*2*Nbd);
end
for k = 2:nvbd
  istart = 3*nv*N + 2*nvbd*Nbd + 3*(k-2) + 1;
  iend = istart + 2;
  RS(:,k) = sigDen(istart:iend);
end
% unstack the tension and density function from the GMRES solver

end % computeSigAndEta

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eta,RS,iter] = computeEta(walls,tt)
% [eta,RS,iter] = computeEta(walls,tt) computes wall densities and RS when
% there is no vesicle but we want to compute the tracers. This is called by
% tstep.m/noVesTracers. It uses the the exact inverse of
% the double-layer potential for stokes flow in a bounded domain which is
% computed by tstep.m/wallsPrecond.

Nbd = walls.N; % points per solid wall
nvbd = walls.nv; % number of solid walls

rhs = zeros(2*Nbd*nvbd + 3*(nvbd-1),1);

U = tt.farField(walls.X);
% no slip boundary condition for velocity on solid walls


for k = 1:nvbd
  istart = 1 + (k-1)*2*Nbd;
  iend = k*2*Nbd;
  rhs(istart:iend) = U(:,k);
end

etaRS = tt.bdiagWall*rhs;
I = [1;1];
iter = I(2);

eta = zeros(2*Nbd,nvbd); % solid wall density
RS = zeros(3,nvbd); % rotlets and stokeslets

for k = 1:nvbd
  eta(:,k) = etaRS((k-1)*2*Nbd+1:k*2*Nbd);
end
for k = 2:nvbd
  istart = 2*nvbd*Nbd + 3*(k-2) + 1;
  iend = istart + 2;
  RS(:,k) = etaRS(istart:iend);
end
% unstack the density and RS function from the GMRES solver

end % computeSigAndEta


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NearSelf,NearOther] = getZone(vesicle1,vesicle2,relate)
% [NearSelf,NearOther] = getZone(vesicle1,vesicle2,relate) constructs
% each vesicle, index of the closest point, nearest point on a local
% interapolant, and argument of that nearest point.  vesicle1 contains
% the source points (which are also target points) and vesicle2 
% contains additional target points.  The 
% values of relate corresond to
% relate == 1  => only require NearSelf  (ie. vesicle to vesicle)
% relate == 2  => only require NearOther (ie. vesicle to wall)
% relate == 3  => require both NearSelf and NearOther
% THIS ROUTINE HAS A LOOP OVER THE TARGET POINTS WHICH SHOULD BE REMOVED
NearSelf = [];
NearOther = [];

N1 = vesicle1.N; % number of source/target points
nv1 = vesicle1.nv; % number of source/target vesicles
X1 = vesicle1.X; % source and target points
oc = curve;
[xsou,ysou] = oc.getXY(X1); 
% separate targets into x and y coordinates

h = vesicle1.length/N1; 
% smallest arclength over all vesicles
ptsperbox = 10; 
% Estimate for number of points per box.  This simply sets the 
% number of uniformly refined boxes we take.  Estimate is not very
% accurate.  What ptsperbox represents is the total number of points
% that could be put in each two-dimensional bin where no two are
% less than distance h from one another.  However, our points live
% on curves and thus will not fill up an entire bin

H = sqrt(ptsperbox)*h;
xmin = min(min(xsou));
xmax = max(max(xsou));
xmin = xmin - H;
xmax = xmax + H;
ymin = min(min(ysou));
ymax = max(max(ysou));
ymin = ymin - H;
ymax = ymax + H;
% Add a buffer around the points so that it is easier to
% work with vesicle2

Nx = ceil((xmax - xmin)/H);
Ny = ceil((ymax - ymin)/H);
% Find bounds for box that contains all points and add a buffer 
% so that all points are guaranteed to be in the box

Nbins = Nx * Ny; % Total number of bins

ii = ceil((xsou - xmin)/H);
jj = ceil((ysou - ymin)/H);
% Index in x and y direction of the box containing each point
bin = (jj-1)*Nx + ii;
% Find bin of each point using lexiographic ordering (x then y)


%figure(2);
%clf; hold on
%plot(xsou,ysou,'k.')
%axis equal
%axis([xmin xmin+Nx*H ymin ymin+Ny*H])
%set(gca,'xtick',linspace(xmin,xmin+Nx*H,Nx+1))
%set(gca,'ytick',linspace(ymin,ymin+Ny*H,Ny+1))
%grid on
%set(gca,'xticklabel',[])
%set(gca,'yticklabel',[])
%figure(1)
%pause
% DEBUG: This does a simple plot of the points with a grid that 
% aligns with the boundary of the boxes

%whos
%disp([Nbins nv1])
%disp([xmin xmax ymin ymax])
%disp([h H])
%pause
fpt = zeros(Nbins,nv1);
lpt = zeros(Nbins,nv1);
% allocate space for storing first and last points
[binsort,permute] = sort(bin);
% build permute.  Need binsort to find first and last points
% in each box

for k = 1:nv1 % Loop over vesicles
  for j = 1:N1 % Loop over bins
    ibox = binsort(j,k);
    if (fpt(ibox,k) == 0)
      fpt(ibox,k) = j;
      lpt(ibox,k) = 1;
    else
      lpt(ibox,k) = lpt(ibox,k) + 1;
    end
  end
  lpt(:,k) = fpt(:,k) + lpt(:,k) - 1;
end
% Construct first and last point in each box corresponding
% to each vesicle.  The order is based on permute.  For example,
% permute(fpt(ibox,k)),...,permute(lpt(ibox,k)) is the set of 
% all points from vesicle k contained in box ibox

neigh = zeros(Nbins,9);

%Do corners first
neigh(1,1:4) = [1 2 Nx+1 Nx+2]; 
% bottom left corner
neigh(Nx,1:4) = [Nx Nx-1 2*Nx 2*Nx-1]; 
% bottom right corner
neigh(Nbins-Nx+1,1:4) = [Nbins-Nx+1 Nbins-Nx+2 ...
    Nbins-2*Nx+1 Nbins-2*Nx+2];
% top left corner
neigh(Nbins,1:4) = [Nbins Nbins-1 Nbins-Nx Nbins-Nx-1]; 
% top right corner

for j = 2:Nx-1
  neigh(j,1:6) = j + [-1 0 1 Nx-1 Nx Nx+1];
end
% neighbors of bottom row

for j = Nbins-Nx+2:Nbins-1
  neigh(j,1:6) = j + [-1 0 1 -Nx-1 -Nx -Nx+1];
end
% neighbors of top row

for j=Nx+1:Nx:Nbins-2*Nx+1
  neigh(j,1:6) = j + [-Nx -Nx+1 0 1 Nx Nx+1];
end
% neighbors of left column

for j=2*Nx:Nx:Nbins-Nx
  neigh(j,1:6) = j + [-Nx-1 -Nx -1 0 Nx-1 Nx];
end
% neighbors of right column

J = (Nx + 1:Nbins - Nx);
J = J(mod(J-1,Nx)~=0);
J = J(mod(J,Nx)~=0);
% J is the index of boxes that are not on the boundary
for j=J
  neigh(j,:) = j + [-Nx-1 -Nx -Nx+1 -1 0 1 Nx-1 Nx Nx+1];
end
% neighbors of interior points
% TREE STRUCTURE IS COMPLETE


if (relate == 1 || relate == 3)
%  distSS = -ones(N1,nv1,nv1); 
%  % dist(n,k,j) is the distance of point n on vesicle k to
%  % vesicle j
%  zoneSS = -ones(N1,nv1,nv1); 
%  % near or far zone
%  nearestSS = -ones(2*N1,nv1,nv1); 
%  % nearest point using local interpolant
%  icpSS = -ones(N1,nv1,nv1); 
%  % index of closest discretization point
%  argnearSS = -ones(N1,nv1,nv1); 
%  % argument in [0,1] of local interpolant
  for k = 1:nv1
    distSS{k} = spalloc(N1,nv1,0);
    % dist(n,k,j) is the distance of point n on vesicle k to
    zoneSS{k} = spalloc(N1,nv1,0);
    % near or far zone
    nearestSS{k} = spalloc(2*N1,nv1,0);
    % nearest point using local interpolant
    icpSS{k} = spalloc(N1,nv1,0);
    % index of closest discretization point
    argnearSS{k} = spalloc(N1,nv1,0);
    % argument in [0,1] of local interpolant
  end
  % New way of representing near-singular integration structure so that
  % we can use sparse matricies.


  % begin classifying points where we are considering 
  % vesicle to vesicle relationships
  for k = 1:nv1
    boxes = unique(bin(:,k));
    % Find all boxes containing points of vesicle k
    boxes = neigh(boxes,:);
    % Look at all neighbors of boxes containing vesicle k
    boxes = unique(boxes(:));
    % Remove repetition
    boxes = boxes(boxes~=0);
    % Delete non-existent boxes that came up because of neigh

    K = [(1:k-1) (k+1:nv1)];
    for k2 = K
      istart = fpt(boxes,k2);
      iend = lpt(boxes,k2);
      istart = istart(istart ~= 0);
      iend = iend(iend ~= -1);
      % Find index of all points in neighboring boxes of vesicle k
      % that are in vesicle k2
      
      neighpts = zeros(sum(iend-istart+1),1);
      
      % Allocate space to assign possible near points
      is = 1;
      for j=1:numel(istart)
        ie = is + iend(j) - istart(j);
        neighpts(is:ie) = permute(istart(j):iend(j),k2);
        is = ie + 1;
      end
      % neighpts contains all points on vesicle k2 that are in 
      % neighboring boxes to vesicle k

      neighpts = sort(neighpts);
      % sorting should help speedup as we won't be jumping around
      % through different boxes

      n = 0;
      for i=1:numel(neighpts)
        ipt = neighpts(i);
        ibox = bin(ipt,k2);
        % box containing ipt on vesicle k2
        if (ibox ~= n)
          n = ibox;
          % Check if we need to move to a new box
          neighbors = neigh(ibox,:);
          % neighbors of this box
          neighbors = neighbors(neighbors~=0);
          % Remove non-existent neighbors
          istart = fpt(neighbors,k);
          iend = lpt(neighbors,k);
          istart = istart(istart ~= 0);
          iend = iend(iend ~= -1);
          % Find points on vesicle k in neighboring boxes
          neighpts2 = zeros(sum(iend-istart+1),1);
          is = 1;
          for j=1:numel(istart)
            ie = is + iend(j) - istart(j);
            neighpts2(is:ie) = permute(istart(j):iend(j),k);
            is = ie + 1;
          end
          % neighpts2 contains all points on vesicle k that 
          % are in neighboring box of ibox 
        end % decide if we need to switch boxes

        [d0,d0loc] = min((xsou(ipt,k2) - xsou(:,k)).^2 + ...
            (ysou(ipt,k2) - ysou(:,k)).^2);
        % Find minimum distance between ipt on vesicle k2 to
        % possible closest points on vesicle k
        d0 = sqrt(d0);
        % Save on not taking the square root on a vector but instead
        % on a single real number

        icpSS{k}(ipt,k2) = d0loc;
        if (d0 < 2*h);
          [distSS{k}(ipt,k2),nearestx,nearesty,argnearSS{k}(ipt,k2)] = ...
              vesicle1.closestPnt([xsou;ysou],xsou(ipt,k2),...
              ysou(ipt,k2),k,icpSS{k}(ipt,k2));
          nearestSS{k}(ipt,k2) = nearestx;
          nearestSS{k}(ipt+N1,k2) = nearesty;
          % Find closest point along a local interpolant using
          % Newton's method.

          if (distSS{k}(ipt,k2) < h)
            zoneSS{k}(ipt,k2) = 1;
          end
          % Point ipt of vesicle k2 is in the near zone of
          % vesicle k
        end


      end % ipt

    end % k2

  end % k

  NearSelf.dist = distSS;
  NearSelf.zone = zoneSS;
  NearSelf.nearest = nearestSS;
  NearSelf.icp = icpSS;
  NearSelf.argnear = argnearSS;
  % Store everything in the structure NearSelf.  This way it is 
  % much cleaner to pass everything around
end % relate == 1 || relate == 3

% Bin target points with respect to the source points
if (relate == 2 || relate == 3)
  N2 = vesicle2.N; % number of additional targets
  nv2 = vesicle2.nv; % number of additional vesicles
  X2 = vesicle2.X; % additional target points
  [xtar,ytar] = oc.getXY(X2);
  % separate additional target points into x and y coordinates
%  figure(2); clf
%  plot(xtar,ytar,'r.')
%  hold on
%  plot(xsou,ysou,'k.')
%  pause
% DEBUG: FOR SEEING TARGET AND SOURCE POINTS IN THE TREE STRUCTURE
% WHICH CAN BE PLOTTED ABOVE

  for k = 1:nv1
    distST{k} = spalloc(N1,nv2,0);
    % dist(n,k,j) is the distance of point n on vesicle k to
    zoneST{k} = spalloc(N1,nv2,0);
    % near or far zone
    nearestST{k} = spalloc(2*N1,nv2,0);
    % nearest point using local interpolant
    icpST{k} = spalloc(N1,nv2,0);
    % index of closest discretization point
    argnearST{k} = spalloc(N1,nv2,0);
    % argument in [0,1] of local interpolant
  end
  % Represent near-singular integration structure using sparse matricies

  itar = ceil((xtar - xmin)/H);
  jtar = ceil((ytar - ymin)/H);
  [indx,indy] = find((itar >= 1) & (itar <= Nx) & ...
      (jtar >= 1) & (jtar <= Ny));
  % Only have to consider xx(ind),yy(ind) since all other points
  % are not contained in the box [xmin xmax] x [ymin ymax]

  for k = 1:nv1 % loop over sources
    for nind = 1:numel(indx) 
      % loop over points that are not outside the box that surrounds
      % all target points with a sufficiently large buffer
      ii = indx(nind);
      jj = indy(nind);
      binTar = (jtar(ii,jj)-1)*Nx + itar(ii,jj);
      boxesTar = neigh(binTar,:);
      boxesTar = boxesTar(boxesTar~=0);
      istart = fpt(boxesTar,k);
      iend  = lpt(boxesTar,k);
      istart = istart(istart ~= 0);
      iend   = iend(iend ~= -1);
      
      neighpts = zeros(sum(iend-istart+1),1);
  
      % Allocate space to assign possible near points
      if numel(neighpts) > 0
        % it is possible of the neighboring boxes to contain
        % no points.
        is = 1;
        
        for j = 1:numel(istart)
          ie = is + iend(j) - istart(j);
          neighpts(is:ie) = permute(istart(j):iend(j),k);
          is = ie + 1;
        end
        % Set of potentially nearest points to 
        % (xtar(jj),ytar(jj))
        
        [d0,d0loc] = min((xtar(ii,jj) - xsou(neighpts,k)).^2 + ...
          (ytar(ii,jj) - ysou(neighpts,k)).^2);
        % find closest point and distance between (xtar(jj),ytar(jj))
        % and vesicle k.  Only need to look at points in neighboring
        % boxes
        
        d0 = d0.^0.5;
        icpST{k}(ii,jj) = neighpts(d0loc);

        if d0 < 2*h
          [distST{k}(ii,jj),nearestx,nearesty,argnearST{k}(ii,jj)] = ...
            vesicle1.closestPnt([xsou;ysou],xtar(ii,jj),...
            ytar(ii,jj),k,icpST{k}(ii,jj));
          nearestST{k}(ii,jj) = nearestx;
          nearestST{k}(ii+N2,jj) = nearesty;
          
          % DEBUG: CHECK THAT NEWTON'S METHOD HAS DONE A GOOD JOB
          % CONVERGING TO THE NEAREST POINT
          % compute distance and nearest point between 
          % (xtar(ii,jj),ytar(ii,jj)) and vesicle k
          if distST{k}(ii,jj) < h
            zoneST{k}(ii,jj) = 1;
            % (xtar(ii,jj),ytar(ii,jj)) is in the near zone of vesicle k
          end
        end % d0 < 2*h
      end % numel(neighpts) > 0

    end % ii and jj

  end % k

  NearOther.dist = distST;
  NearOther.zone = zoneST;
  NearOther.nearest = nearestST;
  NearOther.icp = icpST;
  NearOther.argnear = argnearST;
  % store near-singluar integration requirements in structure NearOther

end % relate == 2 || relate == 3
end % getZone


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dist,nearestx,nearesty,theta] = ...
    closestPnt(o,X,xtar,ytar,k,icp)
% [dist,nearestx,nearesty,theta] = closestPnt(X,xtar,ytar,k,icp)
% computes the closest point on vesicle k to (xtar,ytar)
% using a Lagrange interpolant.  icp is the index of the closest
% point on the discrete mesh which is used as an initial guess

N = size(X,1)/2; % Number of points per vesicle
A = poten.lagrangeInterp;
interpOrder = size(A,1);
% need interpolation matrix and its size

p = ceil((interpOrder+1)/2);
% Accommodate for either an even or odd number of interpolation points
pn = mod((icp-p+1:icp-p+interpOrder)' - 1,N) + 1;
% band of points around icp.  The -1,+1 combination sets index
% 0 to N as required by the code

px = A*X(pn,k); % polynomial interpolant of x-coordinate
py = A*X(pn+N,k); % polynomial interpolant of y-coordinate
Dpx = px(1:end-1).*(interpOrder-1:-1:1)';
Dpy = py(1:end-1).*(interpOrder-1:-1:1)';
D2px = Dpx(1:end-1).*(interpOrder-2:-1:1)';
D2py = Dpy(1:end-1).*(interpOrder-2:-1:1)';
% To do Newton's method, need two derivatives

theta = 1/2;
% midpoint is a good initial guess
for newton = 1:1
  zx = filter(1,[1 -theta],px);
  zx = zx(end);
  zy = filter(1,[1 -theta],py);
  zy = zy(end);
  Dzx = filter(1,[1 -theta],Dpx);
  Dzx = Dzx(end);
  Dzy = filter(1,[1 -theta],Dpy);
  Dzy = Dzy(end);
  D2zx = filter(1,[1 -theta],D2px);
  D2zx = D2zx(end);
  D2zy = filter(1,[1 -theta],D2py);
  D2zy = D2zy(end);
  % Using filter is the same as polyval, but it is much
  % faster when only requiring a single polyval such as here.

  newtonNum = (zx-xtar)*Dzx + (zy-ytar)*Dzy;
  % numerator of Newton's method
  newtonDen = (zx-xtar)*D2zx + (zy-ytar)*D2zy + ...
      Dzx^2 + Dzy^2;
  % denominator of Newton's method
  theta = theta - newtonNum/newtonDen;
  % one step of Newton's method
end
% Do a few (no more than 3) Newton iterations

nearestx = filter(1,[1,-theta],px);
nearestx = nearestx(end);
nearesty = filter(1,[1,-theta],py);
nearesty = nearesty(end);
dist = sqrt((nearestx - xtar)^2 + (nearesty - ytar)^2);
% Compute nearest point and its distance from the target point

end % closestPnt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [icollisionVes,icollisionWall] = ...
    collision(vesicle,walls,NearV2V,NearV2W,fmm,near,op)
% [icollisionVes,icollisionWall] = ...
%    collision(vesicle,walls,NearV2V,NearV2W,fmm,near)
% uses double-layer potential for Laplace's equation to determing if
% vesicles have crossed (stored in icollisionVes) or if solid wall 
% point is inside a vesicle (icollisionWall).
% Can be done with or without the fmm.  vesicle stores all the
% vesicles, walls stores the solid walls (if they exist), NearV2V
% and NearV2W are near-singular integration structures, fmm and near
% are flags for using the FMM and near-singular integration

oc = curve;
[x,y] = oc.getXY(vesicle.X);

f = [ones(vesicle.N,vesicle.nv);zeros(vesicle.N,vesicle.nv)];
% Density function is constant.  Pad second half of it with zero


if ~fmm
  kernel = @op.exactLaplaceDL;
else
  kernel = @op.exactLaplaceDLfmm;
end
% kernel for laplace's double layer potential.  Only difference
% is if FMM is used or not

if near
  DLP = @(X) zeros(2*size(X,1),size(X,2));
  % can cheat here because we know that the double-layer
  % potential applied to our function f will always be 0
  % This won't work if we are considering density functions
  % that are not one everywhere

  Fdlp = op.nearSingInt(vesicle,f,DLP,...
          NearV2V,kernel,kernel,vesicle,true,false);
  % Apply the double-layer laplace potential with constant boundary
  % condition.  Skips self-vesicle term.  This means that if the vesicles have
  % crossed, the Fdlp will achieve a value of 1.  If they have not crossed, it
  % will be 0
  if ~isempty(walls)
    FDLPwall = op.nearSingInt(vesicle,f,DLP,...
            NearV2W,kernel,kernel,walls,false,false);
  end
  % Evaulate the laplace double-layer potential due to all vesicles with the
  % solid walls acting as the targets.  Will be zero at all the targets if
  % none of the vesicles have crossed the solid walls
else
  Fdlp = kernel(vesicle,f,[]);
  if ~isempty(walls)
    [~,FDLPwall] = kernel(vesicle,f,[],walls.X,(1:vesicle.nv));
    % vesicles are sources and solid walls are targets.  If vesicles have not
    % crossed a solid wall, values is 0
  end
end
Fdlp = Fdlp(1:vesicle.N,:);
% this layer-potential is not vector valued, so we throw away the second half.

bufferVes = 1e-4;
bufferWall = 1e-4;
% can't set buffer too big because near singular integration does not
% assign a value of 1 when near points cross.  This is because I did not
% swtich the direction of the normal for this case.  So, the lagrange
% interpolation points may look something like [0 1 1 1 1 ...] instead
% of [1 1 1 1 1 ...].  The local interpolant is affected by this and
% a value somewhere between 0 and 1 will be returned
icollisionVes = any(abs(Fdlp(:)) > bufferVes);

if ~isempty(walls)
  FDLPwall = FDLPwall(1:walls.N,:);
  icollisionWall = any(abs(FDLPwall(:)) > bufferWall);
else
  icollisionWall = false;
end

end % collision

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InOut = sortPts(vesicle,Xtar,fmm,NearV2T)
% InOut = sortPts(vesicle,Xtar,fmm,nearV2T) determines if the set of
% points in Xtar are inside or outside of a vesicle configuration

density = [ones(vesicle.N,vesicle.nv);...
           zeros(vesicle.N,vesicle.nv)];
% density function that is used to check if a point is inside or
% outside

tracers = capsules(Xtar,[],[],0,0,0);
if nargin == 3 
  [~,NearV2T] = vesicle.getZone(tracers,2);
end

op = poten(vesicle.N);
if ~fmm
  kernel = @op.exactLaplaceDL;
else
  kernel = @op.exactLaplaceDLfmm;
end
% kernel for Laplace's double layer potential

DLP = @(X) zeros(2*size(X,1),size(X,2));
% can cheat here because we know that the double-layer potential
% applied to the constant density function will always return zero

kernelDirect = kernel;
InOut = op.nearSingInt(vesicle,density,DLP,...
    NearV2T,kernel,kernelDirect,tracers,false,false);

InOut = InOut(1:end/2);
% only care about the first half as we are working with laplace as
% opposed to the vector-valued stokes layer potentials

thresh = 1e-4;
InOut(InOut > thresh) = 1;
InOut(InOut <= thresh) = 0;
% for points inside a vesicle, but close to its boundary, we will be
% interpolating the function [0;1;1;1;1...] close to the initial
% point.  Therefore, the value returned by the DLP will be some value
% between 0 and 1, but it can be quite close to 0.  So, we threshold

end % sortPts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function press = pressure(vesicle,f,RS,pressTar,fmm,LP)
% press = pressure(vesicle,f,RS,pressTar,fmm,LP) computes the pressure
% due to vesicle at the locations pressTar with near-singular
% integration with traction jump f and Rotlets stored in RS.  LP is
% either SLP or DLP and tells this routine if it needs to evaluate the
% pressure due to the single-layer potential or the double-layer
% potential.  Note that vesicle may correspond to solid walls rather
% than vesicles

N = vesicle.N; % points per vesicle
nv = vesicle.nv; % number of vesicles
op = poten(N);
Ntar = pressTar.N;
% number of target points where we want to compute the pressure
% and stress

[~,NearV2T] = vesicle.getZone(pressTar,2);
% build near-singular integration structure for vesicle
% to pressure target points
InOutFlag = find(vesicle.sortPts(pressTar.X,fmm,NearV2T) == 1);
%figure(1); clf; hold on
%plot(vesicle.X(1:end/2,:),vesicle.X(end/2+1:end,:),'r');
%oc = curve;
%[xx,yy] = oc.getXY(pressTar.X);
%plot(xx,yy,'b.');
%plot(xx(InOutFlag),yy(InOutFlag),'k.');
% DEBUG: MAKE SURE POINTS INSIDE AND OUTSIDE ARE LABELLED CORRECTLY

%fprintf('PRESSURE\n')
if strcmp(LP,'SLP')
  PdiagIn = op.pressSLmatrix(vesicle);
  PdiagOut = PdiagIn;
  % use odd-even integration to evaluate pressure due to
  % each vesicle independent of the others on its boundary
  for k = 1:nv
    PdiagOut(1:N,1:N,k) = PdiagOut(1:N,1:N,k) + ...
        1/2*diag(vesicle.xt(N+1:2*N,k));
    PdiagOut(1:N,N+1:2*N,k) = PdiagOut(1:N,N+1:2*N,k) + ...
        1/2*diag(-vesicle.xt(1:N,k));
    PdiagIn(1:N,1:N,k) = PdiagIn(1:N,1:N,k) - ...
        1/2*diag(vesicle.xt(N+1:2*N,k));
    PdiagIn(1:N,N+1:2*N,k) = PdiagIn(1:N,N+1:2*N,k) - ...
        1/2*diag(-vesicle.xt(1:N,k));
  end
  % Jump term that comes from pressure of single-layer potential is the
  % dot product of the traction jump with the normal vector multiplied
  % by 1/2
  Pdiag = @(X) op.exactPressureSLdiag(vesicle,[PdiagOut;PdiagIn],X);
  % nearSingInt assumes that input and output are vector-valued Take
  % adavntage of the fact that nearSingInt only works with
  % vector-valued densities by passing in two different jumps---one for
  % the limit from the inside and one from the outside.

  if ~fmm
    kernel = @op.exactPressSL;
  else
    kernel = @op.exactPressSLfmm;
  end
  kernelDirect = @op.exactPressSL;
elseif strcmp(LP,'DLP')
  oc = curve;
  PdiagIn = op.pressDLmatrix(vesicle);
  PdiagOut = PdiagIn;
  Deriv = fft1.D1(N);
  % spectral differentiation matrix
  for k = 1:nv
    [tanx,tany] = oc.getXY(vesicle.xt(:,k));
    % tangent vector
    sa = vesicle.sa(:,k);
    % Jacobian

    jump = -[diag(tanx./sa) diag(tany./sa)] * ...
        [Deriv zeros(N); zeros(N) Deriv];
    PdiagIn(:,:,k) = PdiagIn(:,:,k) - jump;
    PdiagOut(:,:,k) = PdiagOut(:,:,k) + jump;
    % add in the jump term.  Jump term is negative because target
    % points are interior to the solid walls
  end
  % Jump term that comes from pressure of double-layer potential
  % is the dot product of the tangential derivative of the traction 
  % jump with the tangent.  Assuming that all the pressure target
  % points are inside the physical domain so there is no use to
  % consider the limiting value from the other side
  Pdiag = @(X) op.exactPressureDLdiag(vesicle,[PdiagOut;PdiagIn],X);
  % nearSingInt assumes that input and output are vector-valued Take
  % adavntage of the fact that nearSingInt only works with
  % vector-valued densities by passing in two different jumps---one for
  % the limit from the inside and one from the outside.

  kernel = @op.exactPressDL;
  kernelDirect = @op.exactPressDL;
end

pressT = op.nearSingInt(vesicle,f,Pdiag,...
    NearV2T,kernel,kernelDirect,pressTar,false,false);
% compute the pressure of the single- or double-layer 
% potential using near-singular integration.  First row 
% corresponds to using the limiting value for the pressure 
% exterior of the vesicle and the second row corresponds to 
% using the limiting value for the pressure interior of the 
% vesicle
press = pressT(1:Ntar);
press1 = press;
press(InOutFlag) = pressT(Ntar+InOutFlag);
% At the interior points, take the second row

if ~isempty(RS)
  for k = 2:vesicle.nv
    cx = vesicle.center(1,k);
    cy = vesicle.center(2,k);
    % center of interior solid wall k
    xi1 = RS(3*(k-2) + 1,k);
    xi2 = RS(3*(k-2) + 2,k);
    % first and second components of the stokeslet on 
    % interior solid wall k 
    [x,y] = oc.getXY(pressTar.X);
    rx = x - cx;
    ry = y - cy;
    press = press +  2./(rx.^2 + ry.^2).*(rx*xi1 + ry*xi2);
    % add in pressure of stokeslet term
    % rotlet has a vanishing pressure gradient
  end
end

end % pressure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stress1,stress2] = stressTensor(vesicle,f,...
    RS,stressTar,fmm,LP)
% [stress1 stress2] = stressTensor(vesicle,f,RS,stressTar,fmm,LP)
% computes the stress tensor due to vesicle at the locations stressTar
% with or without near-singular integration with traction jump f.
% Returns in two components.  stress1 is the stress tensor applied to
% [1;0] and stress2 is the stress tensor applied to [0;1].  LP is
% either SLP or DLP and tells this routine if it needs to evaluate the
% stress due to the single-layer potential or the double-layer
% potential.  Note that vesicle may correspond to solid walls rather
% than vesicles.  If doing solid walls, there is a stress due to the
% rotlet and stokeslet terms that needs to be added in

N = vesicle.N; % points per vesicle
nv = vesicle.nv; % number of vesicles
op = poten(N);
Ntar = stressTar.N*stressTar.nv;
% number of points where we want to evaluate the stress

tangent = vesicle.xt;
oc = curve;
[tanx,tany] = oc.getXY(tangent);
nx = tany;
ny = -tanx;
% decompose tangent and normal vectors into their x and y
% coordinates

[~,NearV2T] = vesicle.getZone(stressTar,2);
% build near-singular integration structure for vesicle
% to pressure target points
InOutFlag = find(vesicle.sortPts(stressTar.X,fmm,NearV2T) == 1);

if strcmp(LP,'SLP')
  [S1diagIn,S2diagIn] = op.stressSLmatrix(vesicle);
  % use odd-even integration to evaluate stress due to
  % each vesicle independent of the others on its boundary
  % S1diagIn is the matrix that takes the traction jump and returns
  % the stress tensor applied to [1;0] evaluated on the boundary of
  % the vesicle
  % S2diagIn is the matrix that takes the traction jump and returns
  % the stress tensor applied to [0;1] evaluated on the boundary of
  % the vesicle
  S1diagOut = S1diagIn;
  S2diagOut = S2diagIn;
  % Use a second copy for the jump from the outside
  for k = 1:nv
    jump = +0.5*nx(:,k) + tanx(:,k).^2.*tany(:,k);
    S1diagIn(1:N,1:N,k) = S1diagIn(1:N,1:N,k) + diag(jump);
    S1diagOut(1:N,1:N,k) = S1diagOut(1:N,1:N,k) - diag(jump);
    % add in (1,1) jump component to stress applied to [1;0] 

    jump = +0.5*ny(:,k) + tanx(:,k).*tany(:,k).^2;
    S1diagIn(N+1:2*N,1:N,k) = S1diagIn(N+1:2*N,1:N,k) + diag(jump);
    S1diagOut(N+1:2*N,1:N,k) = S1diagOut(N+1:2*N,1:N,k) - diag(jump);
    % add in (2,1) jump component to stress applied to [1;0] 

    jump = +0.5*(tany(:,k).^2 - tanx(:,k).^2).*tanx(:,k);
    S1diagIn(1:N,N+1:2*N,k) = S1diagIn(1:N,N+1:2*N,k) + diag(jump);
    S1diagOut(1:N,N+1:2*N,k) = S1diagOut(1:N,N+1:2*N,k) - diag(jump);
    % add in (1,2) jump component to stress applied to [1;0] 

    jump = +0.5*(tany(:,k).^2 - tanx(:,k).^2).*tany(:,k);
    S1diagIn(N+1:2*N,N+1:2*N,k) = ...
        S1diagIn(N+1:2*N,N+1:2*N,k) + diag(jump);
    S1diagOut(N+1:2*N,N+1:2*N,k) = ...
        S1diagOut(N+1:2*N,N+1:2*N,k) - diag(jump);
    % add in (2,2) jump component to stress applied to [1;0] 

    jump = +0.5*(tany(:,k).^2 - tanx(:,k).^2).*tanx(:,k);
    S2diagIn(1:N,1:N,k) = S2diagIn(1:N,1:N,k) + diag(jump);
    S2diagOut(1:N,1:N,k) = S2diagOut(1:N,1:N,k) - diag(jump);
    % add in (1,1) jump component to stress applied to [0;1] 

    jump = +0.5*(tany(:,k).^2 - tanx(:,k).^2).*tany(:,k);
    S2diagIn(N+1:2*N,1:N,k) = S2diagIn(N+1:2*N,1:N,k) + diag(jump);
    S2diagOut(N+1:2*N,1:N,k) = S2diagOut(N+1:2*N,1:N,k) - diag(jump);
    % add in (2,1) jump component to stress applied to [0;1] 

    jump = +0.5*nx(:,k) - tanx(:,k).^2.*tany(:,k);
    S2diagIn(1:N,N+1:2*N,k) = S2diagIn(1:N,N+1:2*N,k) + diag(jump);
    S2diagOut(1:N,N+1:2*N,k) = S2diagOut(1:N,N+1:2*N,k) - diag(jump);
    % add in (1,2) jump component to stress applied to [0;1] 

    jump = +0.5*ny(:,k) - tanx(:,k).*tany(:,k).^2;
    S2diagIn(N+1:2*N,N+1:2*N,k) = ...
        S2diagIn(N+1:2*N,N+1:2*N,k) + diag(jump);
    S2diagOut(N+1:2*N,N+1:2*N,k) = ...
        S2diagOut(N+1:2*N,N+1:2*N,k) - diag(jump);
    % add in (2,2) jump component to stress applied to [0;1] 
  end
  % Jump term that comes from stress of single-layer potential
  % is the 
  S1diagInFn = @(X) op.exactStressSLdiag(vesicle,S1diagIn,X);
  S2diagInFn = @(X) op.exactStressSLdiag(vesicle,S2diagIn,X);
  S1diagOutFn = @(X) op.exactStressSLdiag(vesicle,S1diagOut,X);
  S2diagOutFn = @(X) op.exactStressSLdiag(vesicle,S2diagOut,X);
  % nearSingInt assumes that input and output are vector-valued
  kernel1 = @op.exactStressSL1;
  kernel2 = @op.exactStressSL2;

elseif strcmp(LP,'DLP')
  [S1diagIn,S2diagIn] = op.stressDLmatrix(vesicle);
  % use odd-even integration to evaluate stress due to
  % each vesicle independent of the others on its boundary
  % S1diagIn is the matrix that takes the traction jump and returns
  % the stress tensor applied to [1;0] evaluated on the boundary of
  % the vesicle
  % S2diagIn is the matrix that takes the traction jump and returns
  % the stress tensor applied to [0;1] evaluated on the boundary of
  % the vesicle

  S1diagOut = S1diagIn;
  S2diagOut = S2diagIn;
  % Use a second copy for the jump from the outside
  Deriv = fft1.D1(N);
  % Fourier differentiation matrix

  for k = 1:nv
    [tanx,tany] = oc.getXY(vesicle.xt(:,k));
    % tangent vector
    sa = vesicle.sa(:,k);
    % Jacobian

    jump = diag(1./sa.*(1+tanx.^2-tany.^2).*tanx)*Deriv;
    S1diagIn(1:N,1:N,k) = S1diagIn(1:N,1:N,k) - jump;
    S1diagOut(1:N,1:N,k) = S1diagOut(1:N,1:N,k) + jump;
    % add in (1,1) jump component to stress applied to [1;0] 

    jump = diag(1./sa.*(1+tanx.^2-tany.^2).*tany)*Deriv;
    S1diagIn(1:N,N+1:2*N,k) = S1diagIn(1:N,N+1:2*N,k) - jump;
    S1diagOut(1:N,N+1:2*N,k) = S1diagOut(1:N,N+1:2*N,k) + jump;
    % add in (1,2) jump component to stress applied to [1;0] 

    jump = diag(1./sa.*(2*tanx.*tany).*tanx)*Deriv;
    S1diagIn(N+1:2*N,1:N,k) = S1diagIn(N+1:2*N,1:N,k) - jump;
    S1diagOut(N+1:2*N,1:N,k) = S1diagOut(N+1:2*N,1:N,k) + jump;
    % add in (2,1) jump component to stress applied to [1;0] 

    jump = diag(1./sa.*(2*tanx.*tany).*tany)*Deriv;
    S1diagIn(N+1:2*N,N+1:2*N,k) = ...
        S1diagIn(N+1:2*N,N+1:2*N,k) - jump;
    S1diagOut(N+1:2*N,N+1:2*N,k) = ...
        S1diagOut(N+1:2*N,N+1:2*N,k) + jump;
    % add in (2,2) jump component to stress applied to [1;0] 


    jump = diag(1./sa.*(2*tanx.*tany).*tanx)*Deriv;
    S2diagIn(1:N,1:N,k) = S2diagIn(1:N,1:N,k) - jump;
    S2diagOut(1:N,1:N,k) = S2diagOut(1:N,1:N,k) + jump;
    % add in (1,1) jump component to stress applied to [0;1] 

    jump = diag(1./sa.*(2*tanx.*tany).*tany)*Deriv;
    S2diagIn(1:N,N+1:2*N,k) = S2diagIn(1:N,N+1:2*N,k) - jump;
    S2diagOut(1:N,N+1:2*N,k) = S2diagOut(1:N,N+1:2*N,k) + jump;
    % add in (1,2) jump component to stress applied to [0;1] 

    jump = diag(1./sa.*(1-tanx.^2+tany.^2).*tanx)*Deriv;
    S2diagIn(N+1:2*N,1:N,k) = S2diagIn(N+1:2*N,1:N,k) - jump;
    S2diagOut(N+1:2*N,1:N,k) = S2diagOut(N+1:2*N,1:N,k) + jump;
    % add in (2,1) jump component to stress applied to [0;1] 

    jump = diag(1./sa.*(1-tanx.^2+tany.^2).*tany)*Deriv;
    S2diagIn(N+1:2*N,N+1:2*N,k) = ...
        S2diagIn(N+1:2*N,N+1:2*N,k) - jump;
    S2diagOut(N+1:2*N,N+1:2*N,k) = ...
        S2diagOut(N+1:2*N,N+1:2*N,k) + jump;
    % add in (2,2) jump component to stress applied to [0;1] 
  end
  % Jump term that comes from stress of double-layer potential
  S1diagInFn = @(X) op.exactStressDLdiag(vesicle,S1diagIn,X);
  S2diagInFn = @(X) op.exactStressDLdiag(vesicle,S2diagIn,X);
  S1diagOutFn = @(X) op.exactStressDLdiag(vesicle,S1diagOut,X);
  S2diagOutFn = @(X) op.exactStressDLdiag(vesicle,S2diagOut,X);
  % nearSingInt assumes that input and output are vector-valued

  kernel1 = @op.exactStressDL1;
  kernel2 = @op.exactStressDL2;
end
% Have built all single- or double-layer potentials on boundary so that
% we can compute diagonal value to use near-singular integration

if numel(InOutFlag) > 0
%    fprintf('Stress1 Inside\n')
 stress1 = op.nearSingInt(vesicle,f,S1diagInFn,...
      NearV2T,kernel1,kernel1,stressTar,false,false);
else
  stress1 = zeros(2*Ntar,1);
end
if Ntar - numel(InOutFlag) > 0
%    fprintf('Stress1 Outside\n')
  stressT = op.nearSingInt(vesicle,f,S1diagOutFn,...
      NearV2T,kernel1,kernel1,stressTar,false,false);
  stress1(InOutFlag) = stressT(InOutFlag);
  stress1(Ntar+InOutFlag) = stressT(Ntar+InOutFlag);
  % correct stress at exterior points
end
% Use near-singular integration to compute first component of
% the stress due to the single-layer potential

if numel(InOutFlag) > 0
%    fprintf('Stress2 Inside\n')
  stress2 = op.nearSingInt(vesicle,f,S2diagInFn,...
      NearV2T,kernel2,kernel2,stressTar,false,false);
else
  stress2 = zeros(2*Ntar,1);
end
if Ntar - numel(InOutFlag) > 0
%    fprintf('Stress2 Outside\n')
  stressT = op.nearSingInt(vesicle,f,S2diagOutFn,...
      NearV2T,kernel2,kernel2,stressTar,false,false);
  stress2(InOutFlag) = stressT(InOutFlag);
  stress2(Ntar+InOutFlag) = stressT(Ntar+InOutFlag);
  % correct stress at exterior points
end
% Use near-singular integration to compute second component of
% the stress due to the single-layer potential

if ~isempty(RS)
  for k = 2:vesicle.nv
    cx = vesicle.center(1,k);
    cy = vesicle.center(2,k);
    % center of interior solid wall k
    xi = RS(3*(k-2) + 3,k);
    % rotlet on interior solid wall k
    xi1 = RS(3*(k-2) + 1,k);
    xi2 = RS(3*(k-2) + 2,k);
    % first and second components of the stokeslet on 
    % interior solid wall k 
    [x,y] = oc.getXY(stressTar.X);
    rx = x - cx;
    ry = y - cy;
    rho4 = (rx.^2 + ry.^2).^2;
    stress1(1:Ntar) = stress1(1:Ntar) - 4*xi*rx.*ry./rho4;
    % add in (1,1) component of stress due to Rotlet
    stress1(1+Ntar:2*Ntar) = stress1(1+Ntar:2*Ntar) + ...
        2*xi*(rx.^2 - ry.^2)./rho4;
    % add in (1,2) component of stress due to Rotlet
    stress2(1:Ntar) = stress2(1:Ntar) + ...
        2*xi*(rx.^2 - ry.^2)./rho4;
    % add in (2,1) component of stress due to Rotlet
    stress2(1+Ntar:2*Ntar) = stress2(1+Ntar:2*Ntar) + ...
        4*xi*rx.*ry./rho4;
    % add in (2,2) component of stress due to Rotlet


    stress1(1:Ntar) = stress1(1:Ntar) + ...
        2*(rx*xi1 + ry*xi2).*(-rx.^2 + ry.^2)./rho4;
    stress1(1:Ntar) = stress1(1:Ntar) - ...
        2*(rx*xi1 + ry*xi2)./(rx.^2+ry.^2);
    % add in (1,1) component of stress due to Stokeslet
    % need to add in -1*pressure term as well
    stress1(1+Ntar:2*Ntar) = stress1(1+Ntar:2*Ntar) - ...
        4*rx.*ry.*(rx*xi1 + ry*xi2)./rho4;
    % add in (1,2) component of stress due to Stokeslet
    stress2(1:Ntar) = stress2(1:Ntar) - ...
        4*rx.*ry.*(rx*xi1 + ry*xi2)./rho4;
    % add in (1,2) component of stress due to Stokeslet
    stress2(1+Ntar:2*Ntar) = stress1(1+Ntar:2*Ntar) + ...
        2*(rx*xi1 + ry*xi2).*(rx.^2 - ry.^2)./rho4;
    stress2(1+Ntar:2*Ntar) = stress2(1+Ntar:2*Ntar) - ...
        2*(rx*xi1 + ry*xi2)./(rx.^2+ry.^2);
    % add in (2,2) component of stress due to Stokeslet
    % need to add in -1*pressure term as well
  end
end
% Add in contributions due to Rotlets and Stokeslets

end % stressTensor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [shearStress,normalStress] = ...
    computeShearStress(vesicle,options,prams)

N = vesicle.N; % points per vesicle
nv = vesicle.nv; % number of vesicles

X = vesicle.X;
sigma = vesicle.sig;
sa = vesicle.sa;

f = vesicle.tracJump(X,sigma);
% traction jump of the vesicle shape and tension

if options.adhesion
  f = f + vesicle.adhesionTerm(prams.adStrength,prams.adRange);
end
% add on the adhesion contribution, if adhesion is included in the model

oc = curve;
[fx,fy] = oc.getXY(f);
% extract x and y coordinates of the traction jump
T11 = zeros(N,nv);
T12 = zeros(N,nv);
T22 = zeros(N,nv);
shearStress = zeros(N,nv);
normalStress = zeros(N,nv);
% initialize variables for stress tensor, shear stress, and normal
% stress

for ktar = 1:nv % loop over all vesicles
  ind = [(1:ktar-1) (ktar+1:nv)]; % all other vesicles
  if numel(ind) > 0
    for jtar = 1:N
      for ksou = ind
        rx = X(jtar,ktar) - X(1:end/2,ksou);
        ry = X(jtar+N,ktar) - X(end/2+1:end,ksou);
        rho2 = rx.^2 + ry.^2;
        rdotf = rx.*fx(:,ksou) + ry.*fy(:,ksou);
        coeff = -1/pi*rdotf./rho2.^2.*vesicle.sa(:,ksou)*2*pi/N;
        T11(jtar,ktar) = T11(jtar,ktar) + sum(coeff.*rx.*rx);
        T12(jtar,ktar) = T12(jtar,ktar) + sum(coeff.*rx.*ry);
        T22(jtar,ktar) = T22(jtar,ktar) + sum(coeff.*ry.*ry);
      end
    end
    % loop over all points on all vesicles except self
  end

  for jtar = 1:2:N  % odd-indexed points
    ind = (2:2:N); % use even indexed points
    rx = X(jtar,ktar) - X(ind,ktar);
    ry = X(jtar+N,ktar) - X(ind+N,ktar);
    rho2 = rx.^2 + ry.^2;
    rdotf = rx.*fx(ind,ktar) + ry.*fy(ind,ktar);
    coeff = -1/pi*rdotf./rho2.^2.*vesicle.sa(ind,ktar)*2*pi/N*2;
    T11(jtar,ktar) = T11(jtar,ktar) + sum(coeff.*rx.*rx);
    T12(jtar,ktar) = T12(jtar,ktar) + sum(coeff.*rx.*ry);
    T22(jtar,ktar) = T22(jtar,ktar) + sum(coeff.*ry.*ry);
  end

  for jtar = 2:2:N  % even-indexed points
    ind = (1:2:N); % use odd indexed points
    rx = X(jtar,ktar) - X(ind,ktar);
    ry = X(jtar+N,ktar) - X(ind+N,ktar);
    rho2 = rx.^2 + ry.^2;
    rdotf = rx.*fx(ind,ktar) + ry.*fy(ind,ktar);
    coeff = -1/pi*rdotf./rho2.^2.*vesicle.sa(ind,ktar)*2*pi/N*2;
    T11(jtar,ktar) = T11(jtar,ktar) + sum(coeff.*rx.*rx);
    T12(jtar,ktar) = T12(jtar,ktar) + sum(coeff.*rx.*ry);
    T22(jtar,ktar) = T22(jtar,ktar) + sum(coeff.*ry.*ry);
  end
end

[tx,ty] = oc.getXY(vesicle.xt);
nx = ty;
ny = -tx;

jump11 = +1/2*nx.*fx + 1/2*tx.*(2*tx.*ty.*fx + (ty.^2 - tx.^2).*fy);
jump12 = +1/2*nx.*fy + 1/2*tx.*((ty.^2-tx.^2).*fx - 2*tx.*ty.*fy);
jump22 = +1/2*ny.*fy + 1/2*ty.*((ty.^2-tx.^2).*fx - 2*tx.*ty.*fy);

T11 = T11 - 1*jump11;
T12 = T12 - 1*jump12;
T22 = T22 - 1*jump22;

tt = tstep(options,prams);
[~,bg11,bg12,bg22] = tt.farField(X);

T11 = 1*T11 + 1*bg11;
T12 = 1*T12 + 1*bg12;
T22 = 1*T22 + 1*bg22;

%xtar = 0.3;
%ytar = 1.2;
%rx = xtar - X(1:end/2);
%ry = ytar - X(end/2+1:end);
%rho2 = rx.^2 + ry.^2;
%rdotf = rx.*fx + ry.*fy;
%coeff = -1/pi*rdotf./rho2.^2.*vesicle.sa*2*pi/N;
%T11_tar = sum(coeff.*rx.*rx);
%T11_exact = -(2*xtar^5 + 3*xtar*ytar^2-xtar^3+2*xtar^3*ytar^2)/...
%    (xtar^2 + ytar^2)^3;
%T12_tar = sum(coeff.*rx.*ry);
%T12_exact = -(2*xtar^4*ytar+ytar^3-3*xtar^2*ytar+2*xtar^2*ytar^3)/...
%    (xtar^2 + ytar^2)^3;
%T22_tar = sum(coeff.*ry.*ry);
%press = 1/2/pi*sum(rdotf./rho2.*vesicle.sa)*2*pi/N;
%exact_press = xtar/(xtar^2 + ytar^2);
%T11_tar 
%T11_exact
%T12_tar 
%T12_exact
%press 
%exact_press
%pause

for ksou = 1:nv
  shearStress(:,ksou) = ...
    (T11(:,ksou).*nx(:,ksou) + T12(:,ksou).*ny(:,ksou)).*...
        tx(:,ksou) + ...
    (T12(:,ksou).*nx(:,ksou) + T22(:,ksou).*ny(:,ksou)).*...
        ty(:,ksou);
  normalStress(:,ksou) = ...
    (T11(:,ksou).*nx(:,ksou) + T12(:,ksou).*ny(:,ksou)).*...
        nx(:,ksou) + ...
    (T12(:,ksou).*nx(:,ksou) + T22(:,ksou).*ny(:,ksou)).*...
        ny(:,ksou);
end

%x = X(1:end/2,:);
%y = X(end/2+1:end,:);
%tracx = zeros(N,nv);
%tracy = zeros(N,nv);
%tracx = T11.*nx + T12.*ny;
%tracy = T12.*nx + T22.*ny;
%figure(1); clf; hold on
%plot(x,y,'r')
%quiver(x,y,tracx,tracy)
%axis equal
%clf
%z = (tracx.*y - tracy.*x).*vesicle.sa;
%plot(z(:,1))
%sum(z)
%pause

end % computeShearStress

end % methods

end %capsules



