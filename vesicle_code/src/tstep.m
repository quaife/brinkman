classdef tstep < handle
% This class defines the functions required to advance the geometry
% forward in time.  Handles both implicit and explicit vesicle-vesicle
% interactions, different inextensibility conditions, viscosity
% contrast, solid walls vs. unbounded flows.  This class also
% implements the adaptive time stepping strategy where the errors in
% length, area, and the residual are monitored to adaptively choose a
% time step size.


properties
Xcoeff        % Coefficients that are used in the discretization of the
              % time derivative when doing IMEX.  For example, if doing
              % IMEX Euler, Xcoeff = [1]
rhsCoeff      % Coefficients that are used in the discretization of the
              % explicit terms of IMEX.  In our case, this is
              % configuration of the arclength and layer potentials are
              % discretized.  For example, if doing IMEX Euler,
              % rhscoeff = [1]
beta          % Term multiplying implicit term in discretization of
              % derivative
order         % Time stepping order
expectedOrder % order that is used to pick new time step size
              % (parameter used in adaptive time stepping)
dt            % Time step size
currentTime   % current time needed for adaptive time stepping
finalTime     % time horizon
solver        % method1, method2, method3, or method4.  Changes how
              % inextensiblity condition is treated
Galpert       % Single-layer stokes potential matrix using Alpert
D             % Double-layer stokes potential matrix
lapDLP        % Double-layer laplace potential matrix
wallDLP       % Double-layer potential due to solid walls
wallN0        % Modificiation of double-layer potential to 
              % remove the rank deficiency on outer boundary
farField      % Far field boundary condition
confined      % whether geometry is bounded or not
bdiagVes      % precomputed inverse of block-diagonal precondtioner
              % only for vesicle-vesicle interactions
bdiagTen
bdiagWall     % precomputed inverse of block-diagonal precondtioner
              % only for wall-wall interactions
vesves        % Discretization of vesicle-vesicle interactions.
              % Can be implicit or explicit
fmm           % with or without the FMM
fmmDLP        % use FMM for the double-layer potential
near          % with or without near-singular integration
profile       % display profile times throughout simulation
gmresTol      % GMRES tolerance
gmresMaxIter  % maximum number of gmres iterations
collDete      % decides if we do collision detection
timeAdap      % Using adaptive time stepping
orderGL       % Order of Gauss-Lobatto points
GLpts         % Gauss-Lobatto points 
SDCcorrect    % timeStep changes depending on if we are doing
              % an sdc correction or forming a provisional
              % solution
rtolArea      % maximum allowable error in area
rtolLength    % maximum allowable error in length
nsdc          % number of sdc iterations to take
NearV2V       % near-singular integration for vesicle to vesicle
NearW2V       % near-singular integration for wall to vesicle
NearV2W       % near-singular integration for vesicle to wall 
NearW2W       % near-singular integration for wall to wall 
op            % class for evaluating potential so that we don't have
              % to keep building quadrature matricies
opWall        % class for walls              
betaUp        % maximum amount time step is scaled up
betaDown      % maximum amount time step is scaled down
alpha         % safety factor scaling for new time step size
adhesion      % use adhesion in the model
adStrength    % strength of adhesion
adRange       % range of adhesion
fmmPrecision  % precision of fmm
end % properties

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = tstep(options,prams)
% o=tstep(options,prams: constructor.  Initialize class.  Take all
% elements of options and prams needed by the time stepper

o.order = options.order; % Time stepping order
o.expectedOrder = options.expectedOrder; % Expected time stepping order
o.dt = prams.T/prams.m; % Time step size
[o.Xcoeff,o.rhsCoeff,o.beta] = o.getCoeff(o.order);
% Get coefficients for doing time integration
o.currentTime = 0;
% Method always starts at time 0
o.finalTime = prams.T;
% Need the time horizon for adaptive time stepping
o.solver = options.inextens;
% Discretization of inextensibility
o.vesves = options.vesves;
% Vesicle-vesicle interactions
o.fmm = options.fmm;
o.fmmDLP = options.fmmDLP;
% fast-multipole method
o.near = options.near;
% near-singular integration
o.profile = options.profile;
% disp profile times
o.gmresTol = prams.gmresTol;
% GMRES tolerance
o.gmresMaxIter = prams.gmresMaxIter;
% maximum number of GMRES iterations
o.farField = @(X) o.bgFlow(X,options.farField,options.farFieldSpeed);
% Far field boundary condition built as a function handle
o.confined = options.confined;
% Confined or unbounded geometry
o.collDete = options.collision;
% Collision detection
o.timeAdap = options.timeAdap;
% Time adaptivity flag
o.orderGL = options.orderGL;
% Gauss-Lobatto order
o.GLpts = o.gaussLobatto(o.orderGL);
% load Gauss-Lobatto points
% if options.correctShape
  o.rtolArea = prams.rtolArea;
  o.rtolLength = prams.rtolLength;
% else
%   o.rtolArea = prams.rtolArea/prams.T;
%   o.rtolLength = prams.rtolLength/prams.T;
% end
% if the shape is being corrected after each time step, then the
% rtolArea and rtolLength are exactly the tolerance for the error
% accumulated at each time step.  However, if the shape is not being
% corrected, then the tolerance should be achieved at the time horizon,
% so we need to nondimensionalize by dividing by the time horizon
o.nsdc = options.nsdc;
% number of sdc iterations to take
if o.order > 1
  o.SDCcorrect = false;
end

o.betaUp = prams.betaUp;
o.betaDown = prams.betaDown;
o.alpha = prams.alpha;
% safety factors for adaptive time stepping

o.adhesion = options.adhesion;
o.adStrength = prams.adStrength;
% strength of adhesion
o.adRange = prams.adRange;
% range of adhesion

o.fmmPrecision = options.fmmPrecision;
o.op = poten(prams.N,o.fmmPrecision);

if options.confined
  o.opWall = poten(prams.Nbd,o.fmmPrecision);
else
  o.opWall = [];
end

% Build class with quadrature points and weights as well as lagrange
% interpolation matrix
end % tstep: constructor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xcoeff,rhsCoeff,beta] = getCoeff(o,order)
% [Xcoeff,rhsCoeff,beta] = getCoeff(order) generates the coefficients
% required to discretize the derivative.  First-order time  derivatives
% are approximated by beta*x^{N+1} + rhscoeff.*[x^{N} x^{N-1} ...]
% Explicit terms (operators) are discretized at Xcoeff.*[x^{N} x^{N-1}
% ...] All rules are from Ascher, Ruuth, Wetton 1995.

if (order == 4) % fourth-order
  beta = 25/12;
  Xcoeff = [-1; 4; -6; 4]; 
  rhsCoeff = [-1/4; 4/3; -3; 4];
elseif (order == 3) % third-order
  beta = 11/6;
  Xcoeff = [1; -3; 3];
  rhsCoeff = [1/3; -3/2; 3];
elseif (order == 2) % second-order
  beta = 1.5;
  Xcoeff = [-1; 2];
  rhsCoeff = [-0.5; 2];
else % first-order
  beta = 1;
  Xcoeff = 1;
  rhsCoeff = 1;
end

end % getCoeff

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [walls,wallsCoarse] = initialConfined(o,prams,Xwalls)
% walls = initialConfined(prams,Xwalls) builds an object of class
% capsules for the solid walls
potWall = o.opWall;

uwalls = o.farField(Xwalls);
% velocity on solid walls coming from no-slip boundary condition
walls = capsules(Xwalls,[],uwalls,zeros(prams.nvbd,1),...
    zeros(prams.nvbd,1));

wallsCoarse = walls;

o.wallDLP = potWall.stokesDLmatrix(wallsCoarse);
o.wallN0 = potWall.stokesN0matrix(wallsCoarse);
o.bdiagWall = o.wallsPrecond(wallsCoarse);
% block diagonal preconditioner for solid walls build the double-layer
% potential and the matrix N0 that removes the rank one null space from
% the double-layer potential for the solid walls

end % initialConfined


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xstore,sigStore,uStore,etaStore,RSstore,Xtra] = ...
  firstSteps(o,options,prams,Xinit,sigInit,uInit,...
  walls,wallsCoarse,om,Xtra,pressTar)
% [Xstore,sigStore,uStore,etaStore,RSstore] = ...
%   firstSteps(options,prams,Xinit,sigInit,uInit,...
%   walls,wallsCoarse,om,pressTar)
% refines the first time step [0,dt] and uses a first-order, then
% second-order, then third-order, etc time stepping to find the vesicle
% position, tension, velocity, and density function eta defined on the
% solid walls and the rotlets and stokeslets at t = dt.  Returns
% Xstore, sigStore, uStore etaStore, so that they are immediatly ready
% to use for higher-order time stepping.  This routine also computes
% the tension and density functions of the initial condition.  This is
% needed in SDC since we require these functions at all Gauss-Lobatto
% points in order to compute the residual

N = size(Xinit,1)/2; % points per vesicle
nv = size(Xinit,2); % number of vesicles
if o.confined 
  Xwalls = walls.X; % discretization points of solid walls
  Nbd = size(Xwalls,1)/2; % points per wall
  nvbd = size(Xwalls,2); % number of wall components
else
  Xwalls = [];
  Nbd = 0;
  nvbd = 0;
end

prams.T = prams.T/prams.m*(o.order-1);
% time horizon is enough time steps so that we can
% continue with higher-order time integrator
if o.order ~=1
  mR = min(ceil(prams.m/32),100)*o.order^2;
  mR = mR * (o.order - 1);
else
  mR = 1;
end
% For second-order time stepping, this keeps the local error from time
% step 0 to dt on the order of dt^3
prams.m = mR*(o.order-1);
% number of time steps to take in range [0,dt]

Xstore = zeros(2*N,nv,prams.m+1);
sigStore = zeros(N,nv,prams.m+1);
uStore = zeros(2*N,nv,prams.m+1);
etaStore = zeros(2*Nbd,nvbd,prams.m+1);
RSstore = zeros(3,nvbd,prams.m+1);

Xstore(:,:,1) = Xinit;
vesicle = capsules(Xinit,zeros(N,nv),[],...
    prams.kappa,prams.viscCont);

[sigStore(:,:,1),etaStore(:,:,1),RSstore(:,:,1),u,iter] = ...
    vesicle.computeSigAndEta(o,walls);

% need intial tension, density function, rotlets, and stokeslets so
% that we can do SDC updates, couling with advection-diffusion system,
% etc

om.initializeFiles(Xinit,sigStore(:,:,1),etaStore(:,:,1),...
    RSstore(:,:,1),Xwalls,Xtra,pressTar);
% delete previous versions of files and write some initial
% options/parameters to files and the console

message = ['GMRES took ',num2str(iter),...
    ' iterations to find intial tension and density function'];
om.writeMessage(message,'%s\n');
om.writeMessage(' ','%s\n');

for n = 1:prams.m
  time = n*prams.T/prams.m;

  options.order = min(n,o.order);
  tt = tstep(options,prams);
  % take the highester order possible
  tt.SDCcorrect = false;
  % no SDC corrections
  tt.wallDLP = o.wallDLP;
  tt.wallN0 = o.wallN0;
  tt.bdiagWall = o.bdiagWall;
  % build inverse of the wall-to-wall intearctions includes diagonal and
  % off-diagonal terms.  That is, it could be used to solve the stokes
  % equations in a multiply-connected domain with no vesicles
 
  updatePreco = true;
  [X,sigma,u,eta,RS,iter,iflag] = tt.timeStep(...
      Xstore(:,:,n-tt.order+1:n),...
      sigStore(:,:,n-tt.order+1:n),...
      uStore(:,:,n-tt.order+1:n),...
      etaStore(:,:,n-tt.order+1:n),...
      RSstore(:,:,n-tt.order+1:n),...
      [],[],[],[],[],[],[],...
      prams.kappa,prams.viscCont,walls,wallsCoarse,...
      updatePreco,vesicle);
  % take one time step

  if numel(Xtra) > 1
    vel = o.tracersVel(X,sigma,u,...
        prams.kappa,prams.viscCont,walls,eta,RS,Xtra);
    Xtra = Xtra + tt.dt*vel;
  end

  accept = true;
  dtScale = 0;
  res = 0;
  % Required for adaptive time stepping which have not been
  % implemented yet for high-order time stepping

  Xstore(:,:,n+1) = X;
  sigStore(:,:,n+1) = sigma;
  uStore(:,:,n+1) = u;
  etaStore(:,:,n+1) = eta;
  RSstore(:,:,n+1) = RS;
  % save and update position, tension, velocity, and
  % density function

  terminate = om.outputInfo(X,sigma,u,eta,RS,...
      Xwalls,Xtra,time,iter,dtScale,res,iflag);
  % save data, write to log file, write to console as
  % requested
end 
% end of using small time steps to get the simulation far enough
% to use the desired order time integrator

if o.order > 1
  if options.pressure
    op = poten(N,o.fmmPrecision);
    [press,stress1,stress2] = op.pressAndStress(...
        X,sigma,u,prams.kappa,prams.viscCont,...
        walls,pressTar,eta,RS,options.confined,...
        options.fmm);
    % compute the pressure and stress due to the vesicles
    % and the solid walls
    if options.saveData
      om.writePressure(press);
      om.writeStress(stress1(1:end/2),'11');
      om.writeStress(stress1(1+end/2:end),'12');
      om.writeStress(stress2(1:end/2),'21');
      om.writeStress(stress2(1+end/2:end),'22');
    end
  end % ~options.pressure
  % write pressure and stress to output file at time dt

  Xstore = Xstore(:,:,1:mR:end);
  sigStore = sigStore(:,:,1:mR:end);
  uStore = uStore(:,:,1:mR:end);
  etaStore = etaStore(:,:,1:mR:end);
  RSstore = RSstore(:,:,1:mR:end);
end
% Only want every mR^{th} solution as these correspond to
% time steps that are multiples of dt
% Only need to do this if the time stepping order is not 1

end % firstSteps


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,sigma,u,eta,RS,iter,accept,dtScale,normRes,iflag] = timeStepGL(o,Xstore,sigStore,uStore,...
    etaStore,RSstore,kappa,viscCont,walls,wallsCoarse,om,time,accept)
% [X,sigma,u,eta,RS,iter,dt,accept,dtScale,normRes,iflag] = ...
%    timeStepGL(Xstore,sigStore,uStore,...
%    etaStore,RSstore,kappa,viscCont,walls,wallsCoarse,a0,l0)
% takes the desired number of time steps at Gauss-Lobatto points to
% find solution at dt given solution at 0.  Calls o.timeStep which is
% what we use for the constant sized time step integration.  Returns a
% flag that says if the solution was accepted or not and the amount the
% time step was scaled.  errors is the norm of the residual and iflag
% is tells if any of the gmres runs failed

N = size(Xstore,1)/2; % number of points per vesicle
nv = size(Xstore,2); % number of vesicles
Nbd = size(etaStore,1)/2; % number of points per solid wall
nvbd = size(etaStore,2); % number of solid walls
a0 = om.area; % input area
l0 = om.length; % input length
oc = curve;

Xprov = zeros(2*N,nv,abs(o.orderGL));
sigmaProv = zeros(N,nv,abs(o.orderGL));
uProv = zeros(2*N,nv,abs(o.orderGL));
etaProv = zeros(2*Nbd,nvbd,abs(o.orderGL));
RSprov = zeros(3,nvbd,abs(o.orderGL));
% need to save the position, tension, velocity, density
% functions, rotlets/stokeslets at the Gauss-Lobatto
% time steps
for k = 1:nv
  Xprov(:,k,1) = Xstore(:,k);
  sigmaProv(:,k,1) = sigStore(:,k);
  uProv(:,k,1) = uStore(:,k);
end
% Store the initial conditions in the provisional
% solution which will store the solution at all
% Gauss-Lobatto points

for k = 1:nvbd
  etaProv(:,k,1) = etaStore(:,k);
  RSprov(:,k,1) = RSstore(:,k);
end

Galpert = zeros(2*N,2*N,nv,abs(o.orderGL));
% need the single-layer potential at all levels of the provisional
% solution
if any(viscCont ~= 1)
  D = zeros(2*N,2*N,nv,abs(o.orderGL));
else
  D = [];
end
% need the double-layer potential at all levels
% of the provisional solution

deltaX = zeros(2*N,nv,abs(o.orderGL));
deltaSigma = zeros(N,nv,abs(o.orderGL));
deltaEta = zeros(2*Nbd,nvbd,abs(o.orderGL));
deltaRS = zeros(3,nvbd,abs(o.orderGL));
% Errors that are solved for in each iteration.  The first column of
% deltaX will always be zero since we assume the inital condition is
% exact

X = Xprov(:,:,1);
sigma = sigmaProv(:,:,1);
u = uProv(:,:,1);

vesicle(1) = capsules(X,sigma,u,kappa,viscCont);
% vesicle configuration at first Gauss-Lobatto point

if o.timeAdap
  Xinit_tadap = X;
  [~,aInit,lInit] = oc.geomProp(X);
  eaInit = abs(aInit - a0)./a0;
  elInit = abs(lInit - l0)./l0;
end
% compute intial reduced area, area, and length if using adaptive
% time stepping

dt = o.dt;
% need to save the time step size
dtimes = diff(o.GLpts)*dt/2;
% time step sizes use for Gauss-Lobatto points

% START OF FORMING PROVISIONAL SOLUTION
iflag = 0; % initialize gmres flag as everything okay
iter = 0; % initialize number of gmres iterations
updatePreco = true; 
% Form preconditioner at initial configuartion
for k = 1:abs(o.orderGL)-1
  o.dt = dtimes(k);
  o.SDCcorrect = false;

  [X,sigma,u,eta,RS,subIter,iflagTemp] = o.timeStep(...
      Xprov(:,:,k-o.order+1:k),sigmaProv(:,:,k-o.order+1:k),...
      uProv(:,:,k-o.order+1:k),etaProv(:,:,k-o.order+1:k),...
      RSprov(:,:,k-o.order+1:k),[],[],[],[],[],[],[],...
      kappa,viscCont,walls,wallsCoarse,updatePreco,vesicle(k));
  % form provisional solution at next Gauss-Lobatto point
  iter = iter + subIter;
  % running counter for the total number of gmres iterations
  % per time step.  GMRES is used at each Gauss-lobatto time
  % substep which we need to add up for a fair comparison to
  % the traditional time stepping method

  updatePreco = false;
  % don't want to precompute the preconditioner at the next time
  % step to save CPU time.  Since the configuration doesn't 
  % change by much, the preconditioner shouldn't change
  % by much either and the old preconditioner should
  % work fine

  if iflagTemp ~= 0
    iflag = iflagTemp;
  end
  % if any of the gmres iterations fail, assign the failure to
  % flag to iflag for monitor to output

  if o.nsdc > 0
    Galpert(:,:,:,k) = o.Galpert;
    if any(viscCont ~= 1)
      D(:,:,:,k) = o.D;
    end
    % want to save single-layer potential for computing the residual.
    % This will not save the  one, but this is computed in
    % computeResidual
    NearV2V{k} = o.NearV2V;
    NearW2V{k} = o.NearW2V;
    NearV2W{k} = o.NearV2W;
    NearW2W{k} = o.NearW2W;

    vesicle(k+1) = capsules(X,sigma,u,kappa,viscCont);
    % need to save if doing SDC corrections 
  end

  Xprov(:,:,k+1) = X;
  sigmaProv(:,:,k+1) = sigma;
  uProv(:,:,k+1) = u;
  etaProv(:,:,k+1) = eta;
  RSprov(:,:,k+1) = RS;
  % save the provisional solution at the intermediate 
  % Gauss-Lobatto point
end % k


% need to save the near-singular integration structure and
% layer-potential matricies at the final Gauss-Lobatto point for future
% SDC sweeps and computing the residual, if there are any SDC sweeps
if o.nsdc > 0
  if o.near
    if o.confined
      [NearV2V{abs(o.orderGL)},NearV2W{abs(o.orderGL)}] = ...
          vesicle(abs(o.orderGL)).getZone(walls,3);
      %--------------------------------------------------------------------
      % Need vesicle to vesicle and vesicle to wall interactions
      if nvbd == 1 
        [NearW2W{abs(o.orderGL)},NearW2V{abs(o.orderGL)}] = ...
            walls.getZone(vesicle(abs(o.orderGL)),2);
      else
        [NearW2W{abs(o.orderGL)},NearW2V{abs(o.orderGL)}] = ...
            walls.getZone(vesicle(abs(o.orderGL)),3);
      end
      % Only need wall to vesicle interactions.  Wall to wall 
      % interactions should also use near-singular integration since
      % they may be close to one another
    else
      NearV2V{abs(o.orderGL)} = ...
          vesicle(abs(o.orderGL)).getZone([],1);
      % no solid walls, so only need vesicle-vesicle intearactions
      NearV2W{abs(o.orderGL)} = [];
      NearW2V{abs(o.orderGL)} = [];
      NearW2W{abs(o.orderGL)} = [];
    end
  else
    NearV2V{abs(o.orderGL)} = [];
    NearV2W{abs(o.orderGL)} = [];
    NearW2V{abs(o.orderGL)} = [];
    NearW2W{abs(o.orderGL)} = [];
  end
  % Need near-singular integration strucutures for final state
  % for sdc updates and computing the residual
  Galpert(:,:,:,abs(o.orderGL)) = ...
      o.op.stokesSLmatrix(vesicle(o.orderGL));
  if any(viscCont ~= 1)
    D(:,:,:,abs(o.orderGL)) = ...
        o.op.stokesDLmatrix(vesicle(o.orderGL));
  end
  % need single-layer potential at the final state of the 
  % provisional solution to form the residual
end
% END OF FORMING PROVISIONAL SOLUTION AND THE SINGLE-LAYER
% POTENTIAL AT THESE INTERMEDIATE STATES

o.dt = dt;
% change back to original time step
 
%color = ['r' 'g' 'b' 'k' 'c' 'm' 'y'];
%figure(2); clf; hold on
%for k = 1:abs(o.orderGL)
%  plot(sigmaProv(:,1,k),color(k))
%end
%figure(3); clf; hold on
%for k = 1:abs(o.orderGL)
%  plot(squeeze(etaProv(1:end/2,1,k)),color(k))
%end
%figure(5);clf;hold on;
%for k = 1:abs(o.orderGL)
%  plot(etaProv(end/2+1:end,1,k),color(k))
%end
%pause
%% DEBUG: TO MAKE SURE THAT THE TENSION IS CONTINUOUS FROM BETWEEN
%% THE DIFFERENT GAUSS-LOBATTO POINTS

Xinit = Xprov;
uInit = uProv;
sigmaInit = sigmaProv;
etaInit = etaProv;
RSinit = RSprov;
% save the very first provisional solution so that it can be taken
% as the solution if sdc corrections increase the error but this
% original provisional solution has the smallest error
if o.nsdc > 0
  [vesVel,divVesVel,wallVel,residual] = ...
      o.computeVelocities(vesicle,Galpert,D,walls,...
        etaProv,RSprov,NearV2V,NearV2W,NearW2V,NearW2W);
  % form the residual as well as the velocities due to the different
  % components which are necessary to form the right-hand sides when
  % doing sdc corrections

  normRes = max(abs(residual(:,:,end)));
  normRes = max(normRes);
  % save the size of the residual at the final Gauss-Lobatto
  % point.  We may have more of these later if we do a full
  % deferred correction method.  Use the maximum L2 error where
  % the maximum is taken over all the vesicles
else
  normRes = 0;
end
% compute the integrand that is in the residual and the residual 
% of the picard integral formulation of the dynamic equation
% If the residual is not desired, then errors is set to 0 and
% this will never refine the time step size

[~,a,l] = oc.geomProp(Xprov(:,:,end));
eaBefore = abs(a - a0)./abs(a0);
elBefore = abs(l - l0)./abs(l0);
message = ['sdcCount   ' num2str(0,'%2d') ...
  ': Residual = ' num2str(normRes(end),'%5.2e') ...
  ', eA = ' num2str(max(eaBefore),'%5.4e') ...
  ', eL = ' num2str(max(elBefore),'%5.4e')];
om.writeMessage(message,'%s\n');
% print the provisional solution's residual and errors

% Start doing SDC corrections
for sdcCount = 1:o.nsdc
  updatePreco = false; 
  % use preco that was formed at the provisional solution
  for n = 1:abs(o.orderGL) - 1
    o.dt = dtimes(n);

    o.SDCcorrect = true;
    o.Galpert = Galpert(:,:,:,n+1);
    if any(viscCont ~= 1)
      o.D = D(:,:,:,n+1);
    end
    o.NearV2V = NearV2V{n+1};
    o.NearV2W = NearV2W{n+1};
    o.NearW2V = NearW2V{n+1};
    o.NearW2W = NearW2W{n+1};
    [X,sigma,u,eta,RS,subIter,iflagTemp] = o.timeStep(...
        Xprov(:,:,n+1),sigmaProv(:,:,n+1),...
        uStore,etaProv(:,:,n+1),RSprov(:,:,n+1),...
        deltaX(:,:,n),deltaSigma(:,:,n),...
        deltaEta(:,:,n),deltaRS(:,:,n),...
        residual(:,:,n+1) - residual(:,:,n),...
        vesVel(:,:,n+1),wallVel(:,:,n+1),...
        kappa,viscCont,walls,wallsCoarse,updatePreco,vesicle(n+1),...
        vesicle(1).sa,vesicle(1).IK);
    % Form the sdc update
    iter = iter + subIter;

    if iflagTemp ~= 0
      iflag = iflagTemp;
    end
    o.SDCcorrect = false;
    % turn correct off since it is not used when forming the
    % provisional solution
    deltaX(:,:,n+1) = X;
    deltaSigma(:,:,n+1) = sigma;
    deltaEta(:,:,n+1) = eta;
    deltaRS(:,:,n+1) = RS;
    % approximations of the error
  end
  o.dt = dt;
  % go back to original time step
  
  alpha = 1;
%  eaOld = ea;
%  elOld = el;
%  ea = 2*eaOld; el = 2*elOld;
%  while ((ea > eaOld || el > elOld) && 2*alpha > 1e-2)
%    [~,a,l] = oc.geomProp(Xprov(:,:,end) + alpha*deltaX(:,:,end));
%    ea = abs(a - a0)./abs(a0);
%    el = abs(l - l0)./abs(l0);
%    alpha = alpha/2;
%  end
%  alpha = 2*alpha;
%  % do a line search so that the updated solution has a smaller area
%  % and length
%  if alpha < 1e-2
%    alpha = 0;
%  end
%  % if the line search was unsuccessful, don't update the provisional
%  % solution

  Xprov = Xprov + alpha*deltaX;
  sigmaProv = sigmaProv + alpha*deltaSigma;
  etaProv = etaProv + alpha*deltaEta;
  RSprov = RSprov + alpha*deltaRS;
  % update provision solution
  

  [~,a,l] = oc.geomProp(Xprov(:,:,end));
  ea = abs(a - a0)./abs(a0);
  el = abs(l - l0)./abs(l0);

  normRes = max(abs(residual(:,:,end)));
  normRes = max(normRes);

  message = ['sdcCount   ' num2str(sdcCount,'%2d') ...
    ': Residual = ' num2str(normRes(end),'%5.2e') ...
    ', eA = ' num2str(max(ea),'%5.4e') ...
    ', eL = ' num2str(max(el),'%5.4e')];
  om.writeMessage(message,'%s\n');
  
  if max(ea) > max(eaBefore) || max(el) > max(elBefore)
    message = ['SDC increases the error, do not accept corrected solution'];
    om.writeMessage(message,'%s\n');

    Xprov = Xprov - alpha*deltaX;
    sigmaProv = sigmaProv - alpha*deltaSigma;
    etaProv = etaProv - alpha*deltaEta;
    RSprov = RSprov - alpha*deltaRS;    
    % Do not accept SDC correction
  end
  
  if sdcCount < o.nsdc
    for n = 1:abs(o.orderGL)
      vesicle(n) = capsules(Xprov(:,:,n),sigmaProv(:,:,n),...
          [],kappa,viscCont);
      % find correct upsampling rate
      Galpert(:,:,:,n) = o.op.stokesSLmatrix(vesicle(n));
      if any(viscCont ~= 1)
        D(:,:,:,n) = o.op.stokesDLmatrix(vesicle(n));
      else
        D = [];
      end
    end
    % update the vesicle objects and the single-layer potentials
    
    [vesVel,divVesVel,wallVel,residual] = ...
      o.computeVelocities(vesicle,Galpert,D,walls,...
          etaProv,RSprov,NearV2V,NearV2W,NearW2V,NearW2W);
  end
  % if we are only recording the error in area and length, don't need
  % to compute the final residual
  
end
% End of doing SDC corrections
% update solution with an SDC iteration

if o.timeAdap
  [accept,dtScale] = o.newTimeStepSize(a,l,...
      aInit,lInit,accept,om);
  % if doing adaptive time stepping, get new time step
  % size and time step scaling
else
  accept = true;
  dtScale = 1;
  selfDist = [];
  wallDist = [];
  collUprate = 1;
  % if not doing adaptive time stepping, keep the same
  % time step size and always accept solution
end

if accept
  % take the new solution
  X = Xprov(:,:,end);
  sigma = sigmaProv(:,:,end);
  u = uProv(:,:,end);
  eta = etaProv(:,:,end);
  RS = RSprov(:,:,end);
else
  % revert to the old solution
  X = Xstore;
  sigma = sigStore;
  u = uStore;
  eta = etaStore;
  RS = RSstore;
end

end % timeStepGL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [accept,dtScale] = newTimeStepSize(o,...
      aFinal,lFinal,aInit,lInit,accept,om)
% [accept,dtScale] = newTimeStepSize(...
%      aFinal,lFinal,aInit,lInit,accept,om)
% finds a new time step size based on the change in the area and length
% from the first to final Gauss- Lobatto points.  The output is a flag
% indicating acceptance or rejection, and the amount the time step is
% scaled.  To increase the likelihood that the a new time step size is
% accepted, the time step size is never increased if the previous time
% step was rejected

alpha = o.alpha;
% buffer so that we aren't trying to keep error exactly at 1.  This is a
% safeguard so that the next time step is accepted with higher
% probability
betaUp = o.betaUp;
% allowable upscaling change in time step size
betaDown = o.betaDown;
% allowable downscaling change in time step size

errArea = abs(aInit - aFinal);
errLength = abs(lInit - lFinal);
% absolute errors in area and length
tauArea = max(aInit,aFinal)*o.rtolArea*o.dt;
tauLength = max(lInit,lFinal)*o.rtolLength*o.dt;
% Tolerance for errArea and errLength
err = max(max(errArea./tauArea),max(errLength./tauLength));
% Maximum of relative errors in area and length
% Want this quantity to be as close to 1 as possible

actualOrder = o.expectedOrder;
% order of time stepping method
dtOPT = err^(-1/(actualOrder))*o.dt;
% optimal time step size

dtOld = o.dt;
if accept
  o.dt = alpha^(1/actualOrder) * ...
      min(betaUp*o.dt,max(dtOPT,betaDown*o.dt));
else
  o.dt = alpha^(1/(actualOrder)) * ...
      min(o.dt,max(dtOPT,betaDown*o.dt));
  % don't want to scale up this time step if it was previously rejected.
  % This hopefully gets rid of pattern where a solution alternates
  % between being accepted and rejected.
end
% safety factor added to the optimal time step size also, time step size
% is not scaled up or down too fast For safety factor, take 1/p root.
% In our formulation, this makes alpha the desired value for err
% regardless of the order of the method.
dtScale = o.dt/dtOld;
% time step scaling

if err > 1
  accept = false;
  % reject time step because the error is too large
  message = ['Time Step REJECTED with error ' ...
      num2str(err,'%4.2e')];
  om.writeMessage(message,'%s\n')
  message = ['Time step scaled by           ' ... 
      num2str(dtScale,'%4.2e')];
  om.writeMessage(message,'%s\n')
  message = ['New time step size is         ' ...
      num2str(o.dt,'%4.2e')];
  om.writeMessage(message,'%s\n')
  om.writeMessage(' ','%s\n')
else
  accept = true;
  % accept the solution because the error is small
  message = ['Time Step ACCEPTED with error ' ...
      num2str(err,'%4.2e')];
  om.writeMessage(message,'%s\n')
  message = ['Time step scaled by           ' ... 
      num2str(dtScale,'%4.2e')];
  om.writeMessage(message,'%s\n')
  message = ['New time step size is         ' ...
      num2str(o.dt,'%4.2e')];
  om.writeMessage(message,'%s\n')
end

end % newTimeStepSize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,sigma,u,eta,RS,iter,iflag] = timeStep(o,...
    Xstore,sigStore,uStore,etaStore,RSstore,...
    deltaX,deltaSig,deltaEta,deltaRS,...
    diffResidual,vesVel,wallVel,kappa,viscCont,walls,wallsCoarse,...
    updatePreco,vesicle,sa,IK)
% [X,sigma,u,eta,RS,iter,iflag] = timeStep(o,...
%    Xstore,sigStore,uStore,etaStore,RSstore,...
%    deltaX,deltaSig,deltaEta,deltaRS,...
%    diffResidual,vesVel,wallVel,kappa,viscCont,walls,wallsCoarse,...
%    updatePreco,sa,IK)
% uses explicit or implicit vesicle-vesicle interactions and
% discretizes the inextensibility condition in three different way
% (method1, method2, or method 3).  Must pass in the vesicle positions,
% tension, and velocity from enough previous time steps (depends on
% o.order).  Returns a new positions, tension, velocity, density
% function defined on the solid walls and the number of required GMRES
% iterations if o.SDCcorrect=true, then it uses deltaX, deltaSig, etc
% to compute the right-hand sides needed for sdc updates updatePreco is
% a flog that decides if the block-diagonal preconditioner should be
% updated or not
% NOTE THAT A LOT OF THE FEATURES ARE NOT IMPLEMENTED WITH EXPLICIT
% VESICLE-VESICLE INTERACTIONS.

N = size(Xstore,1)/2; % Number of points per vesicle
nv = size(Xstore,2); % Number of vesicles
if o.confined
  Xwalls = walls.X; % discretization points of solid walls
else
  Xwalls = [];
end
Nbd = size(Xwalls,1)/2; % Number of points on the solid walls
nvbd = size(Xwalls,2); % number of solid wall components
alpha = (1 + viscCont)/2; 
% constant that appears in front of time derivative in
% vesicle dynamical equations


Xm = zeros(2*N,nv);
sigmaM = zeros(N,nv);
uM = zeros(2*N,nv);
Xo = zeros(2*N,nv);
etaM = zeros(2*Nbd,nvbd);
RSm = zeros(3,nvbd);
for k = 1:o.order
  Xm = Xm + Xstore(:,:,k)*o.Xcoeff(k);
  sigmaM = sigmaM + sigStore(:,:,k)*o.Xcoeff(k);
  uM = uM + uStore(:,:,k)*o.Xcoeff(k);
  if o.confined
    etaM = etaM + etaStore(:,:,k)*o.Xcoeff(k);
    RSm = RSm + RSstore(:,:,k)*o.Xcoeff(k);
  end
  Xo = Xo + Xstore(:,:,k)*o.rhsCoeff(k);
end
% Form linear combinations of previous time steps needed for Ascher,
% Ruuth, and Wetton IMEX methods


if o.order ~= 1
  vesicle = capsules(Xm,sigmaM,uM,kappa,viscCont);
end

% build an object vesicle that contains tangent vector, jacobian, etc.

op = o.op;
if ~o.SDCcorrect
  o.Galpert = op.stokesSLmatrix(vesicle);
end
% Build single layer potential matrix and put it in current object
% If we are doing an sdc update, this is already precomputed and 
% stored from when we formed the provisional solution

if ~o.SDCcorrect
  if any(viscCont ~= 1)
    o.D = op.stokesDLmatrix(vesicle);
  else
    o.D = [];
  end
end
% Compute double-layer potential matrix due to each vesicle
% independent of the others.  Matrix is zero if there is no
% viscosity contrast

if ~o.SDCcorrect
  if o.near
    if o.confined
      [o.NearV2V,o.NearV2W] = vesicle.getZone(walls,3);

      % Need vesicle to vesicle and vesicle to wall interactions
      if nvbd == 1 
        [o.NearW2W,o.NearW2V] = walls.getZone(vesicle,2);
      else
        [o.NearW2W,o.NearW2V] = walls.getZone(vesicle,3);
      end
      % Only need wall to vesicle interactions.  Wall to wall
      % interactions should also use near-singular integration since
      % they may be close to one another
    else
      o.NearV2V = vesicle.getZone([],1);
      % no solid walls, so only need vesicle-vesicle intearactions
      o.NearV2W = [];
      o.NearW2V = [];
      o.NearW2W = [];
    end
  else
    o.NearV2V = [];
    o.NearV2W = [];
    o.NearW2V = [];
    o.NearW2W = [];
  end
  % If using near-singular integration, need structures for deciding who
  % is close, how close it is, who is closest, etc.
end
% Only form near-singular integration structure if not doing an SDC
% update.  Otherwise this was formed and saved when forming the
% provisional solution

if o.collDete
  o.lapDLP = op.laplaceDLmatrix(vesicle);
  % build matrix that evaluates double-layer laplace potential for each
  % vesicle

  oc = curve;
  [icollisionVes,icollisionWall] = ...
    oc.collision(vesicle,walls,o.NearV2V,o.NearV2W,o.fmm,o.near);
  % Check for collisions 
  if icollisionVes
    fprintf('VESICLES HAVE CROSSED\n')
    pause
  end
  if icollisionWall
    fprintf('VESICLES HAVE CROSSED SOLID WALL\n')
    pause
  end
end
% check for collisions

if ~o.SDCcorrect
  rhs1 = Xo;
  rhs2 = zeros(N,nv);
  if o.confined
    rhs3 = walls.u;
  else
    rhs3 = [];
  end
else
  rhs1 = deltaX + diffResidual;
  if any(vesicle.viscCont ~= 1)
    z = zeros(2*N*nv,1);
    for k = 1:nv
      z(2*(k-1)*N+1:2*k*N) = rhs1(:,k);
    end
    z = o.IminusD(z,vesicle);
    for k = 1:nv
      rhs1(:,k) = z(2*(k-1)*N+1:2*k*N)/alpha(k);
    end
  end
  if strcmp(o.solver,'method1')
    rhs2 = ones(N,nv);
  else
    rhs2 = -vesicle.surfaceDiv(vesVel);
  end
  if o.confined
    rhs3 = -wallVel + walls.u;
  else
    rhs3 = [];
  end
end
% Parts of rhs from previous solution.  The right-hand-side depends on
% whether we are doing an SDC correction or forming the provisional
% solution.

% START TO COMPUTE RIGHT-HAND SIDE DUE TO VESICLE TRACTION JUMP
if strcmp(o.vesves,'explicit')
  if o.profile
    tic
  end
  if ~o.SDCcorrect
    f = vesicle.tracJump(Xm,sigmaM);
  else
    f = vesicle.tracJump(deltaX,deltaSig);
  end
  if ~o.fmm
    kernel = @op.exactStokesSL;
    kernelDirect = @op.exactStokesSL;
  else
    kernel = @op.exactStokesSLfmm;
    kernelDirect = @op.exactStokesSL;
  end

  if ~o.near
    Fslp = kernel(vesicle,f,[]);
    % Evaulate single-layer potential on all vesicles but itself without
    % near-singular integration
    if o.confined
      [~,FSLPwall] = kernel(vesicle,f,[],walls.X,(1:nv));
      % Evaluate single-layer potential on solid walls due to all
      % vesicles
    else
      FSLPwall = [];
    end

  else
    SLP = @(X) op.exactStokesSLdiag(vesicle,o.Galpert,X);
    Fslp = op.nearSingInt(vesicle,f,SLP,...
        o.NearV2V,kernel,kernelDirect,vesicle,true,false);
    % Use near-singular integration to compute single-layer potential
    % due to all other vesicles.  Need to pass function
    % op.exactStokesSL so that the layer potential can be computed at
    % far points and Lagrange interpolation points

    if o.confined
      FSLPwall = op.nearSingInt(vesicle,f,SLP,...
          o.NearV2W,kernel,kernelDirect,walls,false,false);
      % Evaluate the velocity on the walls due to the vesicles
    else
      FSLPwall = [];
    end

  end
  if o.profile
    fprintf('Build right-hand side V2V and V2W   %5.1e\n',toc);
  end
else
  Fslp = zeros(2*N,nv);
  FSLPwall = zeros(2*Nbd,nvbd);
  % vesicle-vesicle and vesicle-wall interactions are handled
  % implicitly in TimeMatVec
end
for k = 1:nv
  rhs1(:,k) = rhs1(:,k) + o.dt/alpha(k)*Fslp(:,k);
end

rhs3 = rhs3 - FSLPwall;
% END TO COMPUTE RIGHT-HAND SIDE DUE TO VESICLE TRACTION JUMP

% START TO COMPUTE RIGHT-HAND SIDE DUE TO VISCOSITY CONTRAST
if any(vesicle.viscCont ~= 1)
  if o.near
    jump = 1/2*(1-vesicle.viscCont);
    DLP = @(X) X*diag(jump) + op.exactStokesDLdiag(vesicle,o.D,X);
  end
  % Need to add jump to double-layer potential if using near-singular
  % integration so that we can compute boundary values for the
  % near-singular integration algorithm

  if ~o.fmmDLP
    kernel = @op.exactStokesDL;
    kernelDirect = @op.exactStokesDL;
  else
    kernel = @op.exactStokesDLnewfmm;
    kernelDirect = @op.exactStokesDL;
  end
  if strcmp(o.vesves,'implicit')
    if ~o.SDCcorrect
      density = Xo;
    else
      density = deltaX;
    end

    if ~o.near
      Fdlp = kernel(vesicle,density);
      if o.confined
        [~,FDLPwall] = kernel(vesicle,density,walls.X,(1:nv));
      else
        FDLPwall = [];
      end
    else
      Fdlp = op.nearSingInt(vesicle,density,DLP,...
          o.NearV2V,kernel,kernelDirect,vesicle,true,false);
      % Use near-singular integration to compute double-layer
      % potential from previous solution
      if o.confined
        FDLPwall = op.nearSingInt(vesicle,density,DLP,...
          o.NearV2W,kernel,kernelDirect,walls,false,false);
      else
        FDLPwall = [];
      end
    end
  elseif strcmp(o.vesves,'explicit')
    if ~o.near
      Fdlp = -o.dt * kernel(vesicle,uM);
      % Evaulate the velocity due to the viscosity contrast
      % on all vesicles but itself WITHOUT near-singular
      % integration
      if o.confined
        [~,FDLPwall] = kernel(vesicle,uM,walls.X,(1:nv));
        FDLPwall = -o.dt*FDLPwall;
        % Evaulate the velocity due to the viscosity contrast on the
        % walls WITHOUT near-singulation integration
      else
        FDLPwall = [];
      end
    else
      kernelDirect = @op.exactStokesDL;
      Fdlp = -o.dt * op.nearSingInt(vesicle,uM,DLP,...
          o.NearV2V,kernel,kernelDirect,vesicle,true,false);
      % Evaulate the velocity due to the viscosity contrast on all
      % vesicles but itself WITH near-singular integration

      if o.confined
        FDLPwall = -o.dt * op.nearSingInt(vesicle,uM,DLP,...
            o.NearV2W,kernel,kernelDirect,walls,false,false);
        % Evaulate the velocity due to the viscosity contrast on the
        % walls WITH near-singulation integration
      else
        FDLPwall = [];
      end
    end
  end
else
  Fdlp = zeros(2*N,nv);
  FDLPwall = zeros(2*Nbd,nvbd);
  % If no viscosity contrast, there is no velocity induced due to a
  % viscosity contrast
end

if (any(viscCont ~= 1) && ~o.SDCcorrect)
  DXo = op.exactStokesDLdiag(vesicle,o.D,Xo);
  rhs1 = rhs1 - (Fdlp + DXo) * diag(1./alpha);
end
% add in viscosity contrast term due to each vesicle independent of the
% others (o.D * Xo) from the previous solution followed by the term due
% to all other vesicles (Fdlp)

rhs3 = rhs3 + FDLPwall/o.dt;
% compute the double-layer potential due to all other vesicles from the
% appropriate linear combination of previous time steps.  Depends on
% time stepping order and vesicle-vesicle discretization

if o.adhesion 
  adhesion = vesicle.adhesionTerm(o.adStrength,o.adRange);
  if ~o.fmm
    kernel = @op.exactStokesSL;
    kernelDirect = @op.exactStokesSL;
  else
    kernel = @op.exactStokesSLfmm;
    kernelDirect = @op.exactStokesSL;
  end

  Fadhesion = op.exactStokesSLdiag(vesicle,o.Galpert,adhesion);
  % diagonal term of adhesion

  if ~o.near
    Fadhesion = Fadhesion + kernel(vesicle,adhesion);
    % Evaulate single-layer potential on all vesicles but itself without
    % near-singular integration
  else
    SLP = @(X) op.exactStokesSLdiag(vesicle,o.Galpert,X);
    Fadhesion = Fadhesion + ...
        op.nearSingInt(vesicle,adhesion,SLP,...
        o.NearV2V,kernel,kernelDirect,vesicle,true,false);
    % Use near-singular integration to compute single-layer potential
    % due to all other vesicles.  Need to pass function
    % op.exactStokesSL so that the layer potential can be computed at
    % far points and Lagrange interpolation points
  end

  rhs1 = rhs1 + o.dt*Fadhesion*diag(1./alpha);
end
% Compute velocity due to adhesion


% START TO COMPUTE RIGHT-HAND SIDE DUE TO SOLID WALLS
if o.confined
  if strcmp(o.vesves,'explicit')
    if o.profile
      tic
    end
    if o.near
      jump = -1/2;
      DLP = @(X) jump*X + op.exactStokesDLdiag(walls,o.wallDLP,X);
    end
    % compute the matrix for doing evaluating the double-layer potential
    % on the solid walls required for near-singular integration

    if ~o.fmmDLP
      kernel = @op.exactStokesDL;
      kernelDirect = @op.exactStokesDL;
    else
      kernel = @op.exactStokesDLnewfmm;
      kernelDirect = @op.exactStokesDL;
    end
    if ~o.near
      if ~o.SDCcorrect
        charge = etaM;
      else
        charge = deltaEta;
      end
      [~,Fwall2Ves] = kernel(walls,charge,...
          vesicle.X,1:nvbd);
      % velocity field due to the walls evaluated on the vesicle
    else
      if ~o.SDCcorrect
        charge = etaM;
      else
        charge = deltaEta;
      end
      Fwall2Ves = op.nearSingInt(walls,charge,DLP,...
          o.NearW2V,kernel,kernelDirect,vesicle,false,false);
    end

    for k = 2:nvbd
      if ~o.SDCcorrect
        stokeslet = RSm(1:2,k);
        rotlet = RSm(3,k);
      else
        stokeslet = deltaRS(1:2,k);
        rotlet = deltaRS(3,k);
      end
      Fwall2Ves = Fwall2Ves + ...
        o.RSlets(vesicle.X,walls.center(:,k),stokeslet,rotlet);
    end
    if o.profile
      fprintf('Build right-hand side W2V           %5.1e\n',toc);
    end

  elseif strcmp(o.vesves,'implicit')
    Fwall2Ves = zeros(2*N,nv);
  end
else
  Fwall2Ves = zeros(2*N,nv);
  if ~o.SDCcorrect
    rhs1 = rhs1 + o.dt*o.farField(Xm)*diag(1./alpha);
    % Add in far-field condition (extensional, shear, etc.)
  end
end
rhs1 = rhs1 + o.dt*Fwall2Ves*diag(1./alpha);
% right-hand side of the velocity evaluated on the solid walls
% END TO COMPUTE RIGHT-HAND SIDE DUE TO SOLID WALLS


% START TO COMPUTE THE RIGHT-HAND SIDE FOR THE INEXTENSIBILITY CONDITION
if (strcmp(o.solver,'method1'))
  % If using method1, vesicle-vesicle interactions and the presence of a
  % viscosity contrast is irrelevent
  if ~o.SDCcorrect
    rhs2 = rhs2 + vesicle.surfaceDiv(Xo); 
  else
    divf = zeros(N,nv);
    for k = 1:nv
      divf(:,k) = curve.arcDeriv(Xm(1:N,k),1,1./sa(:,k),...
          IK(:,k)).^2 + ...
                  curve.arcDeriv(Xm(N+1:2*N,k),1,1./sa(:,k),...
          IK(:,k)).^2;
    end
    rhs2 = 1/2*(rhs2 - divf);
  end
else 
  % If using method2, method3, or method4, vesicle-vesicle interaction
  % affects the right-hand side of the inextensibility condition

  if ~o.confined && ~o.SDCcorrect
    if any(viscCont ~= 1)
      if o.profile
        tic
      end
      rhs2 = rhs2 - vesicle.surfaceDiv(...
          o.solveIminusD(o.farField(Xm),vesicle));
      if o.profile
        fprintf('Solve system alpha*I - DLP          %5.1e\n',toc);
      end
    else
      rhs2 = rhs2 - vesicle.surfaceDiv(o.farField(Xm));
    end
  end
  % add in term from farfield
  
  if (strcmp(o.vesves,'explicit'))
    rhs2 = rhs2 - vesicle.surfaceDiv(Fslp); 
    rhs2 = rhs2 - vesicle.surfaceDiv(Fwall2Ves);
    if o.adhesion
      rhs2 = rhs2 + vesicle.surfaceDiv(Fadhesion);
    end
  end
end
% rhs2 is the right-hand side for the inextensibility condition
% END TO COMPUTE THE RIGHT-HAND SIDE FOR THE INEXTENSIBILITY CONDITION

if (any(vesicle.viscCont ~= 1) && ...
      strcmp(o.vesves,'implicit') && o.confined)
  rhs3 = rhs3 * o.dt;
end
% This makes sure that the rhs is all order one rather than have rhs3
% being order 1/o.dt and other two parts (rhs1 and rhs2) being order 1.
% This of course needs to be compensated in the TimeMatVec routine

rhs = [rhs1; rhs2];
rhs = rhs(:);
rhs = [rhs; rhs3(:)];
% Stack the right-hand sides in an alternating with respect to the
% vesicle fashion
% Add on the no-slip boundary conditions on the solid walls
rhs = [rhs; zeros(3*(nvbd-1),1)];
% Rotlet and Stokeslet equations

usePreco = true;
% usePreco = false;


% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% BUILD THE WHOLE VESICLE-VESICLE INTERACTION MATRIX
% This is just to look at the properties of the whole matrix and test the
% preconditioner, and investigate if we improve the preconditioner.

% if 1
% [Ben,Ten,Div] = vesicle.computeDerivs;
% INTER = zeros(3*N*nv,3*N*nv);
% PRECO = zeros(3*N*nv,3*N*nv);
% if any(vesicle.viscCont~=1)
%   for k1 = 1:nv
%     for k2 = 1 : nv
%       is = 3*N*(k1-1)+1;
%       ie = 3*N*(k1);
%       js = 3*N*(k2-1)+1;
%       je = 3*N*(k2);
%       if k2 == k1
%         prevTerms = o.beta*eye(2*N);
%         extenTerm = o.beta*Div(:,:,k1);
%         SLmat = o.Galpert(:,:,k1);
%         DLmat = o.D(:,:,k1);
%         PRECO(is:ie,js:je) =  [o.beta*(eye(2*N) - o.D(:,:,k1)/alpha(k1)) + ...
%             o.dt/alpha(k1)*vesicle.kappa*o.Galpert(:,:,k1)*Ben(:,:,k1) ...
%         -o.dt/alpha(k1)*o.Galpert(:,:,k1)*Ten(:,:,k1); ...
%         o.beta*Div(:,:,k1) zeros(N)];
%       else
%         prevTerms = zeros(2*N);  
%         extenTerm = zeros(N,2*N);
%         SLmat = o.op.stokesSLinterMatrix(vesicle,k1,k2);
%         DLmat = o.op.stokesDLinterMatrix(vesicle,k1,k2);
%       end
% 
%       INTER(is:ie,js:je) = [(prevTerms - o.beta*DLmat/alpha(k1)) + ...
%             o.dt/alpha(k1)*vesicle.kappa*SLmat*Ben(:,:,k2) ...
%         -o.dt/alpha(k1)*SLmat*Ten(:,:,k2); ...
%         extenTerm zeros(N)];
%     end
%   end
% else
%   for k1 = 1:nv
%     for k2 = 1 : nv
%       is = 3*N*(k1-1)+1;
%       ie = 3*N*(k1);
%       js = 3*N*(k2-1)+1;
%       je = 3*N*(k2);
% 
%       if k2 == k1
%         prevTerms = o.beta*eye(2*N);  
%         extenTerm = o.beta*Div(:,:,k1);
%         SLmat = o.Galpert(:,:,k1);
%         PRECO(is:ie,js:je) = [o.beta*eye(2*N) + ...
%             o.dt*vesicle.kappa*o.Galpert(:,:,k1)*Ben(:,:,k1) ...
%         -o.dt/alpha(k1)*o.Galpert(:,:,k1)*Ten(:,:,k1); ...
%         o.beta*Div(:,:,k1) zeros(N)];
%       else
%         prevTerms = zeros(2*N);  
%         extenTerm = zeros(N,2*N);
%         SLmat = o.op.stokesSLinterMatrix(vesicle,k1,k2);
%       end
% 
%       INTER(is:ie,js:je) = [prevTerms + ...
%             o.dt*vesicle.kappa*SLmat*Ben(:,:,k2) ...
%         -o.dt/alpha(k1)*SLmat*Ten(:,:,k2); ...
%         extenTerm zeros(N)];
% 
%     end
%   end
% end
% 
% save TGinterMatrixVC1 INTER PRECO
% disp('Interaction matrix is saved')
% pause
% end
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
% START BUILDING BLOCK-DIAGONAL PRECONDITIONER
if usePreco
  if updatePreco
    % only build the preconditioner if updatePreco == true
    if o.profile
      tic
    end
    [Ben,Ten,Div] = vesicle.computeDerivs;
    if o.profile
      fprintf('Build differential operators        %5.1e\n',toc);
    end
    % Compute bending, tension, and surface divergence of current
    % vesicle configuration

    if o.profile
      tic
    end

    if (strcmp(o.solver,'method1') || strcmp(o.solver,'method2'))
      bdiagVes.L = zeros(3*N,3*N,nv);
      bdiagVes.U = zeros(3*N,3*N,nv);
    elseif strcmp(o.solver,'method3')
      bdiagVes.GT = zeros(2*N,N,nv); % SLP of tension
      bdiagVes.DGT = zeros(N,N,nv); % divergence of SLP of tension
      bdiagVes.DGB = zeros(N,2*N,nv); % divergence of SLP of bending
      bdiagVes.schur = zeros(2*N,2*N,nv); 
      % schur complement of the lower right N by N block
    elseif strcmp(o.solver,'method4')
      bdiagVes.GT = zeros(2*N,N,nv); % SLP of tension
      bdiagVes.IpBen = zeros(2*N,2*N,nv); 
      % inverse of identity plus SLP of Bending
      bdiagVes.DGB = zeros(N,2*N,nv); % divergence of SLP of bending
      bdiagVes.schur = zeros(N,N,nv); 
      % schur complement of the upper left 2*N by 2*N block
    end
    
    
    for k=1:nv
      if strcmp(o.solver,'method1')
        if any(vesicle.viscCont ~= 1)
          [bdiagVes.L(:,:,k),bdiagVes.U(:,:,k)] = lu(...
            [o.beta*(eye(2*N) - o.D(:,:,k)/alpha(k)) + ...
                o.dt/alpha(k)*vesicle.kappa*o.Galpert(:,:,k)*Ben(:,:,k) ...
            -o.dt/alpha(k)*o.Galpert(:,:,k)*Ten(:,:,k); ...
            o.beta*Div(:,:,k) zeros(N)]);
        else
          [bdiagVes.L(:,:,k),bdiagVes.U(:,:,k)] = lu(...
            [o.beta*eye(2*N) + ...
                o.dt*vesicle.kappa*o.Galpert(:,:,k)*Ben(:,:,k) ...
            -o.dt/alpha(k)*o.Galpert(:,:,k)*Ten(:,:,k); ...
            o.beta*Div(:,:,k) zeros(N)]);
        end
      elseif strcmp(o.solver,'method2')
        if any(vesicle.viscCont ~= 1)
          [bdiagVes.L(:,:,k),bdiagVes.U(:,:,k)] = lu( ...
            [o.beta*(eye(2*N) - o.D(:,:,k)/alpha(k)) + ...
                o.dt/alpha(k)*vesicle.kappa*o.Galpert(:,:,k)*Ben(:,:,k) ...
            -o.dt/alpha(k)*o.Galpert(:,:,k)*Ten(:,:,k); ...
            -vesicle.kappa*Div(:,:,k)*(alpha(k)*eye(2*N) - o.D(:,:,k))\...
                (o.Galpert(:,:,k)*Ben(:,:,k)) ...
            Div(:,:,k)*(alpha(k)*eye(2*N) - o.D(:,:,k))\...
                (o.Galpert(:,:,k)*Ten(:,:,k))]);
        else
          [bdiagVes.L(:,:,k),bdiagVes.U(:,:,k)] = lu(...
            [o.beta*eye(2*N) + ...
                o.dt*vesicle.kappa*o.Galpert(:,:,k)*Ben(:,:,k) ...
            -o.dt*o.Galpert(:,:,k)*Ten(:,:,k); ...
            -vesicle.kappa*Div(:,:,k)*o.Galpert(:,:,k)*Ben(:,:,k) ...
            Div(:,:,k)*o.Galpert(:,:,k)*Ten(:,:,k)]);
        end
      elseif strcmp(o.solver,'method3')
        % schur complement of lower right block of method2
        bdiagVes.GT(:,:,k) = o.Galpert(:,:,k)*Ten(:,:,k);
        bdiagVes.DGT(:,:,k) = (Div(:,:,k)*o.Galpert(:,:,k)*...
            Ten(:,:,k))\eye(N);
        bdiagVes.DGB(:,:,k) = ...
          vesicle.kappa*Div(:,:,k)*o.Galpert(:,:,k)*Ben(:,:,k);
        bdiagVes.schur(:,:,k) = ...
          inv((o.beta*eye(2*N) + vesicle.kappa*o.dt*...
          o.Galpert(:,:,k)*Ben(:,:,k)) - ...
          o.dt*o.Galpert(:,:,k)*Ten(:,:,k)*bdiagVes.DGT(:,:,k)*...
          vesicle.kappa*Div(:,:,k)*o.Galpert(:,:,k)*Ben(:,:,k));

      elseif strcmp(o.solver,'method4')
        % schur complement of upper left block of method2
        bdiagVes.GT(:,:,k) = o.Galpert(:,:,k)*Ten(:,:,k);
        bdiagVes.IpBen(:,:,k) = inv(o.beta*eye(2*N) + ...
            vesicle.kappa*o.dt*o.Galpert(:,:,k)*Ben(:,:,k));
        bdiagVes.DGB(:,:,k) = ...
          vesicle.kappa*Div(:,:,k)*o.Galpert(:,:,k)*Ben(:,:,k);
        bdiagVes.schur(:,:,k) = ...
          inv(Div(:,:,k)*(-vesicle.kappa*o.dt*o.Galpert(:,:,k)*...
          Ben(:,:,k)*bdiagVes.IpBen(:,:,k) + eye(2*N))*...
          bdiagVes.GT(:,:,k));
      end
    end
    o.bdiagVes = bdiagVes;
    % Build block-diagonal preconditioner of self-vesicle 
    % intearctions in matrix form
    if o.profile
      fprintf('Build block-diagonal preconditioner %5.1e\n\n',toc);
    end
  end % updatePreco
end % usePreco

if o.profile
  tGMRES = tic;
end

% USE LU DECOMPOSITION
if 0
  Xn = o.preconditionerBD(rhs);
  iflag = 0;
  R = 0;
  I = [1;1];
  resvec = [0];
  iter = I(2);
end
 

if 1
  warning off
  % any warning is printed to the terminal and the log file so
  % don't need the native matlab version
  initGMRES = [Xm;sigmaM];
  initGMRES = initGMRES(:);
  if o.confined
    RS = RSm(:,2:end);
    initGMRES = [initGMRES;etaM(:);RS(:)];
  end
  
  if usePreco
    % when doing corrections, expect that these guys will be small
    % GMRES
    
    [Xn,iflag,R,I,resvec] = gmres(@(X) o.TimeMatVec(X,vesicle,walls),...
        rhs,[],o.gmresTol,o.gmresMaxIter,...
        @o.preconditionerBD,[],initGMRES);
    iter = I(2);
  else
    [Xn,iflag,R,I,resvec] = gmres(@(X) o.TimeMatVec(X,vesicle,walls),...
        rhs,[],o.gmresTol,o.gmresMaxIter);
    iter = I(2);
  end
  warning on
end
% Use GMRES to solve for new positions, tension, density
% function defined on the solid walls, and rotlets/stokeslets


if o.profile
  fprintf('Time to do one time step is         %5.1e\n\n',toc(tGMRES))
end

% iter = I(2);
% Number of GMRES iterations

X = zeros(2*N,nv);
sigma = zeros(N,nv);
eta = zeros(2*Nbd,nvbd);
RS = zeros(3,nvbd);

for k=1:nv
  X(:,k) = Xn((3*k-3)*N+1:(3*k-1)*N);
  sigma(:,k) = Xn((3*k-1)*N+1:3*k*N);
end
% unstack the positions and tensions

% Unstack the positions and tensions
Xn = Xn(3*nv*N+1:end);
for k = 1:nvbd
  eta(:,k) = Xn((k-1)*2*Nbd+1:2*k*Nbd);
end
% unstack the density function
otlets = Xn(2*nvbd*Nbd+1:end);
for k = 2:nvbd
  istart = (k-2)*3+1;
  iend = 3*(k-1);
  RS(:,k) = otlets(istart:iend);
end
% unstack the rotlets and stokeslets

u = (o.beta*X - Xo)/o.dt;
% Compute the velocity using the differencing stencil

end % timeStep

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = TimeMatVec(o,Xn,vesicle,walls)
% val = TimeMatVec(Xn,vesicle,walls) 
% MATVEC for GMRES in the IMEX scheme.
% Evaluations vesicle-vesicle and vesicle-boundary interaction formulas
% to the function Xn which contains both the position, tension, and
% density
% 
% - Xn : all state variables in the following order
%   ves1:x,y,sigma,  ves2:x,y,sigma, ... vesNv:x,y,sigma, (outer_solidwall1:fx,fy,
%   inner_solidwall_1:fx,fy; inner_solid_wall2:fx,fy; ...; inner_solid_walln:fx,fy;
%   stokeslet_rotlet_innerwall1, stokeslet_rolet_innerwall2....
%
% - vesicle: class capsules used to evaluate the operators for the GMRES 
% - walls: same thing as "vesicle" but for the confined walls geometry


global matvecs  
matvecs = matvecs + 1;
% counter for the number of matrix-vector multiplications
% that are required for the entire simulation

op = o.op; % poten class
N = vesicle.N; % Number of points
nv = vesicle.nv; % Number of vesicles
if o.confined
  Nbd = walls.N; % Number of points on walls
  nvbd = walls.nv; % Number of components to walls
else
  Nbd = 0;
  nvbd = 0;
end

valPos = zeros(2*N,nv);
% right-hand side that corresponds to position equation
valTen = zeros(N,nv);
% right-hand side that corresponds to inextensibilty equation
valWalls = zeros(2*Nbd,nvbd);
% right-hand side that corresponds to solid wall equation
valLets = zeros(3*(nvbd-1),1);
% right-hand side corresponding to the rotlets and stokeslets

Xm = zeros(2*N,nv);
sigmaM = zeros(N,nv);
for k=1:nv
  Xm(1:2*N,k) = Xn((3*k-3)*N+1:(3*k-1)*N);
  sigmaM(:,k) = Xn((3*k-1)*N+1:3*k*N);
end
% Unstack the position and tension from the input

eta = Xn(3*nv*N+1:end);
etaM = zeros(2*Nbd,nvbd);
for k = 1:nvbd
  etaM(:,k) = eta((k-1)*2*Nbd+1:2*k*Nbd);
end
% Unstack the density function from the input
otlets = Xn(3*nv*N+2*nvbd*Nbd+1:end);
% stokeslets and rotlets of each vesicle.  Ordered as
% [stokeslet1(component 1);stokeslet1(component 2);rotlet1;...
%  stokeslet2(component 1);stokeslet2(component 2);rotlet2;...];

f = vesicle.tracJump(Xm,sigmaM);
% f is the traction jump stored as a 2N x nv matrix

alpha = (1+vesicle.viscCont)/2; 
% constant that multiplies the time derivative in the 
% vesicle position equation

Gf = op.exactStokesSLdiag(vesicle,o.Galpert,f);
% Gf is the single-layer potential applied to the traction jump. 

if any(vesicle.viscCont ~= 1)
  DXm = op.exactStokesDLdiag(vesicle,o.D,Xm);
else
  DXm = zeros(2*N,nv);
end
% DXm is the double-layer potential applied to the position

% START COMPUTING REQUIRED SINGLE-LAYER POTENTIALS
if strcmp(o.vesves,'implicit')
  if o.profile
    tic
  end
  if ~o.fmm
    kernel = @op.exactStokesSL;
    kernelDirect = @op.exactStokesSL;
  else
    kernel = @op.exactStokesSLfmm;
    kernelDirect = @op.exactStokesSL;
  end

  if ~o.near
    Fslp = kernel(vesicle,f,[]);
    % Evaulate single-layer potential due to all other vesicles
    % WITHOUT near singular integration.  FMM is optional
    if o.confined
      [~,FSLPwall] = kernel(vesicle,f,[],walls.X,(1:nv));
      % Evaluate single-layer potential due to all vesicles on
      % the solid walls WITHOUT near-singular integration
    else
      FSLPwall = [];
    end
  else
    SLP = @(X) op.exactStokesSLdiag(vesicle,o.Galpert,X);
    Fslp = op.nearSingInt(vesicle,f,SLP,...
        o.NearV2V,kernel,kernelDirect,vesicle,true,false);
    % Evaulate single-layer potential due to all other vesicles
    % WITH near-singular integration.  FMM is optional
    
    if o.confined
      FSLPwall = op.nearSingInt(vesicle,f,SLP,...
        o.NearV2W,kernel,kernelDirect,walls,false,false);
      % Evaluate single-layer potential due to all vesicles on
      % the solid walls WITH near-singular integration
    else
      FSLPwall = [];
    end
  end
  if o.profile
    fprintf('Apply V2V and V2W interactions      %5.1e\n',toc) 
  end
else
  Fslp = zeros(2*N,nv);
  FSLPwall = zeros(2*Nbd,nvbd); 
  % These terms is handled explicitly if o.vesves is 'explicit'
end
% END COMPUTING REQUIRED SINGLE-LAYER POTENTIALS
% Evaluate single-layer potential due to all vesicles except itself and
% the single-layer potential due to all vesicles evaluated on the solid
% walls.  Sets the layer-potential to zero if using explicit
% interactions


% START COMPUTING REQUIRED DOUBLE-LAYER POTENTIALS FOR VISCOSITY
% CONTRAST
if any(vesicle.viscCont ~= 1)
  if o.near
    jump = 1/2*(1-vesicle.viscCont);
    DLP = @(X) X*diag(jump) + op.exactStokesDLdiag(vesicle,o.D,X);
  end

  if strcmp(o.vesves,'implicit')
    if ~o.fmmDLP
      kernel = @op.exactStokesDL;
      kernelDirect = @op.exactStokesDL;
    else
      kernel = @op.exactStokesDLnewfmm;
      kernelDirect = @op.exactStokesDL;
    end
    if ~o.near
      Fdlp = kernel(vesicle,Xm,[]);
      if o.confined
        [~,FDLPwall] = kernel(vesicle,Xm,walls.X,(1:nv));
      else
        FDLPwall = [];
      end
    else
      Fdlp = op.nearSingInt(vesicle,Xm,DLP,...
          o.NearV2V,kernel,kernelDirect,vesicle,true,false);
      % Use near-singular integration to compute double-layer
      % potential 
      if o.confined
        FDLPwall = op.nearSingInt(vesicle,Xm,DLP,...
          o.NearV2W,kernel,kernelDirect,walls,false,false);
      else
        FDLPwall = [];
      end
    end
  elseif strcmp(o.vesves,'explicit')
    Fdlp = zeros(2*N,nv);
    FDLPwall = zeros(2*Nbd,nvbd);
  end
else
  Fdlp = zeros(2*N,nv);
  FDLPwall = zeros(2*Nbd,nvbd);
end
% END COMPUTING REQUIRED DOUBLE-LAYER POTENTIALS FOR VISCOSITY CONTRAST


% START OF EVALUATING DOUBLE-LAYER POTENTIALS DUE TO SOLID WALLS
if o.confined
  if strcmp(o.vesves,'implicit')
    potWall = o.opWall;
    if ~o.fmmDLP
      kernel = @potWall.exactStokesDL;
      kernelDirect = @potWall.exactStokesDL;
    else
      kernel = @potWall.exactStokesDLnewfmm;
      kernelDirect = @potWall.exactStokesDL;
    end
    if ~o.near
      [~,Fwall2Ves] = kernel(walls,etaM,[],vesicle.X,1:nvbd);
    else
      jump = -1/2;
      DLP = @(X) jump*X + potWall.exactStokesDLdiag(walls,o.wallDLP,X);
      Fwall2Ves = potWall.nearSingInt(walls,etaM,DLP,...
          o.NearW2V,kernel,kernelDirect,vesicle,false,false);
    end
%     disp('FIX HERE')

    if o.profile
      fprintf('Apply W2V interaction               %5.1e\n',toc) 
    end
    % compute the velocity on the vesicles due to the solid walls
  else
    Fwall2Ves = zeros(2*N,nv);
  end
else
  Fwall2Ves = zeros(2*N,nv);
end
% END OF EVALUATING DOUBLE-LAYER POTENTIALS DUE TO SOLID WALLS

% START OF EVALUATING WALL TO WALL INTERACTIONS
if (o.confined && nvbd > 1)  
  if o.profile
    tic
  end
  % only need to do wall to wall interactions if the domain is multiply
  % connected
  potWall = o.opWall;
  if ~o.fmmDLP
    kernel = @potWall.exactStokesDL;
    FDLPwall2wall = kernel(walls,etaM);
  else
    kernel = @potWall.exactStokesDLnewfmm;
    FDLPwall2wall = zeros(2*walls.N,walls.nv);
    for k = 1:nvbd
      isou = [(1:k-1) (k+1:nvbd)];
      [~,FDLPwall2wall(:,k)] = kernel(walls,etaM,walls.X(:,k),isou);
    end
    % Since the double-layer potential is still expensive even with the
    % FMM, this eliminates the need to do one big FMM followed by a
    % bunch of small ones to subtract off the self-interaction term
    % which is calculated using the precomputed matrix
  end

  if o.profile
    fprintf('Apply W2W interaction               %5.1e\n\n',toc)
  end
end
% END OF EVALUATING WALL TO WALL INTERACTIONS

% START OF EVALUATING POTENTIAL DUE TO STOKESLETS AND ROTLETS
if nvbd > 1
  LetsWalls = zeros(2*Nbd,nvbd);
  LetsVes = zeros(2*N,nv);
  for k = 2:nvbd
    stokeslet = otlets(3*(k-2)+1:3*(k-2)+2);
    rotlet = otlets(3*(k-1));
    if strcmp(o.vesves,'implicit')
      LetsVes = LetsVes + o.RSlets(vesicle.X,walls.center(:,k),...
          stokeslet,rotlet);
    end
    % compute velocity due to rotlets and stokeslets on the vesicles

    LetsWalls = LetsWalls + o.RSlets(walls.X,walls.center(:,k),...
        stokeslet,rotlet);
    % compute velocity due to rotlets and stokeslets on the solid walls
  end
  valLets = o.letsIntegrals(otlets,etaM,walls);
  % Integral constraints on the density function eta related
  % to the weights of the stokeslets and rotlets
else
  LetsVes = zeros(2*N,nv);
  LetsWalls = zeros(2*Nbd,nvbd);
  FDLPwall2wall = zeros(2*Nbd,nvbd);
end

% END OF EVALUATING POTENTIAL DUE TO STOKESLETS AND ROTLETS
% START OF EVALUATING VELOCITY ON VESICLES
valPos = valPos - o.dt*Gf*diag(1./alpha);
% self-bending and self-tension terms
valPos = valPos - o.beta*DXm*diag(1./alpha);
% self-viscosity contrast term
valPos = valPos - o.dt*Fslp*diag(1./alpha);
% single-layer potential due to all other vesicles
valPos = valPos - o.beta*Fdlp*diag(1./alpha);
% double-layer potential due to all other vesicles
valPos = valPos - o.dt*Fwall2Ves*diag(1./alpha);
% velocity due to solid walls evaluated on vesicles
valPos = valPos - o.dt*LetsVes*diag(1./alpha);
% velocity on vesicles due to the rotlets and stokeslets
% END OF EVALUATING VELOCITY ON VESICLES
%**************************************************************************
% Repulsion Magnitude Test
% relVel = (Gf+o.beta*DXm+Fslp+Fdlp+LetsVes)*diag(1./alpha);
% integVel = (relVel(1:end/2,1).^2+relVel(end/2+1:end,1).^2).^0.5;
% SurfVelInteg = sum(integVel)*2*pi/N;
% disp('Velocity Integral')
% SurfVelInteg
%**************************************************************************

% START OF EVALUATING VELOCITY ON WALLS
if o.confined
  potWall = o.opWall;
  valWalls = valWalls - 1/2*etaM + ...
    potWall.exactStokesDLdiag(walls,o.wallDLP,etaM);
  valWalls(:,1) = valWalls(:,1) + ...
    potWall.exactStokesN0diag(walls,o.wallN0,etaM(:,1));
end
% evaluate velocity on solid walls due to the density function.
% self solid wall interaction

valWalls = valWalls + FSLPwall;
% velocity on walls due to the vesicle traction jump
valWalls = valWalls + o.beta*FDLPwall/o.dt;
% velocity on walls due to the vesicle viscosity jump
valWalls = valWalls + FDLPwall2wall;
% velocity on walls due to all other walls
valWalls = valWalls + LetsWalls;
% velocity on walls due to the rotlets and stokeslets
% END OF EVALUATING VELOCITY ON WALLS
% START OF EVALUATING INEXTENSIBILITY CONDITION
if (strcmp(o.solver,'method1'))
  valTen = o.beta * vesicle.surfaceDiv(Xm);
  
  % compute surface divergence of the current GMRES iterate
  % method1 sets this equal to the surface divergence of
  % the previous time step
else
  if any(vesicle.viscCont ~= 1)
    if o.profile
      tic
    end
    valTen = vesicle.surfaceDiv(...
        o.solveIminusD(Gf+Fslp+Fwall2Ves+LetsVes,vesicle));
    if o.profile
      fprintf('Solve system alpha*I - DLP          %5.1e\n',toc);
    end
  else
    valTen = -1/o.dt*vesicle.surfaceDiv(valPos);
  end
  % method2, method3, and method4 sets the surface divergence of the sum
  % of single-layer potentials due to bending and tension plus the
  % farField to zero.  The only difference between the two methods is
  % how the preconditioner is used.  method2 uses the possibly
  % ill-conditioned full matrix where as method3 and method 4 use the
  % two schur complements.  Eventually will phase out method2 and it
  % will be fully replaced by method3 and method4
end
% Two possible discretizations of the inextensibility condition
% END OF EVALUATING INEXTENSIBILITY CONDITION

valPos = valPos + o.beta*Xm;
% beta times solution coming from time derivative

val = zeros(3*N*nv,1);
% Initialize output from vesicle and inextensibility equations to zero
for k=1:nv
  val((k-1)*3*N+1:3*k*N) = [valPos(:,k);valTen(:,k)];
end
% Stack val as [x-coordinate;ycoordinate;tension] repeated
% nv times for each vesicle

if (any(vesicle.viscCont ~= 1) && ...
      strcmp(o.vesves,'implicit') && o.confined)
  % This combination of options causes problems with
  % the scaling of the preconditioner.  Need to
  % get rid of the potentially small value o.dt
  valWalls = valWalls * o.dt;
end

val = [val;valWalls(:);valLets];
% Stack velocity along the solid walls in same manner as above
% Stack the stokeslets and rotlet componenets at the end

end % TimeMatVec


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vesVel,divVesVel,wallVel,residual] = ...
      computeVelocities(o,vesicle,Galpert,D,walls,...
          etaProv,RSprov,NearV2V,NearV2W,NearW2V,NearW2W)
% [vesVel,divVesVel,wallVel,residual] = ...
%      computeVelocities(vesicle,Galpert,D,walls,...
%       etaProv,RSprov,NearV2V,NearV2W,NearW2V,NearW2W);
% computes the velocity induced by the provisional vesicle position,
% tension, and density function on the vesicles and the walls.  Also
% returns the vesicle divergence of the velocity field and the residual
% of the picard integral formulation of the vesicle velocity.  These
% quantities are needed to form the modificiations of the right-hand
% side when doing SDC updates

N = vesicle.N;
nv = vesicle.nv;
if o.confined
  Nbd = walls.N;
  nvbd = walls.nv;
else
  Nbd = 0;
  nvbd = 0;
end

vesVel = zeros(2*N,nv,abs(o.orderGL));
% velocity on the vesicles due to the provisional solution 
divVesVel = zeros(N,nv,abs(o.orderGL));
% vesicle divergence of the vesicles due to the provisional
% solution
wallVel = zeros(2*Nbd,nvbd,abs(o.orderGL));
% velocity on the solid walls due to the provisional solution
residual = zeros(2*N,nv,abs(o.orderGL));
% residual of the Picard integral coming from the time
% derivative term of the vesicle position

vesves = o.vesves;
o.vesves = 'implicit';
% need to use implicit so that other vesicles are used to compute
% the integrand z
order = o.order;
o.order = 1;
[o.Xcoeff,o.rhsCoeff,o.beta] = o.getCoeff(o.order);
% need to save the time stepping order
dt = o.dt;
% need to save current time step
o.dt = 1;
% avoid introducing numerical error by multiplying and dividing
% by a potentially small number

if ~o.confined
  z = zeros(3*N*nv,1);
else
  z = zeros(3*N*nv + 2*Nbd*nvbd + 3*(nvbd-1),1);
end
op = o.op;
for n = 1:abs(o.orderGL)
  % loop over the configurations
  o.Galpert = Galpert(:,:,:,n);
  if any(vesicle(n).viscCont ~= 1)
    o.D = D(:,:,:,n);
  end
  % Already computed the single-layer potential when forming the
  % provisional solution Store the single- and double-layer potentials
  % in the object o
  o.NearV2V = NearV2V{n};
  o.NearV2W = NearV2W{n};
  o.NearW2V = NearW2V{n};
  o.NearW2W = NearW2W{n};
  % retrieve near-singular integration strucutre at Gauss-Lobatto time
  % step n

  for k = 1:nv
    z(3*(k-1)*N+1:3*k*N) = [vesicle(n).X(:,k);vesicle(n).sig(:,k)];
  end
  % put positions and tension into z the way that TimeMatVec requires

  for k = 1:nvbd
    istart = 3*nv*N+2*(k-1)*Nbd+1;
    iend = istart + 2*Nbd - 1;
    z(istart:iend) = etaProv(:,k,n);
  end
  % put in density function on solid walls into z the way that
  % TimeMatVec requires them

  for k = 2:nvbd
    istart = 3*nv*N+2*nvbd*Nbd+3*(k-2)+1;
    iend = istart + 2;
    z(istart:iend) = RSprov(:,k,n);
  end
  % put the rotlet and stokeslet coefficients into z

  viscCont = vesicle(n).viscCont;
  vesicle(n).viscCont = ones(1,nv);
  % don't want the double-layer contribution when computing the velocity
  % of the vesicle
  z = o.TimeMatVec(z,vesicle(n),walls);
  % use TimeMatVec to find velocity on vesicles and solid
  % walls due to the provisional solution
  vesicle(n).viscCont = viscCont;
  % set the viscosity contrast back to its original value

  rhs = zeros(2*N,nv);
  for k = 1:nv
    istart = 3*(k-1)*N+1;
    iend = istart + 2*N - 1;
    rhs(:,k) = -z(istart:iend);
  end
  if ~o.confined
    rhs = rhs + o.farField(vesicle(n).X);
  end

  vesVel(:,:,n) = vesicle(n).X + rhs;
  % form the velocity on the vesicle due to the current provisional
  % solution

  if any(vesicle(n).viscCont ~= 1)
    if o.profile
      tic
    end
    vesVel(:,:,n) = o.solveIminusD(vesVel(:,:,n),vesicle(n));
    if o.profile
      fprintf('Solve system alpha*I - DLP          %5.1e\n',toc);
    end
  end 
  % need to apply inv(alpha*I - DLP) if there is a viscosity contrast to
  % obtain the velocity of the vesicles.

  divVesVel(:,:,n) = vesicle(n).surfaceDiv(vesVel(:,:,n));
  % surface divergence of the velocity of the vesicle due
  % to the provisional solution

  for k = 1:nvbd
    istart = 3*nv*N + (k-1)*2*Nbd + 1;
    iend = istart + 2*Nbd - 1;
    wallVel(:,k,n) = z(istart:iend);
  end

  if any(vesicle(n).viscCont ~= 1) && o.confined
    jump = 1/2*(1-vesicle(n).viscCont);
    DLP = @(X) X*diag(jump) + op.exactStokesDLdiag(vesicle(n),o.D,X);
    if ~o.fmmDLP
      kernel = @op.exactStokesDL;
      kernelDirect = @op.exactStokesDL;
    else
      kernel = @op.exactStokesDLnewfmm;
      kernelDirect = @op.exactStokesDL;
    end

    den = vesVel(:,:,n);
    FDLPwall = op.nearSingInt(vesicle(n),den,DLP,...
        o.NearV2W,kernel,kernelDirect,walls,false,false);
    % Evaulate the velocity due to the viscosity contrast
    % on the walls WITH near-singulation integration
    wallVel(:,:,n) = wallVel(:,:,n) + FDLPwall;
  end
  % velocity on the solid walls due to the provisional solution
end

o.dt = dt;
o.order = order;
[o.Xcoeff,o.rhsCoeff,o.beta] = o.getCoeff(o.order);
o.vesves = vesves;
% change back to original time step and vesicle-vesicle interaction

IvesVel = o.lobattoInt(vesVel);
% integrate the vesicles velocity using quadrature rules 
% that are exact for polynomials defined at the 
% Gauss-Lobatto points
for n = 1:abs(o.orderGL)
  residual(:,:,n) = vesicle(1).X - vesicle(n).X + ...
      o.dt/2 * IvesVel(:,:,n);
end
% compute residual by adding the initial vesicle configuartion and
% subtracting the current vesicle configuartion

end % computeVelocities

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zOut = solveIminusD(o,zIn,vesicle)
% zOut = o.solveIminusD(zIn,vesicle) inverts the linear system \alpha*I
% - DLP where alpha=0.5*(1+viscCont) and DLP is the double-layer
% potential.  This is the opeartor that needs to be inverted to find the
% velocity when there is a viscosity contrast involved.  If there is no
% viscosity contrast, alpha = 1 and DLP = 0 so that zOut = zIn


warning off
%tGMRES = tic;
%[zIn,flag,relres,iter] = gmres(@(X) o.IminusD(X,vesicle),zIn(:),...
%    [],1e-2*o.gmresTol,min(2*vesicle.N*vesicle.nv,o.gmresMaxIter));
if 1
[zIn,flag,relres,iter] = gmres(@(X) o.IminusD(X,vesicle),zIn(:),...
    [],1e-2*o.gmresTol,min(2*vesicle.N*vesicle.nv,o.gmresMaxIter),...
    @(z) o.precoIminusD(z,vesicle));
warning on
else
zIn = o.precoIminusD(zIn(:),vesicle);
end
%fprintf('GMRES time is %4.2e\n',toc(tGMRES))
% solve with block-diagonal preconditioned GMRES.  Integral equation is
% of the form identity + compact.  Need a bit more accuracy in this
% gmres solve as it is an inner iteration within an outer GMRES
% iteration.  This hides the fact that the GMRES solver is not linear

zOut = zeros(2*vesicle.N,vesicle.nv);
for k = 1:vesicle.nv
  zOut(:,k) = zIn((k-1)*2*vesicle.N+1:2*k*vesicle.N);
end
% Sort output in appropriate format

end % solveIminusD


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = IminusD(o,X,vesicle)
% val = IminusD(X,vesicle) does a matrix-vector multiply with the matrix
% (alpha*I - D) where alpha = 0.5*(1+viscCont) and D is the double-layer
% potential.  This matrix needs to be inverted when solving for the
% velocity field of a given vesicle configuration and tension

op = o.op;
N = vesicle.N;
nv = vesicle.nv;
alpha = 0.5*(1+vesicle.viscCont);

Xm = zeros(2*N,nv);
for k = 1:nv
  Xm(:,k) = X((k-1)*2*N+1:k*2*N);
end

val = Xm*diag(alpha);
% "jump" term since we are computing alpha * I - DLP

val = val - op.exactStokesDLdiag(vesicle,o.D,Xm);

% if vesicle.N > 32
%   val = val - op.exactStokesDLdiag(vesicle,o.D,Xm);
% else
%   val = val - op.exactStokesDLdiag(vesicle,o.D,Xm);
% end
% self-interaction term

if o.near
  jump = 0.5*(1-vesicle.viscCont);
  DLP = @(X) X*diag(jump) + op.exactStokesDLdiag(vesicle,o.D,X);
end
% compute term for evaluating limiting value along boundary of domain
if ~o.fmmDLP
  kernel = @op.exactStokesDL;
  kernelDirect = @op.exactStokesDL;
else
  kernel = @op.exactStokesDLnewfmm;
  kernelDirect = @op.exactStokesDL;
end
% kernel for evaulating the double-layer potential

%tic
if ~o.near
  Fdlp = kernel(vesicle,Xm,[]);
else
  Fdlp = op.nearSingInt(vesicle,Xm,DLP,...
    o.NearV2V,kernel,kernelDirect,vesicle,true,false);
end
%fprintf('FMM other interactions %4.2e\n',toc)
% potential due to other vesicles

val = val - Fdlp;
% add potential due to self and other vesicles.  Jump is already
% included
val = val(:);

end % IminusD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = sigDenMatVec(o,sigma,vesicle,walls)
% val = sigDenMatVec(sigma,vesicle,walls,NearV2V,NearV2W) does the
% matvec multiply required to find the tension and density function of a
% given vesicle configuration and farfield or solid wall velocity field

N = vesicle.N; % Number of points
nv = vesicle.nv; % Number of vesicles
if o.confined
  Nbd = walls.N; % Number of points per wall
  nvbd = walls.nv; % Number of walls
else
  Nbd = 0;
  nvbd = 0;
end

u = zeros(2*N,nv);
% need to do some manipulations to the TimeMatVec matrix applied only
% to the first component, which is the velocity
for k = 1:nv
  istart = (k-1)*3*N + 1;
  iend = istart + 2*N - 1;
  u(:,k) = sigma(istart:iend);
end

vesves = o.vesves;
o.vesves = 'implicit';
% want interactions to be implicit so that the most accurate tension and
% density functions are found
inextens = o.solver;
o.solver = 'method1';
dt = o.dt;
o.dt = 1;

sigma = o.TimeMatVec(sigma,vesicle,walls);

% Do a matvec but let the incoming postitions be zero since they are
% handled by the initial condition
% Can overwrite sigma as we don't need it from this point onwards
o.vesves = vesves;
o.solver = inextens;
o.dt = dt;
% change back to old vesicle-vesicle and vesicle-boundary interactions

valVel = zeros(2*N,nv);
valTen = zeros(N,nv);
valDen = zeros(2*Nbd,nvbd);
valRS = zeros(3,(nvbd-1));
for k = 1:nv
  valVel(:,k) = sigma((k-1)*3*N+1:(3*k-1)*N);
  valTen(:,k) = sigma((3*k-1)*N+1:3*k*N);
end

% part that corresponds to the velocity and tension
for k = 1:nvbd
  valDen(:,k) = sigma(3*nv*N+(k-1)*2*Nbd+1:3*nv*N+k*2*Nbd);
end

% part that corresponds to the density function
for k = 2:nvbd
  istart = 3*nv*N + 2*nvbd*Nbd + 3*(k-2) + 1;
  iend = istart + 2;  
  valRS(:,k-1) = sigma(istart:iend);
end
% part that corresponds to the Stokeslets and Rotlets

op = o.op;
if ~o.fmm
  kernel = @op.exactStokesSL;
  kernelDirect = @op.exactStokesSL;
else
  kernel = @op.exactStokesSLfmm;
  kernelDirect = @op.exactStokesSL;
end
% kernel for single-layer potential.  Only difference is if the FMM is
% used or not

f = vesicle.tracJump(u,zeros(N,nv));

% bending due to the velocity
if ~o.near
  Fslp = kernel(vesicle,f,[]);
  if o.confined
    [~,FSLPwall] = kernel(vesicle,f,[],walls.X,(1:nv));
  else
    FSLPwall = [];
  end
else
  SLP = @(X) op.exactStokesSLdiag(vesicle,o.Galpert,X);
  Fslp = op.nearSingInt(vesicle,f,SLP,...
      o.NearV2V,kernel,kernelDirect,vesicle,true,false);
      
  if o.confined
    FSLPwall = op.nearSingInt(vesicle,f,SLP,...
      o.NearV2W,kernel,kernelDirect,walls,false,false);
    % Evaluate single-layer potential due to all vesicles on
    % the solid walls WITH near-singular integration
  else
    FSLPwall = [];
  end
end

alpha = (1+vesicle.viscCont)/2; 
valVel = valVel * diag(alpha);
% multiply top row of matrix by alpha

valVel = valVel + op.exactStokesSLdiag(vesicle,o.Galpert,f) + Fslp;
valDen = valDen - FSLPwall;
% subtract off terms that TimeMatVec introduces but we do not have in
% this linear system

val = [valVel;valTen];
val = [val(:);valDen(:);valRS(:)];
% stack the different components coming from the inextensibility, solid
% walls, and rotlets/stokeslets

end % sigDenMatVec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = letsIntegrals(o,otlets,etaM,walls)
% z = letsIntegrals(stokeslet,rotlet,etaM,walls) integrates the density
% function to enforce constraints on stokeslets and rotlets

Nbd = walls.N;
nvbd = walls.nv;
z = zeros(3*(nvbd-1),1);

for k = 2:nvbd
  stokeslet = otlets(3*(k-2)+1:3*(k-2)+2);
  % two stokeslet terms per inner boundary
  rotlet = otlets(3*(k-1));
  % one rotlet term per inner boundary
  ind = 3*(k-2)+1;
  z(ind) = -2*pi*stokeslet(1) + ...
    sum(etaM(1:Nbd,k).*walls.sa(:,k))*2*pi/Nbd;
  % integral of density function dotted with [1;0]
  % is one stokeslet
  z(ind+1) = -2*pi*stokeslet(2) + ...
    sum(etaM(Nbd+1:2*Nbd,k).*walls.sa(:,k))*2*pi/Nbd;
  % integral of density fuction dotted with [0;1]
  % is the other stokeslet
  z(ind+2) = -2*pi*rotlet + sum(...
    ((walls.X(Nbd+1:2*Nbd,k)).*etaM(1:Nbd,k) - ...
    (walls.X(1:Nbd,k)).*etaM(Nbd+1:2*Nbd,k)).*...
    walls.sa(:,k))*2*pi/Nbd;
  % integral of density function dotted with (-y,x)
  % is the rotlet
end % k

end % letsIntegrals

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vel = RSlets(o,X,center,stokeslet,rotlet)
% vel = RSlets(o,X,center,stokeslet,rotlet) evaluates the velocity due
% to the stokeslet and rotlet terms.  Center of the rotlet and
% stokeslet is contained in center

oc = curve;
[x,y] = oc.getXY(X);
% set of points where we are evaluating the velocity
[cx,cy] = oc.getXY(center);
% the center of the rotlet/stokeslet terms

rho2 = (x-cx).^2 + (y-cy).^2;
% distance squared

LogTerm = -0.5*log(rho2)*stokeslet(1);
rorTerm = 1./rho2.*((x-cx).*(x-cx)*stokeslet(1) + ...
    (x-cx).*(y-cy)*stokeslet(2));
RotTerm = (y-cy)./rho2*rotlet;
velx = 1/4/pi*(LogTerm + rorTerm) + RotTerm;
% x component of velocity due to the stokeslet and rotlet

LogTerm = -0.5*log(rho2)*stokeslet(2);
rorTerm = 1./rho2.*((y-cy).*(x-cx)*stokeslet(1) + ...
    (y-cy).*(y-cy)*stokeslet(2));
RotTerm = -(x-cx)./rho2*rotlet;
vely = 1/4/pi*(LogTerm + rorTerm) + RotTerm;
% y component of velocity due to the stokeslet and rotlet

vel = [velx;vely];
% velocity

end % RSlets

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START OF DIFFERENT PRECONDITIONERS INCLUDING BLOCK-DIAGONAL, ONE FOR
% THE SYSTEM THAT SOLVES FOR THE TENSION AND DENSITY GIVEN A POSITION,
% MULTIGRID IDEAS, SCHUR COMPLEMENTS, AND ANALYTIC (BASED ON A
% CIRCLE).  THE BLOCK-DIAGONAL PRECONDITIONER IS THE MOST ROBUST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = preconditionerBD(o,z)
% val = preconditionBD(z) applies the block diagonal preconditioner
% required by preconditioned-GMRES to the vector z

if (strcmp(o.solver,'method1') || strcmp(o.solver,'method2'))
  nv = size(o.bdiagVes.L,3); % number of vesicles
  N = size(o.bdiagVes.L,1)/3; % number of points
elseif strcmp(o.solver,'method3')
  nv = size(o.bdiagVes.schur,3); % number of vesicles
  N = size(o.bdiagVes.schur,1)/2; % number of vesicles
elseif strcmp(o.solver,'method4')
  nv = size(o.bdiagVes.schur,3); % number of vesicles
  N = size(o.bdiagVes.schur,1); % number of vesicles
end
nvbd = size(o.wallDLP,3); % number of solid walls
Nbd = size(o.wallDLP,1)/2; % number of points per solid wall

zves = z(1:3*N*nv);
% extract the position and tension part.  Solid walls is
% handled in the next section of this routine
valVes = zeros(3*N*nv,1);

if (strcmp(o.solver,'method1') || strcmp(o.solver,'method2'))
  for k=1:nv
    valVes((k-1)*3*N+1:3*k*N) = o.bdiagVes.U(:,:,k)\...
      (o.bdiagVes.L(:,:,k)\zves((k-1)*3*N+1:3*k*N));
  end % k
  % precondition with the block diagonal preconditioner for the
  % vesicle position and tension
elseif strcmp(o.solver,'method3')
  for k = 1:nv
    istart = (k-1)*3*N+1;
    iend = istart + 2*N - 1;
    bx = zves(istart:iend); 
    istart = iend + 1;
    iend = istart + N - 1;
    bsig = zves(istart:iend);
    % seperate position and tension terms

    rhs = bx + o.dt*o.bdiagVes.GT(:,:,k)*o.bdiagVes.DGT(:,:,k)*bsig;
    pos = o.bdiagVes.schur(:,:,k)*rhs;
    rhs = bsig + o.bdiagVes.DGB(:,:,k)*pos;
    sig = o.bdiagVes.DGT(:,:,k)*rhs;
    valVes((k-1)*3*N+1:3*k*N) = [pos;sig];
    % use schur decomposition to operate preconditioner
  end % k
elseif strcmp(o.solver,'method4')
  for k = 1:nv
    istart = (k-1)*3*N+1;
    iend = istart + 2*N - 1;
    bx = zves(istart:iend); 
    istart = iend + 1;
    iend = istart + N - 1;
    bsig = zves(istart:iend);
    % seperate position and tension terms

    rhs = bsig + o.bdiagVes.DGB(:,:,k)*o.bdiagVes.IpBen(:,:,k)*bx;
    sig = o.bdiagVes.schur(:,:,k)*rhs;
    rhs = bx + o.dt*o.bdiagVes.GT(:,:,k)*sig;
    pos = o.bdiagVes.IpBen(:,:,k)*rhs;
    valVes((k-1)*3*N+1:3*k*N) = [pos;sig];
    % use schur decomposition to operate preconditioner
  end % k

end % o.solver

zwalls = z(3*N*nv+1:end);
% part of z from solid walls
valWalls = o.bdiagWall * zwalls;
% this matrix is well-coditioned since it is in the form I + DLP

val = [valVes;valWalls];
% stack the two componenets of the preconditioner

end % preconditionerBD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Mat = wallsPrecond(o,walls)
% wallsPrecond(walls) computes the matrix which is the exact inverse of
% the double-layer potential for stokes flow in a bounded domain.  Used
% in the preconditioner for vesicle simulations and capsules.m/computeEta
% which computes eta and RS when there is no vesicle.

Nbd = walls.N;
nvbd = walls.nv;
oc = curve;
[x,y] = oc.getXY(walls.X);
[nory,norx] = oc.getXY(walls.xt);
nory = -nory;
sa = walls.sa;
[cx,cy] = oc.getXY(walls.center);

M11 = zeros(2*Nbd*nvbd,2*Nbd*nvbd);
M12 = zeros(2*Nbd*nvbd,3*(nvbd-1));
M21 = zeros(3*(nvbd-1),2*Nbd*nvbd);
% Allocate space for blocks of matrix that carries the double- layer
% potential, rotlets, and stokeslets to the velocity and the conditions
% in (A4) and (A5) in Rahimian et al.

M11(1:2*Nbd,1:2*Nbd) = M11(1:2*Nbd,1:2*Nbd) + o.wallN0(:,:,1);
jump = - 1/2; 
for k = 1:nvbd
  istart = (k-1)*2*Nbd+1;
  iend = 2*k*Nbd;
  M11(istart:iend,istart:iend) = M11(istart:iend,istart:iend) + ...
      jump*eye(2*Nbd) + o.wallDLP(:,:,k);
end

% Self interaction terms with the jump coming from the double layer
% potential

% temporary space for while we build each off-diagonal component of the
% double-layer potential
for ktar = 1:nvbd % loop over targets
  itar = 2*(ktar-1)*Nbd + 1;
  jtar = 2*ktar*Nbd;
  K = [(1:ktar-1) (ktar+1:nvbd)];
  
  D = zeros(2*Nbd,2*Nbd);
  for ksou = K % loop over all other walls
    isou = 2*(ksou-1)*Nbd + 1;
    jsou = 2*ksou*Nbd;
    
    xtar = x(:,ktar); ytar = y(:,ktar);
    xtar = xtar(:,ones(Nbd,1)); 
    ytar = ytar(:,ones(Nbd,1));
    
    xsou = x(:,ksou); ysou = y(:,ksou);
    xsou = xsou(:,ones(Nbd,1))';
    ysou = ysou(:,ones(Nbd,1))';
    
    norxtmp = norx(:,ksou); norytmp = nory(:,ksou);
    norxtmp = norxtmp(:,ones(Nbd,1))';
    norytmp = norytmp(:,ones(Nbd,1))';
    
    satmp = sa(:,ksou);
    satmp = satmp(:,ones(Nbd,1))';
    
    rho2 = (xtar-xsou).^2 + (ytar-ysou).^2;
    
    coeff = 1/pi*((xtar-xsou).*norxtmp + ...
        (ytar-ysou).*norytmp).*satmp./rho2.^2;
    
    D(1:Nbd,:) = 2*pi/Nbd*[coeff.*(xtar-xsou).^2 ...
        coeff.*(xtar-xsou).*(ytar-ysou)];
    D(Nbd+1:end,:) = 2*pi/Nbd*[coeff.*(ytar-ysou).*(xtar-xsou) ...
        coeff.*(ytar-ysou).^2];
    
    M11(itar:jtar,isou:jsou) = D;
  end %end ktar
end %end ksou

for k = 1:nvbd-1
  icol = 3*(k-1)+1;
  istart = 2*k*Nbd+1;
  iend = istart + Nbd - 1;
  M21(icol,istart:iend) = 2*pi/Nbd*sa(:,k+1)';
  M21(icol+2,istart:iend) = 2*pi/Nbd*sa(:,k+1)'.*y(:,k+1)';
  istart = istart + Nbd;
  iend = iend + Nbd;
  M21(icol+1,istart:iend) = 2*pi/Nbd*sa(:,k+1)';
  M21(icol+2,istart:iend) = -2*pi/Nbd*sa(:,k+1)'.*x(:,k+1)';
end % k
% These compute the integral of the density function around each of the
% inner componenents of the geometry


for k = 1:nvbd - 1
  for ktar = 1:nvbd
    rho2 = (x(:,ktar) - cx(k+1)).^2 + (y(:,ktar) - cy(k+1)).^2;
    istart = (ktar-1)*2*Nbd + 1;
    iend = istart + Nbd - 1;

    icol = 3*(k-1)+1;
    M12(istart:iend,icol) = ...
      M12(istart:iend,icol) + ...
      1/4/pi*(-0.5*log(rho2) + (x(:,ktar)-cx(k+1))./rho2.*...
          (x(:,ktar)-cx(k+1)));
    M12(istart + Nbd:iend + Nbd,icol) = ...
      M12(istart + Nbd:iend + Nbd,icol) + ...
      1/4/pi*((x(:,ktar)-cx(k+1))./rho2.*(y(:,ktar)-cy(k+1)));

    icol = 3*(k-1)+2;
    M12(istart:iend,icol) = ...
      M12(istart:iend,icol) + ...
      1/4/pi*((y(:,ktar)-cy(k+1))./rho2.*(x(:,ktar)-cx(k+1)));
    M12(istart + Nbd:iend + Nbd,icol) = ...
      M12(istart + Nbd:iend + Nbd,icol) + ...
      1/4/pi*(-0.5*log(rho2) + (y(:,ktar)-cy(k+1))./rho2.*...
          (y(:,ktar)-cy(k+1)));

    icol = 3*(k-1)+3;
    M12(istart:iend,icol) = ...
      M12(istart:iend,icol) + ...
      (y(:,ktar)-cy(k+1))./rho2;
    M12(istart + Nbd:iend + Nbd,icol) = ...
      M12(istart + Nbd:iend + Nbd,icol) - ...
      (x(:,ktar)-cx(k+1))./rho2;
  end
end
% This is the evaluation of the velocity field due to the stokeslet
% and rotlet terms

M22 = -2*pi*eye(3*(nvbd-1));
% different combinations of the density functions have to multiply to
% 2*pi multiplied by rotlet or stokeslet terms

Mat = ([M11 M12; M21 M22])\eye(2*nvbd*Nbd + 3*(nvbd-1));
% invert the matrix


end % wallsPrecond


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = precoIminusD(o,z,vesicle)
% z = o.precoIminusD(z,vesicle) is the block-diagonal preconditioner
% for the system alpha*I + DLP.  This did not significantly reduce the
% number of GMRES steps, so it is currently not being called

N = vesicle.N;
nv = vesicle.nv;
alpha = 0.5*(1+vesicle.viscCont);
val = zeros(2*N*nv,1);

for k = 1:nv
  istart = (k-1)*2*N + 1;
  iend = k*2*N;
  z(istart:iend) = (alpha(k)*eye(2*N) - o.D(:,:,k))\z(istart:iend);
end

end % precoIminusD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = preconditionerTen(o,z)
% val = preconditionerTen(z) applies the preconditioner to the tension
% term for when we are solving for just the tension and density
% function given a configuration.  Configuration is eliminated using
% the Schur complelent

nv = size(o.bdiagTen,3); % number of vesicles
N = size(o.bdiagTen,1); 
% number of points per vesicle

nvbd = size(o.wallDLP,3); % number of solid walls
Nbd = size(o.wallDLP,1)/2; 
% number of points per solid wall

zves = z(1:N*nv);
% part of z correpsonding to the vesicles
zwall = z(N*nv+1:end-3*(nvbd-1));
% part of z corresponding to the solid walls
zrot = z(end-3*(nvbd-1)+1:end);
% part of z corresonding to rotlets and stokeslets

valVes = zeros(N*nv,1);
for k=1:nv
  valVes((k-1)*N+1:k*N) = o.bdiagTen(:,:,k)*...
    zves((k-1)*N+1:k*N);
end

valWall = o.bdiagWall*[zwall;zrot];

val = [valVes;valWall];
% stack the preconditioned values due to the tension term and
% the solid wall term

end % preconditionerTen


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sig = invDST(o,vesicle,z)
% sig = invDST(o,vesicle,z) solves the equation DST*sig = z using pGMRES

N = vesicle.N; nv = vesicle.nv;

[z,flag,relres,iter] = gmres(@(X) o.DSTMatVec(vesicle,X),...
    z(:),[],1e-2*o.gmresTol,min(o.gmresMaxIter,N),...
    @(X) o.precoDST(vesicle,X));
% need a smaller tolerance here to hide that fact that this does not
% make a liner preconditioner.  Krylov stuff breaks down in theory, but
% this shouldn't come up until the error of the outer iteration is
% smaller than the requested tolerance

sig = zeros(N,nv);
for k = 1:nv
  sig(:,k) = z((k-1)*N+1:k*N);
end

end % invDST

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sig = DSTMatVec(o,vesicle,sig)
% sig = DSTMatVec(vesicle,sig) applies the operator Div*SLP*Tension to
% sig

N = vesicle.N; nv = vesicle.nv;
sigCols = zeros(N,nv);
for k = 1:nv
  sigCols(:,k) = sig((k-1)*N+1:k*N); 
end

tension = vesicle.tensionTerm(sigCols);
% compute Tension * sig
for k = 1:nv
  tension(:,k) = o.Galpert(:,:,k)*tension(:,k);
end
% apply the single-layer potential
sig = vesicle.surfaceDiv(tension);
% compute the surface divergence
sig = sig(:);

end % DSTMatVec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sig = precoDST(o,vesicle,sig)
% sig = precoDST(vesicle,sig) inverts the operator Div*SLP*Tension where
% the geometry is assumed to be a circle that has the same radius as the
% vesicle

N = vesicle.N; nv = vesicle.nv;
rad = vesicle.length/2/pi;
% radius of circle with same length as the vesicles

imodes = 1./[(1:N/2-1)';(N/2:-1:1)'];
for k = 1:nv
  sig((k-1)*N+1:k*N) = fft(sig((k-1)*N+1:k*N));
  sig((k-1)*N+2:k*N) = -4*rad*sig((k-1)*N+2:k*N).*imodes;
  sig((k-1)*N+1:k*N) = ifft(sig((k-1)*N+1:k*N));
end

sig = real(sig);

end % precoDST

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = invIplusSB(o,vesicle,z)
% x = invIplusSB(vesicle,z) solves the equation (I + dt*SLP*Bending)x =
% z using pGMRES 

N = vesicle.N; nv = vesicle.nv;

rad = vesicle.length/2/pi;
% radius of circle with same length as the vesicles
alpha = -o.dt*vesicle.kappa/8/rad^3;
% re-occuring constant
scaling = 1./(1-2*alpha*[(0:N/2-1) (N/2:-1:1)]'.^3);
% scaling that happens to the modes from modulus greater than or equal
% to 2 in the preconditioner
[z,flag,relres,iter] = gmres(@(X) o.IplusSBMatVec(vesicle,X),...
    z(:),[],1e-2*o.gmresTol,min(o.gmresMaxIter,N),...
    @(X) o.precoIplusSB(vesicle,alpha,scaling,X));

x = zeros(2*N,nv);
for k = 1:nv
  x(:,k) = z(2*(k-1)*N+1:k*2*N);
end

end % invIplusSB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = IplusSBMatVec(o,vesicle,x)
% x = IplusSBMatVec(vesicle,x) applies the operator identity + dt *
% SLP*Bending where \kappa is absored into the term Bending

N = vesicle.N; nv = vesicle.nv;
xCols = zeros(2*N,nv);
for k = 1:nv
  xCols(:,k) = x(2*(k-1)*N+1:k*2*N); 
end

xCols = -vesicle.bendingTerm(xCols);
% compute Bending * x
for k = 1:nv
  xCols(:,k) = o.Galpert(:,:,k)*xCols(:,k);
end
x = x + o.dt*xCols(:);

end % IplusSBMatVec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = precoIplusSB(o,vesicle,alpha,scaling,x)
% val = precoIplusSB(vesicle,alpha,scaling,x) solves the equation (I +
% dt*kappa*SLP*Ben)*val = x analytically for a circle

N = vesicle.N; nv = vesicle.nv;
val = zeros(2*N,nv);

for k = 1:nv
  val(:,k) = x((k-1)*2*N+1:k*2*N);
end
val(1:N,:) = fft(val(1:N,:));
val(N+1:2*N,:) = fft(val(N+1:2*N,:));

for k = 1:nv
  val(3:N-1,k) = val(3:N-1,k).*scaling(3:N-1);
  val(N+3:2*N-1,k) = val(N+3:2*N-1,k).*scaling(3:N-1);
end
% all modes with modulus greater than 1 are diagonal

plusOneMode = [val(2,:) ; val(N+2,:)];
minusOneMode = [val(N,:) ; val(2*N,:)];
% need to save these values since the plus and minus one modes
% communicate ie. isn't diagonal

val(N,:) = 1/(4*alpha-1)/(2*alpha-1) *...
  ((2*alpha^2-4*alpha+1)*val(N,:) - ...
  2*1i*alpha^2*val(2*N,:)) + ...
  alpha/(4*alpha-1) * ...
  (plusOneMode(1,:) + 1i*plusOneMode(2,:));
% -1 mode of the x component

val(2*N,:) = 1/(4*alpha-1)/(2*alpha-1) *...
  (2*1i*alpha^2*val(N,:) + ...
  (2*alpha^2-4*alpha+1)*val(2*N,:)) + ...
  alpha/(4*alpha-1) * ...
  (1i*plusOneMode(1,:) - plusOneMode(2,:));
% -1 mode of the y component

val(2,:) = 1/(4*alpha-1)/(2*alpha-1) *...
  ((2*alpha^2-4*alpha+1)*val(2,:) + ...
  2*1i*alpha^2*val(N+2,:)) + ...
  alpha/(4*alpha-1) * ...
  (minusOneMode(1,:) - 1i*minusOneMode(2,:));
% 1 mode of the x component

val(N+2,:) = 1/(4*alpha-1)/(2*alpha-1) *...
  (-2*1i*alpha^2*val(2,:) + ...
  (2*alpha^2-4*alpha+1)*val(N+2,:)) + ...
  alpha/(4*alpha-1) * ...
  (-1i*minusOneMode(1,:) - minusOneMode(2,:));
% 1 mode of the y component

val(1:N,:) = real(ifft(val(1:N,:)));
val(N+1:2*N,:) = real(ifft(val(N+1:end,:)));
val = val(:);

end % precoIplusSB


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vInf = bgFlow(o,X,varargin);
% vInf = bgFlow(X,varargin) computes the velocity field at the points X.
% Default flow is shear with magnitude 1.  Other flows are relaxation,
% extensional, parabolic, Taylor-Green, cubic, invParabolic.  Some flows
% have an option to input a strength of the flow (ie. shear rate,
% extensional rate, etc).  Flows are given by:
%     cubic:          (y^3,x)
%     relaxation:     (0,0)
%     extensional:    (-x,y)
%     parabolic:      (k(1-y/r)^2,0)
%     invParabolic:   (ky^2,0)
%     rotate:         (y,-x)
%     taylorGreen:    (sin(x)cos(y),-cos(x)sin(y))
%     shear:          (ky,0)
%     choke:          poeusille-like flow at intake and outtake
%     doublechoke:    same as choke
%     couette:        rotating boundary
%     doubleCouette   two rotating boundaries
%     quadCouette     four rotation boundaries
%     doubleFlower    same as doubleCouette
%     figureEight     tangential boundary condition
%     tube            shear boundary condition
% k is a 'flow strength', and r controls the width of the parabola
% in the parabolic flow

N = size(X,1)/2; % number of points per vesicle
nv = size(X,2); % number of vesicles

oc = curve;
[x,y] = oc.getXY(X);
% Separate out x and y coordinates

if any(strcmp(varargin,'cubic'))
  vInf = [y.^3;zeros(N,nv)];

elseif any(strcmp(varargin,'relaxation'))
  vInf = zeros(2*N,nv); 

elseif any(strcmp(varargin,'extensional'))
  vInf = [-x;y];

elseif any(strcmp(varargin,'parabolic'))
  R = find(strcmp(varargin,'R'));
  if isempty(R)
%    R = 100;
      R = 10;
  else
    R = varargin{R+1};
  end
  UM = find(strcmp(varargin,'Umax'));
  if isempty(UM)
%    UM = 1e5; % default value for strength of flow
    UM = 2;
  else
    UM = varargin{UM+1};
  end
  vInf = [UM*(1-(y/R).^2);zeros(N,nv)];

elseif any(strcmp(varargin,'invParabolic'))
  k = find(strcmp(varargin,'k'));
  if isempty(k)
    k = 1; % default value for strength of flow
  else
    k = varargin{k+1};
  end
  vInf = [k*y.^2;zeros(N,nv)];

elseif any(strcmp(varargin,'rotate'))
  vInf = [y;-x];
%  ./[x.^2+y.^2;x.^2+y.^2];

elseif any(strcmp(varargin,'taylorGreen'))
  vInf = [sin(x).*cos(y);-cos(x).*sin(y)];

elseif any(strcmp(varargin,'shear'))
  k = find(strcmp(varargin,'k'));
  if isempty(k)
    k = 1; % default value for strength of flow
  else
    k = varargin{k+1};
  end
  vInf = [k*y;zeros(N,nv)];

elseif (any(strcmp(varargin,'choke')) || ...
      any(strcmp(varargin,'doublechoke')) || ...
      any(strcmp(varargin,'choke2')))
  vInf = zeros(2*N,nv);
  ind = abs(x)>7;
  vx = exp(1./((y(ind)/max(y)).^2-1))/exp(-1);
  % typical mollifer so that velocity decays smoothly to 0
  vx(vx==Inf) = 0;
  vInf(ind,:) = vx;
  
elseif any(strcmp(varargin,'couette'));
  vInf = [zeros(2*N,1) 1*[-y(:,2)+mean(y(:,2));x(:,2)-mean(x(:,2))]];
  
elseif any(strcmp(varargin,'couetteOuter'));
  vInf = [1*[-y(:,1)+mean(y(:,1));x(:,1)-mean(x(:,1))] zeros(2*N,1)];  

elseif any(strcmp(varargin,'couette10'));
  vInf = [zeros(2*N,1) 10*[-y(:,2)+mean(y(:,2));x(:,2)-mean(x(:,2))]];
  
elseif any(strcmp(varargin,'couette100'));
  vInf = [zeros(2*N,1) 100*[-y(:,2)+mean(y(:,2));x(:,2)-mean(x(:,2))]];
  
elseif (any(strcmp(varargin,'doubleCouette')) || ...
      any(strcmp(varargin,'doubleFlower')));
  vInf = [zeros(2*N,1) 1*[-y(:,2)+mean(y(:,2));x(:,2)-mean(x(:,2))] ...
      -[y(:,3)-mean(y(:,3));-x(:,3)+mean(x(:,3))]];

elseif (any(strcmp(varargin,'quadCouette')));
  vInf = [zeros(2*N,1) ...
      -[y(:,2)-mean(y(:,2));-x(:,2)+mean(x(:,2))] ...
      -[y(:,3)-mean(y(:,3));-x(:,3)+mean(x(:,3))] ...
      +[y(:,4)-mean(y(:,4));-x(:,4)+mean(x(:,4))] ...
      +[y(:,5)-mean(y(:,5));-x(:,5)+mean(x(:,5))]];

elseif any(strcmp(varargin,'cylinder'))
%  theta = (0:N-1)'*2*pi/N;
%  vInf = [cos(10*theta);sin(2*theta)];
  vInf = 1*[-y+mean(y);x-mean(x)];

elseif any(strcmp(varargin,'figureEight'))
  oc = curve;
  [~,vInf,~] = oc.diffProp([x;y]);
%  vInf(1:end/2,:) = 3*vInf(1:end/2,:);

  sup = find(abs(x)<=1 & y>0);
  sdown = find(abs(x)<=1 & y<0);
  omega = linspace(-1,1,numel(sup)+2)';
  omega = omega(2:end-1);
  mollifier = 4*exp(1./(omega.^2 - 1))+1;
  vInf(sup,:) = vInf(sup,:) .* mollifier;
  vInf(sdown,:) = vInf(sdown,:) .* mollifier;
  % increase the velocity in a smooth fashion near the middle
  % of the solid walls

elseif any(strcmp(varargin,'shear'))
  vInf = [y;zeros(N,nv)];

elseif any(strcmp(varargin,'diffuser'));
  vInf = zeros(2*N,nv);
  ind = abs(x(:,1))>9;
  vx = exp(1./((y(ind,1)/max(y(ind,1))).^2-1))/exp(-1);
  % typical mollifer so that velocity decays smoothly to 0
  vx(vx==Inf) = 0;
  vInf(ind,1) = vx;

elseif any(strcmp(varargin,'microfluidic'));
  oc = curve;
  [~,tangent,~] = oc.diffProp(X); 
  vInf = tangent;
  vInf(:,1) = 0*vInf(:,1);
  vInf(:,2) = +1*vInf(:,2);
  vInf(:,3) = -1*vInf(:,3);
  vInf(:,4) = -1*vInf(:,4);
  vInf(:,5) = +1*vInf(:,5);
  
elseif any(strcmp(varargin,'DlateralD'))
    
  vInf = zeros(2*N,nv);
  vOut = zeros(2*N,1);
  ind = abs(x(:,1))>max(abs(x(:,1)))-max(abs(x(:,1)))*0.1;
  vx = exp(1./((y(ind,1)/max(y(:,1))).^2-1))/exp(-1);
  % typical mollifer so that velocity decays smoothly to 0
  vx(vx==Inf) = 0;
  vOut(ind,1) = vx;
  vInf(:,1) = vOut;
else 
  vInf = [y;zeros(N,nv)];
  % default flow is shear
end

speed = varargin{2};
% speed of the background velocity
vInf = vInf * speed;

end % bgFlow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GLpts = gaussLobatto(o,orderGL)
% GLpts = gaussLobatto(orderGL) loads the Guass- Lobatto points for the
% desired number of points 

if (orderGL == 2 || orderGL == -2)
  GLpts = [-1 1];
elseif (orderGL == 3 || orderGL == -3)
  GLpts = [-1 0 1];
elseif orderGL == 4
  GLpts = [-1 -0.447213595499958 0.447213595499958 1];
elseif orderGL == -4
  GLpts = [-1 -0.333333333333333 0.333333333333333 1];
elseif orderGL == 5
  GLpts = [-1 -0.654653670707977 0 0.654653670707977 1];
elseif orderGL == -5
  GLpts = [-1 -0.5 0 0.5 1];
elseif orderGL == 6
  GLpts = [-1 -0.765055323929465 -0.285231516480645 ...
      0.285231516480645 0.765055323929465 1];
elseif orderGL == -6
  GLpts = [-1 -0.6 -0.2 0.2 0.6 1];
elseif orderGL == 7
  GLpts = [-1 -0.830223896278567 -0.468848793470714 ...
    0 0.468848793470714 0.830223896278567 1];
elseif orderGL == -7
  GLpts = [-1 -0.6666666666666666  -0.3333333333333333 ...
    0 0.3333333333333333 0.6666666666666667 1];
else
  fprintf('**************************************\n')
  fprintf('NO GAUSS-LOBATTO POINTS FOR THIS ORDER\n')
  fprintf('**************************************\n')
  pause
end

end % gaussLobatto

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fint = lobattoInt(o,f)
% fint = lobattoInt(f) returns the integral of f from [-1,t] where t
% ranges from -1 to 1.  If orderGL > 0, the interior points are
% Gauss-Lobato points.  Otherwise, they are equispaced points.  The
% quadrature rules were found by first interpolating f with a
% polynomial and then integrating the polynomial exactly.  All
% coefficients are precomputed f is of size (N,nv,order) where N is the
% number of points that need to be integrated, nv is a second index
% which in this case corresponds to number of vesicles, and order is
% the number of time steps

nv = size(f,2);
% number of vesicles
fint = zeros(size(f));
% initialize the integral to zero

if (o.orderGL == 2 || o.orderGL == -2)
  t = [-1 1];

  for n = 1:2
    for k = 1:nv
      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^2*0.25*(-f(:,k,1)+f(:,k,2));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n)*0.5*(f(:,k,1)+f(:,k,2));

      fint(:,k,n) = fint(:,k,n) + ...
          (0.75*f(:,k,1) + 0.25*f(:,k,2));
    end
  end
  % order 1 Gauss-Lobatto or equispaced


elseif (o.orderGL == 3 || o.orderGL == -3)
  t = [-1 0 1];

  for n = 1:3
    for k = 1:nv
      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^3*(0.1666666666666667*(f(:,k,1)+f(:,k,3)) - ...
          0.3333333333333333*f(:,k,2));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^2*(-0.25*(f(:,k,1)-f(:,k,3)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n)*f(:,k,2);

      fint(:,k,n) = fint(:,k,n) + ...
          0.4166666666666667*f(:,k,1) + ...
          0.6666666666666667*f(:,k,2) - ...
          0.0833333333333333*f(:,k,3);
    end
  end
  % order 2 Gauss-Lobatto or equispaced


elseif (o.orderGL == 4)
  t = [-1 -0.447213595499958 0.447213595499958 1];

  for n = 1:4
    for k = 1:nv
      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^4*(-0.15625*(f(:,k,1)-f(:,k,4)) + ...
          0.3493856214843422*(f(:,k,2)-f(:,k,3)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^3*(0.2083333333333333*(f(:,k,1)+f(:,k,4)) - ...
          0.2083333333333333*(f(:,k,2)+f(:,k,3)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^2*(0.0625*(f(:,k,1)-f(:,k,4)) - ...
          0.6987712429686845*(f(:,k,2)-f(:,k,3)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n)*(-0.125*(f(:,k,1)+f(:,k,4)) + ...
          0.625*(f(:,k,2)+f(:,k,3)));

      fint(:,k,n) = fint(:,k,n) + ...
          (0.1770833333333333*f(:,k,1) + ...
          0.7660522881510090*f(:,k,2) + ...
          0.0672810451823247*f(:,k,3) - ...
          0.0104166666666666*f(:,k,4));
    end
  end
  % order 3 Gauss-Lobatto

elseif (o.orderGL == -4)
  t = [-1 -0.333333333333333 0.333333333333333 1];

  for n = 1:4
    for k = 1:nv
      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^4*(-0.140625*(f(:,k,1)-f(:,k,4)) + ...
          0.421875*(f(:,k,2)-f(:,k,3)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^3*(0.1875*(f(:,k,1)+f(:,k,4)) - ...
          0.1875*(f(:,k,2)+f(:,k,3)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^2*(0.03125*(f(:,k,1)-f(:,k,4)) - ...
          0.84375*(f(:,k,2)-f(:,k,3)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n)*(-0.0625*(f(:,k,1)+f(:,k,4)) + ...
          0.5625*(f(:,k,2)+f(:,k,3)));

      fint(:,k,n) = fint(:,k,n) + ...
          (0.234375*f(:,k,1) + ...
          0.796875*f(:,k,2) - ...
          0.046875*f(:,k,3) + ...
          0.015625*f(:,k,4));
    end
  end
  % order 3 equispaced

elseif (o.orderGL == 5)
  t = [-1 -0.654653670707977 0 0.654653670707977 1];

  for n = 1:5
    for k = 1:nv
      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^5*(0.175*(f(:,k,1)+f(:,k,5)) - ...
          0.4083333333333333*(f(:,k,2)+f(:,k,4)) + ...
          0.4666666666666666*f(:,k,3));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^4*(-0.21875*(f(:,k,1)-f(:,k,5)) + ...
          0.3341461444238633*(f(:,k,2)-f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^3*(-0.125*(f(:,k,1)+f(:,k,5)) + ...
          0.6805555555555555*(f(:,k,2)+f(:,k,4)) - ...
          1.1111111111111111*f(:,k,3));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^2*(0.1875*(f(:,k,1)-f(:,k,5)) - ...
          0.6682922888477265*(f(:,k,2)-f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n)*f(:,k,3);

      fint(:,k,n) = fint(:,k,n) + ...
          0.08125*f(:,k,1) + ...
          0.6063683666460855*f(:,k,2) + ...
          0.3555555555555555*f(:,k,3) - ...
          0.0619239222016411*f(:,k,4) + ...
          0.01875*f(:,k,5);
    end
  end
  % order 4 Gauss-Lobatto

elseif (o.orderGL == -5)
  t = [-1 -0.5 0 0.5 1];

  for n = 1:5
    for k = 1:nv
      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^5*(0.1333333333333333*(f(:,k,1)+f(:,k,5)) - ...
          0.53333333333333333*(f(:,k,2)+f(:,k,4)) + ...
          0.8*f(:,k,3));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^4*(-0.16666666666666667*(f(:,k,1)-f(:,k,5)) + ...
          0.3333333333333333*(f(:,k,2)-f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^3*(-0.0555555555555556*(f(:,k,1)+f(:,k,5)) + ...
          0.8888888888888889*(f(:,k,2)+f(:,k,4)) - ...
          1.6666666666666667*f(:,k,3));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^2*(0.0833333333333333*(f(:,k,1)-f(:,k,5)) - ...
          0.6666666666666667*(f(:,k,2)-f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n)*f(:,k,3);

      fint(:,k,n) = fint(:,k,n) + ...
          0.1611111111111111*f(:,k,1) + ...
          0.6888888888888889*f(:,k,2) + ...
          0.1333333333333333*f(:,k,3) + ...
          0.0222222222222222*f(:,k,4) - ...
          0.0055555555555556*f(:,k,5);
    end
  end
  % order 4 equi-spaced points

elseif (o.orderGL == 6)
  t = [-1 -0.765055323929465 -0.285231516480645 ...
      0.285231516480645 0.765055323929465 1];

  for n = 1:6
    for k = 1:nv
      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^6*(-0.21875*(f(:,k,1)-f(:,k,6)) + ...
          0.5212094304495727*(f(:,k,2)-f(:,k,5)) - ...
          0.6310805056491861*(f(:,k,3)-f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^5*(0.2625*(f(:,k,1)+f(:,k,6)) - ...
          0.4785048595772274*(f(:,k,2)+f(:,k,5)) + ...
          0.2160048595772274*(f(:,k,3)+f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^4*(0.21875*(f(:,k,1)-f(:,k,6)) - ...
          0.8454202131918329*(f(:,k,2)-f(:,k,5)) + ...
          1.500687022042463*(f(:,k,3)-f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^3*(-0.2916666666666666*(f(:,k,1)+f(:,k,6)) + ...
          0.8623909800799931*(f(:,k,2)+f(:,k,5)) - ...
          0.5707243134133265*(f(:,k,3)+f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^2*(-0.03125*(f(:,k,1)-f(:,k,6)) + ...
          0.1272121350349483*(f(:,k,2)-f(:,k,5)) - ...
          1.108132527137368*(f(:,k,3)-f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n)*(0.0625*(f(:,k,1)+f(:,k,6)) - ...
          0.1946486423538423*(f(:,k,2)+f(:,k,5)) + ...
          0.6321486423538425*(f(:,k,3)+f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          0.06458333333333333*f(:,k,1) + ...
          0.3862361258562357*f(:,k,2) + ...
          0.5159551992618339*f(:,k,3) + ...
          0.03890317777365201*f(:,k,4) - ...
          0.007761169558388661*f(:,k,5) + ...
          0.002083333333333333*f(:,k,6);
    end
  end
  % order 5 Gauss-Lobatto

elseif (o.orderGL == -6)
  t = [-1 -0.6 -0.2 0.2 0.6 1];

  for n = 1:6
    for k = 1:nv
      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^6*(-0.1356336805555556*(f(:,k,1)-f(:,k,6)) + ...
          0.6781684027777778*(f(:,k,2)-f(:,k,5)) - ...
          1.3563368055555556*(f(:,k,3)-f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^5*(0.1627604166666667*(f(:,k,1)+f(:,k,6)) - ...
          0.48828125*(f(:,k,2)+f(:,k,5)) + ...
          0.3255208333333333*(f(:,k,3)+f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^4*(0.0813802083333333*(f(:,k,1)-f(:,k,6)) - ...
          1.0579427083333333*(f(:,k,2)-f(:,k,5)) + ...
          2.7669270833333333*(f(:,k,3)-f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^3*(-0.108506944444444*(f(:,k,1)+f(:,k,6)) + ...
          0.8463541666666667*(f(:,k,2)+f(:,k,5)) - ...
          0.7378472222222222*(f(:,k,3)+f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^2*(-0.005859375*(f(:,k,1)-f(:,k,6)) + ...
          0.0813802083333333*(f(:,k,2)-f(:,k,5)) - ...
          1.46484375*(f(:,k,3)-f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n)*(0.01171875*(f(:,k,1)+f(:,k,6)) - ...
          0.09765625*(f(:,k,2)+f(:,k,5)) + ...
          0.5859375*(f(:,k,3)+f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          0.1260850694444444*f(:,k,1) + ...
          0.5588107638888889*f(:,k,2) + ...
          0.2278645833333333*f(:,k,3) + ...
          0.1193576388888889*f(:,k,4) - ...
          0.0379774305555556*f(:,k,5) + ...
          0.005859375*f(:,k,6);
    end
  end
  % order 5 equi-spaced points

elseif (o.orderGL == 7)
  t = [-1.000000000000000 -0.830223896278567 ...
      -0.468848793470714 0 0.468848793470714 ...
      0.830223896278567 1.000000000000000];

  for n = 1:7
    for k = 1:nv
      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^7*(0.29464285714285718*(f(:,k,1)+f(:,k,7)) - ...
          0.71040995801242226*(f(:,k,2)+f(:,k,6)) + ...
          0.88719567229813686*(f(:,k,3)+f(:,k,5)) - ...
          0.94285714285714356*f(:,k,4));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^6*(-0.34375*(f(:,k,1)-f(:,k,7)) + ...
          0.68809921051219412*(f(:,k,2)-f(:,k,6)) - ...
          0.48528739061765716*(f(:,k,3)-f(:,k,5)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^5*(-0.375*(f(:,k,1)+f(:,k,7)) + ...
          1.21320038050366995*(f(:,k,2)+f(:,k,6)) - ...
          2.09820038050367082*(f(:,k,3)+f(:,k,5)) + ...
          2.52000000000000176*f(:,k,4));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^4*(0.46875*(f(:,k,1)-f(:,k,7)) - ...
          1.25903493358549612*(f(:,k,2)-f(:,k,6)) + ...
          1.22967339607367386*(f(:,k,3)-f(:,k,5)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^3*(0.10416666666666660*(f(:,k,1)+f(:,k,7)) - ...
          0.36437739881046465*(f(:,k,2)+f(:,k,6)) + ...
          1.42687739881046537*(f(:,k,3)+f(:,k,5)) - ...
          2.3333333333333346*f(:,k,4));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^2*(-0.15625*(f(:,k,1)-f(:,k,7)) + ...
          0.45377223563440987*(f(:,k,2)-f(:,k,6)) - ...
          1.00348462029437623*(f(:,k,3)-f(:,k,5)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n)*(1.0*f(:,k,4));

      fint(:,k,n) = fint(:,k,n) + ...
          0.55059523809523700e-1*f(:,k,1) + ...
          0.25557651111967517*f(:,k,2) + ...
          0.47497130544329094*f(:,k,3) + ...
          0.24380952380952357*f(:,k,4) - ...
          0.4322592423342813e-1*f(:,k,5) + ...
          0.2124953624189091e-1*f(:,k,6) - ...
          0.7440476190476161e-2*f(:,k,7);
    end
  end
  % order 6 Gauss-Lobatto

elseif (o.orderGL == -7)
  t = [-1 -0.6666666666666666  ...
      -0.3333333333333333 0 0.3333333333333333 ...
      0.6666666666666667 1];

  for n = 1:7
    for k = 1:nv
      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^7*(0.1446428571428571*(f(:,k,1)+f(:,k,7)) - ...
          0.8678571428571429*(f(:,k,2)+f(:,k,6)) + ...
          2.1696428571428571*(f(:,k,3)+f(:,k,5)) - ...
          2.8928571428571429*f(:,k,4));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^6*(-0.16875*(f(:,k,1)-f(:,k,7)) + ...
          0.675*(f(:,k,2)-f(:,k,6)) - ...
          0.84375*(f(:,k,3)-f(:,k,5)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^5*(-0.1125*(f(:,k,1)+f(:,k,7)) + ...
          1.35*(f(:,k,2)+f(:,k,6)) - ...
          4.3875*(f(:,k,3)+f(:,k,5)) + ...
          6.3*f(:,k,4));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^4*(0.140625*(f(:,k,1)-f(:,k,7)) - ...
          1.125*(f(:,k,2)-f(:,k,6)) + ...
          1.828125*(f(:,k,3)-f(:,k,5)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^3*(0.0166666666666667*(f(:,k,1)+f(:,k,7)) - ...
          0.225*(f(:,k,2)+f(:,k,6)) + ...
          2.25*(f(:,k,3)+f(:,k,5)) - ...
          4.0833333333333333*f(:,k,4));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^2*(-0.025*(f(:,k,1)-f(:,k,7)) + ...
          0.225*(f(:,k,2)-f(:,k,6)) - ...
          1.125*(f(:,k,3)-f(:,k,5)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n)*(1.0*f(:,k,4));

      fint(:,k,n) = fint(:,k,n) + ...
          0.1019345238095238*f(:,k,1) + ...
          0.4821428571428571*f(:,k,2) + ...
          0.1727678571428571*f(:,k,3) + ...
          0.3238095238095238*f(:,k,4) - ...
          0.1084821428571429*f(:,k,5) + ...
          0.03214285714285714*f(:,k,6) - ...
          0.004315476190476190*f(:,k,7);
    end
  end
  % order 6 equi-spaced points
end % o.orderGL

end % lobattoInt

end % methods

end % classdef
