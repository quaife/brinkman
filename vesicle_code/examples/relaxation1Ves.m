%function[Xfinal] = relaxation1Ves(beta, m, betas, ms)
fprintf('Simple elliptical vesicle in a relaxation flow.\n');
fprintf('First-order semi-implicit time stepping.\n');

beta = 0;
m = 100;
ms = m;
betas = beta;

% Physics parameters
prams.N = 128;               % points per vesicle
prams.nv = 1;               % number of vesicles
prams.T = 10;               % time horizon (two tumbling)
prams.m = m;              % number of time steps
prams.kappa = ones(prams.nv,1); % bending coefficient
prams.viscCont = ones(prams.nv,1);         % viscosity contrast
options.farField = 'relaxation'; % background velocity
options.order = 1;          % time stepping order
options.vesves = 'implicit';
% Discretization of vesicle-vesicle interactions.
% Either 'explicit' or 'implicit'
options.inextens = 'method1';
options.near = true;        % near-singular integration
options.fmm = false;
options.antiAlias = false;
options.semipermeable = false;
options.adhesion = false;
defaultPram.adStrength = 1;
defaultPram.adRange = 4e-1;
prams.gmresMaxIter = 3*prams.N;
prams.gmresTol = 1e-10;
prams.errorTol = 1000;
if ~options.semipermeable
  beta = 0;
end
prams.fluxCoeff = beta*ones(prams.nv,1);

% ADD-ONS
options.alignCenterAngle = false;
options.correctShape = false;
options.reparameterization = false;

% TIME ADAPTIVITY (parameters for new implementation)
options.timeAdap = false;
prams.rtolArea = 1e10;
prams.rtolLength = 1e-2;
if 1
  prams.dtMax = 2;
  prams.dtMin = 1e-4;
  prams.betaInc = 1e-1;
  prams.betaDec = 5e-1;
  prams.betaUp = 1.2;
  prams.betaDown = 0.5;
  prams.alpha = 0.9;
end

options.orderGL = 2;
options.nsdc = 1;
options.expectedOrder = 2;

options.expForce = false;

% Plot on-the-fly
options.usePlot = true;
options.axis = [-3 3 -4 4];
options.track = false;
% Save vesicle information and create a log file
%options.logFile = sprintf('%s','output/nsdc1_relaxation2Ves',ms,'_',betas,'.log');
options.logFile = 'output/relaxation1Ves.log';
% Name of log file for saving messages
%options.dataFile = sprintf('%s','output/nsdc1_relaxation2Ves',ms,'_',betas,'.bin');
options.dataFile = 'output/relaxation1VesData.bin';
% Name of binary data file for storing vesicle information

options.saveError = true;
options.errorFile =sprintf('%s','output/nsdc1_relaxation2Ves',ms,'_',betas,'Error.log');
% Name of binary data file for storing truncation errors after each step

prams.fluxShape = beta*ones(prams.N,1); %flux shape 1

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

oc = curve;

ra = 0.85;
scale = 1;

%centerx = [-5 5];
%centery = [0 0];
%ang = [pi/2 pi/2];
centerx = 0;
centery = 0;
ang = pi/2;

X = oc.initConfig(prams.N,'nv',prams.nv,...
    'reducedArea',ra,...
    'center',[centerx;centery],...
    'angle',ang,...
    'scale',scale);

Xfinal = Ves2D(X,[],prams,options);
% Run vesicle code
%end
