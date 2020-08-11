%function[Xfinal] = relaxation1Ves(beta, m, betas, ms)
fprintf('Simple elliptical vesicle in a relaxation flow.\n');
fprintf('First-order semi-implicit time stepping.\n');

beta = 0;
m = 100000;
ms = m;
betas = beta;

% Physics parameters
prams.N = 96;               % points per vesicle
prams.nv = 1;               % number of vesicles
prams.T = 10;               % time horizon (two tumbling)
prams.m = 100*10;              % number of time steps
prams.kappa = ones(prams.nv,1); % bending coefficient
prams.viscCont = ones(prams.nv,1);         % viscosity contrast
options.farField = 'relaxation'; % background velocity
options.near = true;        % near-singular integration
options.fmm = false;
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

% TIME ADAPTIVITY (parameters for new implementation)
options.timeAdap = false;
prams.rtolArea = 1e10;
prams.rtolLength = 1e-2;
prams.betaUp = 1.2;
prams.betaDown = 0.5;
prams.alpha = 0.9;
prams.dtMax = 1e2;
prams.dtMin = 1e-5;

options.orderGL = 2;
options.nsdc = 0;
options.expectedOrder = 1;

options.expForce = false;

% Plot on-the-fly
options.usePlot = true;
options.axis = [-3 3 -4 4];
% Save vesicle information and create a log file
%options.logFile = sprintf('%s','output/nsdc1_relaxation2Ves',ms,'_',betas,'.log');
options.logFile = 'output/relaxation1Ves.log';
% Name of log file for saving messages
%options.dataFile = sprintf('%s','output/nsdc1_relaxation2Ves',ms,'_',betas,'.bin');
options.dataFile = 'output/relaxation1VesData.bin';
% Name of binary data file for storing vesicle information

options.fluxShape = 1; % constant value
%options.fluxShape = 2; % gating via tension

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

oc = curve;

ra = 0.65;
%ra = 0.99;
scale = sqrt(ra)/2;

centerx = 0;
centery = 0;
ang = pi/2;

X = oc.initConfig(prams.N,'nv',prams.nv,...
    'reducedArea',ra,...
    'center',[centerx;centery],...
    'angle',ang,...
    'scale',scale);

th = (0:prams.N-1)'*2*pi/prams.N;
X = [cos(th);3*sin(th)];

%X = oc.initConfig(prams.N,'star',...
%    'folds',5,...
%    'nv',prams.nv,...
%    'reducedArea',ra,...
%    'center',[centerx;centery],...
%    'angle',ang,...
%    'scale',1.18);


Xfinal = Ves2D(X,[],prams,options);
% Run vesicle code
%end
