clear all; clc

fprintf('Simple elliptical vesicle in a shear flow.\n');
fprintf('First-order semi-implicit time stepping.\n');

% Physics parameters
prams.N = 64;               % points per vesicle
prams.nv = 1;               % number of vesicles
prams.T = 100;               % time horizon (two tumbling)
prams.m = 2000;              % number of time steps
prams.kappa = 1e-1;         % bending coefficient
prams.viscCont = 1;         % viscosity contrast
options.farField = 'shear'; % background velocity
options.order = 1;          % time stepping order
options.vesves = 'implicit';
% Discretization of vesicle-vesicle interactions.
% Either 'explicit' or 'implicit'
options.inextens = 'method1';
options.near = true;        % near-singular integration
options.fmm = false;
options.antiAlias = false;
prams.gmresMaxIter = 3*prams.N;
prams.gmresTol = 1e-10;
prams.errorTol = 1;

% ADD-ONS
options.alignCenterAngle = true;
options.correctShape = true;
options.reparameterization = false;


% TIME ADAPTIVITY (parameters for new implementation)
options.timeAdap = true;

prams.rtolArea = 5e-3;
prams.rtolLength = 5e-3;
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
options.nsdc = 0;
options.expectedOrder = 1;

% Plot on-the-fly
options.usePlot = true;
options.track = false;
% Save vesicle information and create a log file
options.logFile = 'output/shear1Ves.log';
% Name of log file for saving messages
options.dataFile = 'output/shear1VesData.bin';
% Name of binary data file for storing vesicle information

options.saveError = true;
options.errorFile = 'output/shear1VesError.bin';
% Name of binary data file for storing truncation errors after each step

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

oc = curve;
X = oc.initConfig(prams.N,'reducedArea',0.65,'angle',pi/2);
% Initial configuration of reduce area 0.65 and aligned

Ves2D(X,[],prams,options);
% Run vesicle code

