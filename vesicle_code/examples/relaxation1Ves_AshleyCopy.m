clear all; clc

fprintf('Simple elliptical vesicle in a relaxation flow.\n');
fprintf('First-order semi-implicit time stepping.\n');

% Physics parameters
prams.N = 128;               % points per vesicle
prams.nv = 1;               % number of vesicles
prams.T = 20;               % time horizon (two tumbling)
prams.m = 800;              % number of time steps
prams.kappa = 1e-1;         % bending coefficient
prams.viscCont = 1;         % viscosity contrast
options.farField = 'relaxation'; % background velocity
options.farFieldSpeed = 1;
options.order = 1;          % time stepping order
options.vesves = 'implicit';
% Discretization of vesicle-vesicle interactions.
% Either 'explicit' or 'implicit'
options.inextens = 'method1';
options.near = true;        % near-singular integration
options.fmm = false;
options.antiAlias = false;
options.semipermeable = true;
prams.gmresMaxIter = 3*prams.N;
prams.gmresTol = 1e-10;
prams.errorTol = 1000;
prams.PhysBeta = 1;

% ADD-ONS
options.alignCenterAngle = false;
options.correctShape = false;
options.reparameterization = false;

% TIME ADAPTIVITY (parameters for new implementation)
options.timeAdap = true;

prams.rtolArea = 1e+10;
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
options.nsdc = 0;
options.expectedOrder = 1;

options.expForce = false;

% Plot on-the-fly
options.usePlot = true;
options.axis = [-5 5 -5 5];
options.track = false;
% Save vesicle information and create a log file
options.logFile = 'output/1Vesicle_Shear2_beta1_RApt65.log';
% Name of log file for saving messages
options.dataFile = 'output/1Vesicle_Shear2_beta1_RApt65.bin';
% Name of binary data file for storing vesicle information

options.saveError = true;
options.errorFile = 'output/relaxation1VesError.bin';
% Name of binary data file for storing truncation errors after each step

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

oc = curve;
ra = 0.65;
%scale = sqrt(ra);
scale = 1;
X = oc.initConfig(prams.N,...
    'reducedArea',ra,...
    'angle',pi/2,...
    'center',[0;0],...
    'scale',scale);
% Initial configuration of reduce area 0.65 and aligned
%theta = (0:prams.N-1)'*2*pi/prams.N;
%X = [0.99*cos(theta);sin(theta)];



Xfinal = Ves2D(X,[],prams,options);


% Run vesicle code
