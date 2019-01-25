clear all; clc

fprintf('Simple elliptical vesicles in a relaxation flow.\n');

% Physics parameters
prams.N = 256;               % points per vesicle
prams.nv = 6;               % number of vesicles
prams.T = 100;               % time horizon (two tumbling)
prams.m = 2000;             % number of time steps
prams.kappa = ones(1,prams.nv);         % bending coefficient
prams.viscCont = ones(1,prams.nv);         % viscosity contrast
options.farField = 'relaxation'; % background velocity
options.farFieldSpeed = 1;
options.order = 1;          % time stepping order
options.vesves = 'implicit';
% Discretization of vesicle-vesicle interactions.
% Either 'explicit' or 'implicit'
options.inextens = 'method1';
options.near = true;        % near-singular integration
options.fmm = false;
options.verbose = true;
options.antiAlias = false;
prams.gmresMaxIter = 3*prams.N;
prams.gmresTol = 1e-10;
prams.errorTol = 1;

% ADD-ONS
options.correctShape = false;
options.adhesion = true;
prams.adRange = 4e-1;
prams.adStrength = 7e-1;

% TIME ADAPTIVITY (parameters for new implementation)
options.timeAdap = true;

prams.rtolArea = 1e-2;
prams.rtolLength = 1e-2;
prams.betaUp = 1.2;
prams.betaDown = 0.5;
prams.alpha = 0.9;

options.orderGL = 2;
options.nsdc = 1;
options.expectedOrder = 1;

% Plot on-the-fly
options.usePlot = true;
options.axis = [-5 5 -5 5];
options.track = false;
% Save vesicle information and create a log file
options.logFile = 'output/relaxationManyVes.log';
% Name of log file for saving messages
options.dataFile = 'output/relaxationManyVesData.bin';
% Name of binary data file for storing vesicle information

options.saveError = true;
options.errorFile = 'output/relaxationManyVesError.bin';
% Name of binary data file for storing truncation errors after each step

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

oc = curve;
Xouter = oc.initConfig(prams.N,'nv',1,...
  'reducedArea',0.75,...
  'center',[0;0],...
  'angle',pi/3,...
  'scale',2.0);

centerx = [-1.8 1.0 0.8 -2.0 2.2];
centery = [-0.9 -1 1.5 -3.7 3.5];
ang = [pi/4 pi/5 0 pi/3 -pi/3];
ra = 0.8;
scale = 0.6;
Xinner = oc.initConfig(prams.N,'nv',prams.nv-1,...
    'reducedArea',ra,...
    'center',[centerx;centery],...
    'angle',ang,...
    'scale',scale);
X = [Xouter Xinner];

%clf
%plot(X(1:end/2,:),X(end/2+1:end,:),'r')
%axis equal
%pause

Xfinal = Ves2D(X,[],prams,options);
