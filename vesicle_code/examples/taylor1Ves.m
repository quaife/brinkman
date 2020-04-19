clear all; clc

fprintf('Simple elliptical vesicle in a taylor-Green flow.\n');
fprintf('First-order semi-implicit time stepping.\n');

% Physics parameters
prams.N = 64;               % points per vesicle
prams.nv = 9;               % number of vesicles
prams.T = 10;               % time horizon (two tumbling)
prams.m = 500;              % number of time steps
prams.kappa = 1e-1;         % bending coefficient
prams.viscCont = 1e0;         % viscosity contrast
options.farField = 'taylorGreen'; % background velocity
options.order = 1;          % time stepping order
options.vesves = 'implicit';
% Discretization of vesicle-vesicle interactions.
% Either 'explicit' or 'implicit'
options.inextens = 'method1';
options.near = true;        % near-singular integration
options.fmm = false;
prams.gmresMaxIter = 3*prams.N;
prams.gmresTol = 1e-10;
prams.errorTol = 1;

% TIME ADAPTIVITY (parameters for new implementation)
options.timeAdap = false;

prams.rtolArea = 1e-2;
prams.rtolLength = 1e-2;
prams.betaUp = 1.5;
prams.betaDown = 0.5;
prams.alpha = 0.9;

options.orderGL = 2;
options.nsdc = 0;
options.expectedOrder = 1;

% Plot on-the-fly
options.usePlot = true;
% Save vesicle information and create a log file
options.logFile = 'output/taylor1Ves.log';
% Name of log file for saving messages
options.dataFile = 'output/taylor1VesData.bin';
% Name of binary data file for storing vesicle information

options.axis = [0 pi 0 pi];

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

oc = curve;
ang = ones(prams.nv,1);
centerx = [0.5 pi/2 pi-0.5 0.5 pi/2 pi-0.5 0.5 pi/2 pi-0.5];
centery = [0.5 0.5 0.5 pi/2 pi/2 pi/2 pi-0.5 pi-0.5 pi-0.5];
X = oc.initConfig(prams.N,'nv',prams.nv,...
                          'reducedArea',0.65,...
                          'angle',ang,...
                          'center',[centerx;centery],...
                          'scale',0.15);
% Initial configuration of reduce area 0.65 and aligned

Ves2D(X,[],prams,options);
% Run vesicle code

