clear all; clc

fprintf('Simple elliptical vesicle in a taylor-Green flow.\n');
fprintf('First-order semi-implicit time stepping.\n');

% Physics parameters
prams.N = 64;               % points per vesicle
prams.nv = 1;               % number of vesicles
prams.T = 50;               % time horizon (two tumbling)
prams.m = 200;              % number of time steps
prams.kappa = 1e-1;         % bending coefficient
prams.viscCont = 1e1;         % viscosity contrast
options.farField = 'rotate'; % background velocity
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
options.timeAdap = true;

prams.rtolArea = 1e-2;
prams.rtolLength = 1e-2;
prams.dtMax = 2;
prams.dtMin = 1e-4;
prams.betaUp = 1.2;
prams.betaDown = 0.5;
prams.alpha = 0.9;

options.orderGL = 2;
options.nsdc = 1;
options.expectedOrder = 1;

% Plot on-the-fly
options.usePlot = true;
% Save vesicle information and create a log file
options.logFile = 'output/taylor1Ves.log';
% Name of log file for saving messages
options.dataFile = 'output/taylor1VesData.bin';
% Name of binary data file for storing vesicle information

options.axis = [-2 2 -2 2];

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

oc = curve;
X = oc.initConfig(prams.N,'reducedArea',0.65,...
                          'angle',pi/2,...
                          'center',[1;1],...
                          'scale',0.04);
% Initial configuration of reduce area 0.65 and aligned

Ves2D(X,[],prams,options);
% Run vesicle code

