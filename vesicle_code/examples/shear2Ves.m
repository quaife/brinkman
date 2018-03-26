clear all; clc

fprintf('Simple elliptical vesicle in a shear flow.\n');
fprintf('First-order semi-implicit time stepping.\n');

% Physics parameters
prams.N = 96;               % points per vesicle
prams.nv = 2;               % number of vesicles
prams.T = 40;               % time horizon (two tumbling)
prams.m = 2000;             % number of time steps
prams.kappa = [1e-1 1e2];         % bending coefficient
prams.viscCont = [1 1];         % viscosity contrast
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
options.correctShape = false;
options.adhesion = true;
prams.adRange = 4e-1;
prams.adStrength = 2e-1;

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
options.axis = [-10 10 -5 5];
options.track = false;
% Save vesicle information and create a log file
options.logFile = 'output/shear2Ves.log';
% Name of log file for saving messages
options.dataFile = 'output/shear2VesData.bin';
% Name of binary data file for storing vesicle information

options.saveError = true;
options.errorFile = 'output/shear2VesError.bin';
% Name of binary data file for storing truncation errors after each step

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

oc = curve;
ang = pi/2*ones(2,1);
centerx = [-3.0 3.0];
centery = zeros(1,2);
X = oc.initConfig(prams.N,'nv',prams.nv,...
    'reducedArea',0.65,...
    'center',[centerx;centery],...
    'angle',ang);
% Initial configuration of reduce area 0.65 and aligned
%clf;
%plot(X(1:end/2,:),X(end/2+1:end,:),'r')
%axis equal
%pause

Ves2D(X,[],prams,options);
% Run vesicle code

