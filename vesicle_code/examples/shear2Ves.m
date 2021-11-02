clear all; clc

fprintf('Two elliptical vesicles in a shear flow.\n');

% Physics parameters
prams.N = 128;               % points per vesicle
prams.nv = 2;               % number of vesicles
prams.T = 15000;               % time horizon (two tumbling)
prams.m = 10000;             % number of time steps
prams.kappa = [1 1];         % bending coefficient
prams.viscCant = [1 1];         % viscosity contrast
options.farField = 'shear'; % background velocity
options.farFieldSpeed = 0.0013;
options.vesves = 'implicit';
% Discretization of vesicle-vesicle interactions.
% Either 'explicit' or 'implicit'
options.inextens = 'method1';
options.near = true;        % near-singular integration
options.fmm = false;
options.verbose = true;
prams.gmresMaxIter = 3*prams.N;
prams.gmresTol = 1e-10;
prams.errorTol = 1;

% ADD-ONS
options.adhesion = false;
prams.adRange = 1e-1;
prams.adStrength = 7e-1;

% TIME ADAPTIVITY (parameters for new implementation)
options.timeAdap = true;

prams.rtolArea = 1e-3;
prams.rtolLength = 1e-3;
prams.betaUp = 1.2;
prams.betaDown = 0.5;
prams.alpha = 0.9;

options.orderGL = 2;
options.nsdc = 1;
options.expectedOrder = 2;

% Plot on-the-fly
options.usePlot = true;
options.axis = [-20 20 -20 20];
% Save vesicle information and create a log file
options.logFile = 'output/shear2Ves.log';
% Name of log file for saving messages
options.dataFile = 'output/shear2VesData.bin';
% Name of binary data file for storing vesicle information

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

oc = curve;
centerx = [-9 9];
centery = zeros(1,2);
ang = pi/2*ones(2,1);
ra = 0.98;
scale = 3*sqrt(ra);
X = oc.initConfig(prams.N,'nv',prams.nv,...
    'reducedArea',ra,...
    'center',[centerx;centery],...
    'angle',ang,...
    'scale',scale);
% Initial configuration of reduce area 0.65 and aligned
%clf;
%plot(X(1:end/2,:),X(end/2+1:end,:),'r')
%axis equal
%pause

Xfinal = Ves2D(X,[],prams,options);
% Run vesicle code

