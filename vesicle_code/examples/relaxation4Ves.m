clear all; clc

fprintf('Four Simple elliptical vesicles in a relaxation flow.\n');

% Physics parameters
prams.N = 128;               % points per vesicle
prams.nv = 5;               % number of vesicles
prams.T = 10;               % time horizon (two tumbling)
prams.m = 200;             % number of time steps
prams.kappa = 1e-1*ones(prams.nv,1); % bending coefficient
prams.viscCont = ones(prams.nv,1);   % viscosity contrast
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
options.adhesion = false;
prams.adRange = 2e-1;
prams.adStrength = 2;

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
options.axis = [-1 1 -1 1];
options.track = false;
% Save vesicle information and create a log file
options.logFile = 'output/relaxation4Ves.log';
% Name of log file for saving messages
options.dataFile = 'output/relaxation4VesData.bin';
% Name of binary data file for storing vesicle information

options.saveError = true;
options.errorFile = 'output/relaxation4VesError.bin';
% Name of binary data file for storing truncation errors after each step

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

oc = curve;
centerx = 0.44*[1.0 -1.0 -1.0 1.0 0];
centery = 0.44*[1.0 1.0 -1.0 -1.0 0];
ang = [-0.5 pi/2-0.5 -0.5 pi/2-0.5 0.2];
ra = 0.65;
scale = 0.21*sqrt(ra);
X = oc.initConfig(prams.N,'nv',prams.nv,...
    'reducedArea',ra,...
    'center',[centerx;centery],...
    'angle',ang,...
    'scale',scale);
% Initial configuration of reduce area 0.65 and aligned
clf;
plot(X(1:end/2,:),X(end/2+1:end,:),'r')
axis equal
hold on
theta = linspace(0,2*pi,1000);
plot(exp(1i*theta),'k-')
pause


Xfinal = Ves2D(X,[],prams,options);
