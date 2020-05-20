%function [] = shear1Ves(fluxCoeff,farFieldSpeed,fileName)
%clear all; clc

fprintf('Simple elliptical vesicle in a shear flow.\n');
fprintf('First-order semi-implicit time stepping.\n');

% Physics parameters
prams.N = 128;               % points per vesicle
prams.nv = 1;               % number of vesicles
prams.T = 1e4;               % time horizon (two tumbling)
prams.m = 500;              % number of time steps
prams.kappa = ones(prams.nv,1);   % bending coefficient
prams.viscCont = ones(prams.nv,1);         % viscosity contrast
options.farField = 'shear'; % background velocity
options.farFieldSpeed = 1e-1;
%options.farFieldSpeed = farFieldSpeed;
options.vesves = 'implicit';
% Discretization of vesicle-vesicle interactions.
% Either 'explicit' or 'implicit'
options.inextens = 'method1';
options.near = true;        % near-singular integration
options.fmm = false;
options.verbose = true;
options.semipermeable = true;
prams.gmresMaxIter = 3*prams.N;
prams.gmresTol = 1e-8;
prams.errorTol = 1000;
prams.fluxCoeff = 1e-3;
%prams.fluxCoeff = fluxCoeff;
options.fluxShape = 1;

options.adhesion = false;
prams.adRange = 0.4;
prams.adStrength = 100;

% TIME ADAPTIVITY (parameters for new implementation)
options.timeAdap = true;

prams.rtolArea = 1e10;
prams.rtolLength = 1e-1;
prams.dtMax = 1e1;
prams.dtMin = 1e-4;
prams.betaUp = 1.2;
prams.betaDown = 0.5;
prams.alpha = 0.9;

options.orderGL = 2;
options.nsdc = 1;
options.expectedOrder = 2;

options.expForce = false;

% Plot on-the-fly
options.usePlot = true;
options.axis = [-5 5 -5 5];
% Save vesicle information and create a log file
%options.logFile = ['output/' fileName '/shear1Ves.log'];
%options.logFile = ['output/' fileName '.log'];
%options.logFile = ['output/shear1Ves.log.' num2str(10*fluxCoeff)];
options.logFile = 'output/shear1Ves.log';
% Name of log file for saving messages
%options.dataFile = ['output/' fileName '/shear1VesData.bin'];
%options.dataFile = ['output/' fileName
%options.dataFile = ['output/shear1VesData.bin.' num2str(10*fluxCoeff)];
options.dataFile = 'output/shear1VesData.bin';
% Name of binary data file for storing vesicle information

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

%if ~exist(['output/' fileName],'dir')
%  mkdir(['output/' fileName]);
%end

oc = curve;
ra = 0.65;
%centerx = [-2.5 2.5];
%centery = zeros(1,2);
centerx = 0;
centery = 0;
ang = pi/2;
scale =  sqrt(ra);
X = oc.initConfig(prams.N,...
    'nv', prams.nv, ...
    'reducedArea',ra,...
    'angle',ang,...
    'center',[centerx;centery],...
    'scale',scale);

[~,~,L] = oc.geomProp(X);

%X = oc.initConfig(prams.N,...
%    'reducedArea',ra,...
%    'angle',pi/2,...
%    'center',[0;0],...
%    'scale',scale,'folds',15,'star');
% Initial configuration of reduce area 0.65 and aligned
%ymax = max(X(end/2+1:end));
%X = X/ymax*3; % make maximum y value equal to 3

%theta = (0:prams.N-1)'*2*pi/prams.N;
%X = [cos(theta);3*sin(theta)];

Xfinal = Ves2D(X,[],prams,options);
% Run vesicle code

