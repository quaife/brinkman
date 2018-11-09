function [] = extensional2Ves

clear all; clc

fprintf('Two elliptical vesicles in a extensional flow.\n');

% Physics parameters
prams.N = 512;               % points per vesicle
prams.nv = 2;               % number of vesicles
prams.T = 600;               % time horizon (two tumbling)
prams.m = 2000;             % number of time steps
prams.kappa = [1e-1 1e-1];         % bending coefficient
prams.viscCant = [1 1];         % viscosity contrast
options.farField = 'extensional'; % background velocity
options.farFieldSpeed = 1e-2;
aptions.order = 1;          % time stepping order
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
prams.adRange = 2e-1;
prams.adStrength = 1e0;

% TIME ADAPTIVITY (parameters for new implementation)
options.timeAdap = true;

prams.rtolArea = 1e-2;
prams.rtolLength = 1e-2;
prams.betaUp = 1.2;
prams.betaDown = 0.5;
prams.alpha = 0.9;

options.orderGL = 2;
options.nsdc = 1;
options.expectedOrder = 2;

% Plot on-the-fly
options.usePlot = false;
options.axis = [-4 4 -2 2];
options.track = false;
% Save vesicle information and create a log file
options.logFile = 'output/extensional2Ves.log';
% Name of log file for saving messages
options.dataFile = 'output/extensional2VesData.bin';
% Name of binary data file for storing vesicle information

options.saveError = true;
options.errorFile = 'output/extensional2VesError.bin';
% Name of binary data file for storing truncation errors after each step

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

oc = curve;
centerx = [-2 2];
centery = zeros(1,2);
ang = 0*ones(2,1);
ra = 0.6;
scale = 0.5*sqrt(ra);
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

