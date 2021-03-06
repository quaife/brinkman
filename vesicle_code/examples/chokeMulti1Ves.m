clear all;
%clc; clf;
fprintf('One elliptical vesicles in a constricted tube.\n');
fprintf('Second-order semi-implicit time stepping.\n');
fprintf('Implicit vesicle-vesicle interactions.\n');
fprintf('Implicit vesicle-boundary interactions.\n');

%format long
prams.N = 32;                 % points per vesicle
prams.nv = 1;                  % number of vesicles
prams.T = 100;                  %time horizon
prams.m = 10000;                % number of time steps
prams.Nbd = 4*1024;               % number of points on solid wall
prams.nvbd = 1;                % number of components to solid walls
prams.kappa = 1e0;            % bending coefficient
prams.viscCont = 1*ones(prams.nv,1);            % viscosity contrast
prams.gmresTol = 1e-10;        % GMRES tolerance
prams.errorTol = 8e-1;
% Maximum relative error in area and length before the simulation
% is stopped
prams.dtMin = 1e-3;
prams.dtMax = 1e0;

options.semipermeable = false;
prams.fluxCoeff = 0e-3;
options.fluxShape = 1; % constant value

options.farField = 'chokeMulti';      % Constricted domain
options.farFieldSpeed = 1;
% background velocity or solid walls (which comes with a default
% boundary condition)
options.vesves = 'implicit';
% Discretization of vesicle-vesicle and vesicle-boundary 
% intearctions.  Either 'explicit' or 'implicit'
options.near = true;        % near-singular integration
options.fmm = ~true;
options.fmmDLP = ~true;
% use FMM to compute single-layer potential
options.confined = true;   % confined or unbounded geometry
% Used for geometry and boundary conditins of solid walls
options.usePlot = false;     % View vesicles during run
options.axis = [-1.3 11.3 -0.1 5.6]; 
% Axis for plots if usePlot = true
options.saveData = true;    
% Save vesicle information and create a log file
options.logFile = 'output/chokeMulti1Ves2.log';
% Name of log file
options.dataFile = 'output/chokeMulti1Ves2Data.bin';
% Name of binary data file for storing vesicle information
options.verbose = true;
% Decides how much info is written to console
options.profile = false;    % Profile code
options.collision = false;   % Collision detection

options.nsdc = 1;
options.orderGL = 2;

options.timeAdap = true;
prams.rtolArea = 1e-2;
prams.rtolLength = 1e-2;
prams.betaUp = 1.5;
prams.betaDown = 0.5;
prams.alpha = 0.9;

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

oc = curve;
X = oc.initConfig(prams.N,'nv',prams.nv,...
   'angle',0,...
   'scale',0.05,...
   'center',[[0;5]],...
   'reducedArea',0.65);

% load chokeIC
Xwalls = oc.initConfig(prams.Nbd,options.farField);
%clf;
%plot(Xwalls(1:end/2),Xwalls(end/2+1:end),'k');
%hold on;
%axis equal;
%plot(X(1:end/2),X(end/2+1:end),'r');
%pause

%pressx = linspace(-18,18,1000)';
%pressy = 0*ones(size(pressx));
pressx = [-60;60];
pressy = [0;0];
pressTar = [pressx;pressy];

Xfinal = Ves2D(X,Xwalls,prams,options);
% Run vesicle code


