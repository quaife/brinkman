clear all;
%clc; clf;
fprintf('One elliptical vesicles in a constricted tube.\n');
fprintf('Second-order semi-implicit time stepping.\n');
fprintf('Implicit vesicle-vesicle interactions.\n');
fprintf('Implicit vesicle-boundary interactions.\n');

%format long
prams.N = 192;                 % points per vesicle
prams.nv = 1;                  % number of vesicles
prams.T = 0.1;                  %time horizon
prams.m = 12000;                % number of time steps
prams.Nbd = 1024;               % number of points on solid wall
prams.nvbd = 1;                % number of components to solid walls
prams.kappa = 1e0;            % bending coefficient
prams.viscCont = 1*ones(prams.nv,1);            % viscosity contrast
prams.gmresTol = 1e-10;        % GMRES tolerance
prams.errorTol = 8e-1;
% Maximum relative error in area and length before the simulation
% is stopped
prams.dtMin = 1e-5;
prams.dtMax = 1e0;

options.semipermeable = true;
prams.fluxCoeff = 1e-3;
options.fluxShape = 1; % constant value

options.farField = 'chokeLonger';      % Constricted domain
options.farFieldSpeed = 4200;
%options.farFieldSpeed = 4.2;
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
options.usePlot = true;     % View vesicles during run
options.axis = [-110.5 110.5 -12.5 12.5]; 
% Axis for plots if usePlot = true
options.saveData = true;    
% Save vesicle information and create a log file
options.logFile = 'output/choke1Ves.log';
% Name of log file
options.dataFile = 'output/choke1VesData.bin';
% Name of binary data file for storing vesicle information
options.verbose = true;
% Decides how much info is written to console
options.profile = false;    % Profile code
options.collision = false;   % Collision detection

options.nsdc = 1;
options.orderGL = 2;

options.timeAdap = true;
prams.rtolArea = 1e20;
prams.rtolLength = 1e-2;
prams.betaUp = 1.5;
prams.betaDown = 0.5;
prams.alpha = 0.9;

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

oc = curve;
X = oc.initConfig(prams.N,'nv',prams.nv,'angle',-pi/3,...
   'scale',1.44,...
   'center',[[-45;5]],'reducedArea',0.65);
%X = oc.initConfig(prams.N,'nv',prams.nv,'angle',pi/2,...
%   'scale',1,...
%   'center',[[-45;0]],'reducedArea',0.65);

% load chokeIC
Xwalls = oc.initConfig(prams.Nbd,options.farField);

%pressx = linspace(-18,18,1000)';
%pressy = 0*ones(size(pressx));
pressx = [-60;60];
pressy = [0;0];
pressTar = [pressx;pressy];

%Xfinal = Ves2D(X,Xwalls,prams,options,pressTar);
Xfinal = Ves2D(X,Xwalls,prams,options);
% Run vesicle code


