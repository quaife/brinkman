clear all; clc

fprintf('Simple elliptical vesicle in a shear flow.\n');
fprintf('First-order semi-implicit time stepping.\n');

% Physics parameters
prams.N = 128;               % points per vesicle
prams.nv = 1;               % number of vesicles
prams.T = 10;               % time horizon (two tumbling)
prams.m = 5*100;              % number of time steps
prams.kappa = ones(prams.nv,1);   % bending coefficient
prams.viscCont = ones(prams.nv,1);         % viscosity contrast
options.farField = 'shear'; % background velocity
options.farFieldSpeed = 2.0;
options.order = 1;          % time stepping order
options.vesves = 'implicit';
% Discretization of vesicle-vesicle interactions.
% Either 'explicit' or 'implicit'
options.inextens = 'method1';
options.near = true;        % near-singular integration
options.fmm = false;
options.verbose = true;
options.antiAlias = false;
options.semipermeable = true;
prams.gmresMaxIter = 3*prams.N;
prams.gmresTol = 1e-8;
prams.errorTol = 1000;
prams.fluxCoeff = 1.0;
options.adhesion = false;
prams.adRange = 0.4;
prams.adStrength = 100;

% ADD-ONS
options.alignCenterAngle = false;
options.correctShape = false;
options.reparameterization = false;

% TIME ADAPTIVITY (parameters for new implementation)
options.timeAdap = false;

prams.rtolArea = 1e10;
prams.rtolLength = 1e-2;
if 1
  prams.dtMax = 2;
  prams.dtMin = 1e-4;
  prams.betaInc = 1e-1;
  prams.betaDec = 5e-1;
  prams.betaUp = 1.2;
  prams.betaDown = 0.5;
  prams.alpha = 0.9;
end

options.orderGL = 2;
options.nsdc = 1;
options.expectedOrder = 2;

options.expForce = false;

% Plot on-the-fly
options.usePlot = true;
options.axis = [-5 5 -5 5];
options.track = false;
% Save vesicle information and create a log file
options.logFile = 'output/shear1Ves.log';
% Name of log file for saving messages
options.dataFile = 'output/shear1VesData.bin';
% Name of binary data file for storing vesicle information

options.saveError = true;
options.errorFile = 'output/shear1VesError.bin';
% Name of binary data file for storing truncation errors after each step

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

theta = (0:prams.N-1)'*2*pi/prams.N;
%prams.fluxShape = 0*sin(theta);
prams.fluxShape = ones(prams.N,1); %flux shape 1
%fluxWidth = 1;
%prams.fluxShape = exp(-fluxWidth*(theta - pi/2).^2) + ...
%                  exp(-fluxWidth*(theta - 3*pi/2).^2); % flux shape 2
%prams.fluxShape = exp(-fluxWidth*(theta - pi/2).^2); % flux shape 3
% set up the distribution for the flux

oc = curve;
ra = 0.55;
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

