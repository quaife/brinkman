%function [] = shear1Ves(fluxCoeff,farFieldSpeed,concentration,shortax, scaleL, fileName)
%   fluxCoeff = 0;
%   farFieldSpeed = 10;
%   concentration =0.3;
%   fileName ='extensional_RAp3_Conc0p3_Chi10_beta0_n500_dt0en5';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code has been developed using the methods described in "Dynamics of 
%multicomponent vesicles in a viscous fluid" by Sohn, Tseng, Li Voigt, &  
%Lowengrub JCP 2013

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Initialize parameters and options %%%%%%%%%%%%%%%%%%%%%%%
% TODO: SOME OF THESE ARE MORE OPTIONS THAN PARAMETERS

params.N = 512; % points on vesicle
params.Nbd = [];
params.dt = 1e-4; % time step size
params.T = 10; % time horizon
params.saveRate = 1; % ouptut frequency
params.concentra = 0.0; % constant, initial concentration of 
                                  % lipid species
params.oddeven = -1; % flag for initial lipid species profile?
params.shortax = 3.45; % short axis length
params.scaleL = 0.49;
% Available options: 'relaxation', 'shear', 'parabolic, 'extensional'
params.farFieldFlow = 'relaxation'; 
params.shearRate = 0; % shear rate
params.viscosityInside = 1.0;
params.viscosityOutside = 1.0;
params.bendsti = 1; % maximum bending stiffness
params.bendratio = 1; % ratio between max and min bending stiffness
params.consta = 100; % parameter 'a' in the Cahn-Hilliard energy
params.nloop = 20; % number of time steps of Cahn-Hilliard to be taken at 
                   %each time step of the hydrodynamics
params.epsch = 0.04; % small parameter  in the double-well potential 
params.gmresTol = 1e-10; %GMRES tolerance
params.gmresMaxIter = params.N; %maximum number of GMRES iterations
params.SPcoeff = 0.0; %semi-permeable coefficient
params.vesCenter = [0;0];
params.geomCenter = [0;0];
params.vesGeometry = 'ellipse';
params.wallGeometry = 'longchoke';
params.angle = 0; % The tracking point is programmed to be at 0,0.  
                      % Rotate the vesicle counter-clockwise to keep 
                      % desired center.  

options.confined = false; %param for now to pass into tstep, change later
options.verbose = false;  % write data to console
options.saveData = true; % save the data
options.usePlot = true;  % plot the data
options.axis = [-2 2 -2 2];
options.dataFile = true; % data file name
options.logFile = true;  % log file name

fileName = 'relaxation1Ves';
options.logFile = ['output/' fileName '.log'];
options.dataFile = ['output/' fileName '.bin'];

ves = Ves2D(params,options);      

