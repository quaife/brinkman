%function [] = shear1Ves(fluxCoeff,farFieldSpeed,concentration,fileName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code has been developed using the methods described in "Dynamics of 
%multicomponent vesicles in a viscous fluid" by Sohn, Tseng, Li Voigt, &  
%Lowengrub JCP 2013

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Initialize parameters and options %%%%%%%%%%%%%%%%%%%%%%%
% TODO: SOME OF THESE ARE MORE OPTIONS THAN PARAMETERS
params.N = 64; % points on vesicle
params.dt = 1e-3; % time step size
params.T = 0.5; % time horizon
params.outpt = 1e-3; % ouptut frequency
params.concentra = 0.3; %concentration; % constant, initial concentration of lipid species
params.oddeven = 0; % flag for initial lipid species profile?
params.shortax = .8146/3;%8146; % short axis length
params.shearRate = 0; %farFieldSpeed; % shear rate
params.viscosityInside = 1.0;
params.viscosityOutside = 1.0;
params.bendsti = 1; % maximum bending stiffness
params.bendratio = 0.1; % ratio between max and min bending stiffness
params.consta = 100; % parameter 'a' in the Cahn-Hilliard energy
params.nloop = 200; % number of time steps of Cahn-Hilliard to be taken at 
                   %each time step of the hydrodynamics
params.epsch = 0.04; % small parameter  in the double-well potential 
params.gmresTol = 1e-10; %GMRES tolerance
params.gmresMaxIter = params.N; %maximum number of GMRES iterations
params.SPcoeff = 1; %fluxCoeff; %semi-permeable coefficient

options.saveData = false;
options.verbose = true;  % write data to console
options.saveData = false; % save the data
options.usePlot = true;  % plot the data
options.dataFile = true; % data file name
options.logFile = true;  % log file name

options.logFile = ['output/Chi80_ra08146_beta01_conc031e4d4_dt.log'];
options.dataFile = ['output/Chi80_ra08146_beta01_conc03_dt1e3.bin'];
ves = Ves2D(params,options);      

%end
