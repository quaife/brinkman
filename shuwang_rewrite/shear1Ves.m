%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code has been developed using the methods described in "Dynamics of 
%multicomponent vesicles in a viscous fluid" by Sohn, Tseng, Li Voigt, &  
%Lowengrub JCP 2010

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Initialize parameters and options %%%%%%%%%%%%%%%%%%%%%%%
% TODO: SOME OF THESE ARE MORE OPTIONS THAN PARAMETERS
params.N = 96; % points on vesicle
params.dt = 1e-2; % time step size
params.T = 10; % time horizon
params.outpt = 1e-3; % ouptut frequency
params.concentra = 0; % constant, initial concentration of lipid species
params.oddeven = 0; % flag for initial lipid species profile?
params.shortax = 3.0; % short axis length
params.shearRate = 1; % shear rate
params.viscosityInside = 1.0;
params.viscosityOutside = 1.0;
params.bendsti = 1; % maximum bending stiffness
params.bendratio = 1; % ratio between max and min bending stiffness
params.consta = 100; % parameter 'a' in the Cahn-Hilliard energy
params.nloop = 20; % number of time steps of Cahn-Hilliard to be taken at 
                   %each time step of the hydrodynamics
params.epsch = 5e-2; % small parameter  in the double-well potential 
params.gmresTol = 1e-10; %GMRES tolerance
params.gmresMaxIter = params.N; %maximum number of GMRES iterations

options.saveData = true;
options.verbose = true;  % write data to console
options.saveData = true; % save the data
options.usePlot = true;  % plot the data
options.dataFile = true; % data file name
options.logFile = true;  % log file name
options.logFile = 'output/relaxation1VesTest.log';
options.dataFile = 'output/relaxation1VesDataTest.bin';
 
ves = Ves2D(params,options);        
