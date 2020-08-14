%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code has been developed using the methods described in "Dynamics of 
%multicomponent vesicles in a viscous fluid" by Sohn, Tseng, Li Voigt, &  
%Lowengrub JCP 2010

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Initialize parameters and options %%%%%%%%%%%%%%%%%%%%%%%
% TODO: SOME OF THESE ARE MORE OPTIONS THAN PARAMETERS
params.N = 64; % points on vesicle
params.dt = 1e-3*10; % time step size
params.T = 10; % time horizon
params.outpt = 1e-3; % ouptut frequency
params.concentra = 0; % constant, initial concentration of lipid species
params.oddeven = 0; % flag for initial lipid species profile?
params.shortax = 2.0; % short axis length
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

%Form initial shape that has points equally spaced in arclength. alpha
%has the parameter values that give this parameterization
oc = curve;
alpha = (0:params.N-1)'/params.N;
[alpha,X] = oc.initConfig(params.N,false,'ellipse',...
            'shortax',params.shortax,'parameter',alpha);
        
%Define the initial concentration field
rcon = oc.initialConcentration(params.N,alpha,...
      params.concentra,params.oddeven);
  
%build object for the vesicle but without a band-limited opening angle
ves = capsules(X,rcon,params);
        
ves = Ves2D(ves,alpha,params,options);        
