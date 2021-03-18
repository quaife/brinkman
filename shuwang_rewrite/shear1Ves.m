function [] = shear1Ves(fluxCoeff,farFieldSpeed,concentration,fileName)
  fluxCoeff = 0;
  farFieldSpeed = 10;
  concentration =0;
  fileName ='TestingFField';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code has been developed using the methods described in "Dynamics of 
%multicomponent vesicles in a viscous fluid" by Sohn, Tseng, Li Voigt, &  
%Lowengrub JCP 2013

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Initialize parameters and options %%%%%%%%%%%%%%%%%%%%%%%
% TODO: SOME OF THESE ARE MORE OPTIONS THAN PARAMETERS

params.N = 64*4; % points on vesicle
params.dt = 1e-3; % time step size
params.T = 1; % time horizon
params.outpt = 1e-3; % ouptut frequency
params.concentra = concentration; % constant, initial concentration of 
                                  % lipid species
params.oddeven = 0; % flag for initial lipid species profile?
params.shortax = 0.3;%.8146; % 0.75; % short axis length
params.farFieldFlow = 'relaxation'; % Available options: 'relaxation', 
                                    % 'shear', 'parabolic, 'extensional'
params.shearRate = farFieldSpeed; % shear rate
params.viscosityInside = 1.0;
params.viscosityOutside = 1.0;
params.bendsti = 1; % maximum bending stiffness
params.bendratio = 0.1; % ratio between max and min bending stiffness
params.consta = 100; % parameter 'a' in the Cahn-Hilliard energy
params.nloop = 20; % number of time steps of Cahn-Hilliard to be taken at 
                   %each time step of the hydrodynamics
params.epsch = 0.04; % small parameter  in the double-well potential 
params.gmresTol = 1e-10; %GMRES tolerance
params.gmresMaxIter = params.N; %maximum number of GMRES iterations
params.SPcoeff = fluxCoeff; %semi-permeable coefficient
params.center = [0;0];
params.angle = 0; % The tracking point is programmed to be at 0,0.  
                      % Rotate the vesicle counter-clockwise to keep 
                      % desired center.  

options.verbose = false;  % write data to console
options.saveData = true; % save the data
options.usePlot = true;  % plot the data
options.dataFile = true; % data file name
options.logFile = true;  % log file name

options.logFile = ['output/' fileName '.log'];
options.dataFile = ['output/' fileName '.bin'];

ves = Ves2D(params,options);      

end
