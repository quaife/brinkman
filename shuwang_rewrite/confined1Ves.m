%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code has been developed using the methods described in "Dynamics
%of multicomponent vesicles in a viscous fluid" by Sohn, Tseng, Li
%Voigt, &  Lowengrub JCP 2013

fluxCoeff = 0;
farFieldSpeed = 400;
concentration = 0;
shortax = 2.7;
scaleL = 0.6;
%fileName = 'Unconfined_Chi0_Scale55_shortax2p7_conc0p3_testRand_ep28';
fileName = 'longChokeTest';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Initialize parameters and options %%%%%%%%%%%%%%%%%%%%%%%
% TODO: SOME OF THESE ARE MORE OPTIONS THAN PARAMETERS

params.N = 4*128; % points on vesicle
params.Nbd = 8*128; % points on the solid wall
params.dt = 1e-5; % time step size
params.T = 0.2; % time horizon
params.saveRate = 1; % ouptut frequency
params.concentra = concentration; % constant, initial concentration of 
                                  % lipid species
params.oddeven = -1; % flag for initial lipid species profile
params.shortax = shortax; % short axis length
params.scaleL = scaleL;
params.farFieldFlow = 'longchoke';   
% Available options: 'relaxation', 'shear', 'parabolic, 'extensional'
% 'tube', 'choke', 'doublechoke', 'contracting'
params.wallGeometry = 'longchoke';
params.vesGeometry = 'ellipse';
% Available options: 'relaxation', 'shear', 'parabolic, 'extensional'
% 'tube', 'choke', 'doublechoke', 'contracting', 'tube' 
params.shearRate = farFieldSpeed; % shear rate
params.viscosityInside = 1.0;
params.viscosityOutside = 1.0;
params.bendsti = 1; % maximum bending stiffness
params.bendratio = 0.1; % ratio between max and min bending stiffness
params.consta = 100; % parameter 'a' in the Cahn-Hilliard energy
params.nloop = 20; % number of time steps of Cahn-Hilliard to be taken at 
                   %each time step of the hydrodynamics
params.epsch = 0.28; % small parameter  in the double-well potential 
params.gmresTol = 1e-10; %GMRES tolerance
params.gmresMaxIter = params.N; %maximum number of GMRES iterations
params.SPcoeff = fluxCoeff; %semi-permeable coefficient
params.vesCenter = [-25;0];%[0;0.1];
params.geomCenter = [0;0];
params.angle = 0;%pi/6; % The tracking point is programmed to be at 0,0.  
                      % Rotate the vesicle counter-clockwise to keep 
                      % desired center.  

options.confined = true; %param for now to pass into tstep, change later
options.verbose = false;  % write data to console
options.saveData = true; % save the data
options.usePlot = true;  % plot the data
options.axis = [-35 0 -4 4];
options.dataFile = true; % data file name
options.logFile = true;  % log file name

options.logFile = ['output/' fileName '.log'];
options.dataFile = ['output/' fileName '.bin'];
%tic
%profile on
ves = Ves2D(params,options);
%profile off
%profile viewer
%toc
