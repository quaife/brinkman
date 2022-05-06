function [] = confined1Ves_continuation(fluxCoeff,farFieldSpeed,concentration,...
                               shortax, scaleL, fileName, farFieldFlow,...
                               wallGeometry, vesCenter, dt, Nbd, N, ...
                               b_min, b_max, a, oddeven, T, restart,rwhere)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code has been developed using the methods described in "Dynamics
%of multicomponent vesicles in a viscous fluid" by Sohn, Tseng, Li
%Voigt, &  Lowengrub JCP 2013

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Initialize parameters and options %%%%%%%%%%%%%%%%%%%%%%%
% TODO: SOME OF THESE ARE MORE OPTIONS THAN PARAMETERS

params.N = N; % points on vesicle
params.Nbd = Nbd; % points on the solid wall
params.dt = dt; % time step size
params.T = T; % time horizon
params.saveRate = 1; % ouptut frequency
params.concentra = concentration; % constant, initial concentration of 
                                  % lipid species
params.oddeven = oddeven; % flag for initial lipid species profile
params.shortax = shortax; % short axis length
params.scaleL = scaleL;
params.farFieldFlow = farFieldFlow;   
% Available options: 'relaxation', 'shear', 'parabolic, 'extensional'
% 'tube', 'choke', 'doublechoke', 'contracting'
params.wallGeometry = wallGeometry;
params.vesGeometry = 'ellipse';
% Available options: 'relaxation', 'shear', 'parabolic, 'extensional'
% 'tube', 'choke', 'doublechoke', 'contracting', 'tube' 
params.shearRate = farFieldSpeed; % shear rate
params.viscosityInside = 1.0;
params.viscosityOutside = 1.0;
params.bendsti = b_max; % maximum bending stiffness
params.bendratio = b_min/b_max; % ratio between max and min bending stiffness
params.consta = a; % parameter 'a' in the Cahn-Hilliard energy
params.nloop = 20e2*5; % number of time steps of Cahn-Hilliard to be taken at 
                   %each time step of the hydrodynamics
params.epsch = 0.04; % small parameter  in the double-well potential 
params.gmresTol = 1e-10; %GMRES tolerance
params.gmresMaxIter = params.N; %maximum number of GMRES iterations
params.SPcoeff = fluxCoeff; %semi-permeable coefficient
params.vesCenter = vesCenter;%[0;0.1];
params.geomCenter = [0;0];
params.angle = 0;%pi/6; % The tracking point is programmed to be at 0,0.  
                      % Rotate the vesicle counter-clockwise to keep 
                      % desired center.  

options.confined = false; %param for now to pass into tstep, change later
options.verbose = false;  % write data to console
options.saveData = true; % save the data
options.usePlot = true;  % plot the data
options.axis = [-25 25 -4 4];
options.dataFile = true; % data file name
options.logFile = true;  % log file name

if restart
    if rwhere
        options.logFile = ['output/' fileName '_contFromMid.log'];
        options.dataFile = ['output/' fileName '_contFromMid.bin'];
    else
        options.logFile = ['output/' fileName '_contFromEnd.log'];
        options.dataFile = ['output/' fileName '_contFromEnd.bin'];
    end
else
    options.logFile = ['output/' fileName '_cont.log'];
    options.dataFile = ['output/' fileName '_cont.bin'];
end
%tic
%profile on
ves = Ves2D_Continuation(params,options,fileName, restart, rwhere);
%profile off
%profile viewer
%toc
end