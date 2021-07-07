classdef tstep < handle
% This class implements the velocity evaulation schemes and the time
% stepping methods


properties

walls; % outer wall structure
confined; %confined option T or F
shearRate; % background shear rate
farFieldFlow; % far field flow type
bendsti; % maximum bending stiffness
bendratio; % ratio between max and min bending stiffness
SPc; %Semi permeability coefficient
kmatrix; % matrix for doing odd-even integration
viscIn; %viscosity inside the vesicle
viscOut; %viscosity outside the vesicle
gmresTol; %tolerance for GMRES 
gmresMaxIter; %maximum number of iterations GMRES
saveRate; % how often the solution will be saved
DLPmat; %double-layer Stokes potential matrix
end % properties

methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = tstep(params,options,ves,walls)
%o = tstep(prams,ves,walls) is a constructor that initializes the class.
%Take all elements of prams needed by the time stepper
o.confined = options.confined;
o.shearRate = params.shearRate;
o.farFieldFlow = params.farFieldFlow;
o.bendsti = params.bendsti; %maximum bending stiffness
o.bendratio = params.bendratio; %ratio between max and min bending 
                                %stiffness
o.viscIn = params.viscosityInside;
o.viscOut = params.viscosityOutside;
o.SPc = params.SPcoeff; %Semi-permeability coefficient
o.gmresTol = params.gmresTol; %GMRES tolerance
o.gmresMaxIter = params.gmresMaxIter; %maximum number of GMRES iterations
%Build poten class op
op = poten(params.N);
%construct the odd/even matrix
o.kmatrix = op.oddEvenMatrix; 
oc = curve(params.N);
o.saveRate = params.saveRate;
if o.confined
    o.walls = walls;
    o.DLPmat = op.StokesDLP(o.walls);
else
    o.walls = [];
    o.DLPmat = [];
end

end % tstep: constructor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [uinf] = bgFlow(o, X, farFieldSpeed, farFieldFlow)
% function to compute the far field velocity using specified flow type
% varargin. If varagin is not specified, default flow is quiescent
% (relaxation)
%     relaxation:     (0,0)
%     shear:          (ky,0)
%     extensional:    (-x,y)
%     parabolic:      (k(1-y/r)^2,0) - r is hardcoded for now.

N = length(X)/2; % number of points on vesicle or solid wall
oc = curve(N);
[x,y] = oc.getXY(X);

if strcmp(farFieldFlow,'relaxation')
    uinf = zeros(2*N, 1);
elseif strcmp(farFieldFlow,'shear')
    uinf = farFieldSpeed*[y;zeros(N,1)];
elseif strcmp(farFieldFlow, 'parabolic')
    [ra,A,~] = oc.geomProp(X);
    R0 = sqrt(A/pi);
    W = 10*R0;
    W = 3.5;
    uinf = farFieldSpeed*[(1-(y/W).^2);zeros(N,1)];
elseif strcmp(farFieldFlow,'extensional')
    uinf = farFieldSpeed*[-x;y];
elseif strcmp(farFieldFlow, 'tube')
    W = max(y);
    uinf = farFieldSpeed*[(1-(y/W).^2);zeros(N,1)];
elseif strcmp(farFieldFlow, 'choke') || strcmp(farFieldFlow, 'doublechoke')
    W = 3.5;
    uinf = farFieldSpeed*[(1-(y/W).^2);zeros(N,1)];
    uinfx = uinf(1:end/2);
    uinfx(abs(x) < 8) = 0;
    uinf = [uinfx;uinf(end/2+1:end)];
elseif strcmp(farFieldFlow, 'longchoke')
    W = 3.5;
    uinf = farFieldSpeed*[(1-(y/W).^2);zeros(N,1)];
    uinfx = uinf(1:end/2);
    uinfx(abs(x) < 28) = 0;
    uinf = [uinfx;uinf(end/2+1:end)];
elseif strcmp(farFieldFlow, 'contracting')
    uinf = zeros(N,1);
    xmax = max(x);
    xmin = min(x);
    ind = find(x < (xmin + 0.1) | x > (xmax - 0.1));
    ymax = max(y(ind));
    vx = (ymax^2-y(ind).^2)/ymax^2;
    % parabolic flow
    uinf(ind,:) = vx;  
else
    %default flow is relaxed
    uinf = zeros(N, 1);
    
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eta] = etaSolver(o, rhs)
    walls = o.walls;
    op = poten(walls.N);
    N0 = op.StokesN0mat(walls);
    % solve for the density function which satisfies the second-kind
    % Fredholm integral equation
    eta = (0.5*eye(2*walls.N) + o.DLPmat + N0)\rhs;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [uxvel,uyvel,fdotn,eta] = usetself(o,ves,eta)
%Usetself returns the x and y component of the velocity as well as the 
%mechanical force, beta*fdotn.

walls = o.walls; % shorthand for the walls object
x0 = ves.x0; %shorthand for x component of tracking point
y0 = ves.y0; %shorthand for y component of tracking point
N = ves.N; %shorthand for Number of discratization points
L = ves.L; %shorthand for Length of vesicle
theta = ves.theta; %shorthand for opening angle of vesicle 
op = poten(N); %shorthand for poten class
oc = curve(N); %shorthand for curve class
% Compute the curvature
cur = oc.acurv(ves.N,ves.theta,ves.L);
% Reconstruct the vesicle shape
ves.X = oc.recon(ves.N, ves.x0, ves.y0, ves.L, ves.theta);
% Compute Eu and Esigma, equations (13) and (14) 
[Eu,Esigma] = ves.variationsNonStiff;
% Compute the force in eq (33)
%   NOTE: this does not include u_s term.
tau = [[-Esigma.*sin(theta) - Eu.*cos(theta)]; ...
       [Esigma.*cos(theta) - Eu.*sin(theta)]];
% --- Pulling force code - Yuan's project
% x0 = 0.2;
% y0 = 3.0;
% d = 1;
% x = ves.X(1:end/2);
% y = ves.X(end/2+1:end);
% r = sqrt((x - x0).^2 + (y - y0).^2);
% fx = -exp(-r/d).*(x-x0)./r;
% fy = -exp(-r/d).*(y-y0)./r;
% tau = tau + 5*[fx;fy];
% ---------------------------------------- 

%Construct Stokes matrix without the log singularity. 
%ie. StokesMat only contains the rightmost kernel in equation (43)
StokesMat = op.StokesMatrixLogless(ves);
%Evaluate the single-layer poterntial to find the solution of equation (33) 
vel = op.StokesSLP(ves,StokesMat,tau);
%Include the background shear flow which is not stated in equation (33),
%but instead in equations (6) and (7). The background velocity is 0 if the
%geometry is confined.
if ~o.confined
  uinf = o.bgFlow(ves.X, o.shearRate, o.farFieldFlow);
else
  uinf = zeros(size(ves.X));
end
[uinfx, uinfy] = oc.getXY(uinf);
vesVelx = vel(1:N) + uinfx;
vesVely = vel(N+1:end) + uinfy;

if o.confined
  %near singular integration structure for characterizing source points as
  %close or far from target points
  [~,NearW2V] = getZone(walls,ves,2);
  kernel = @op.StokesDLPtar;
  %DLP is the limiting value of the double-layer potential
  DLP = @(X) +0.5*X + o.DLPmat*X;
  %Evaluate the DLP with near singular integration to get the the velocity
  %on the vesicle due to the walls
  wallVel2Ves = op.nearSingInt(walls,eta,DLP,...
                               NearW2V,kernel,kernel,ves,false,false);
  [wallVel2Vesx, wallVel2Vesy] = oc.getXY(wallVel2Ves);
else
  %if the flow is unconfined, there is no velocity on the vesicle due to
  %walls.
  wallVel2Vesx = zeros(size(vesVelx));
  wallVel2Vesy = zeros(size(vesVely));
end

%Calculate v dot n in eq (40), modified for confined flow
vdotn = (vesVelx + wallVel2Vesx).*sin(theta) - ...
        (vesVely + wallVel2Vesy).*cos(theta);
%Calculate v dot s in eq (40), modified for confined flow
vdots = (vesVelx + wallVel2Vesx).*cos(theta) + ...
        (vesVely + wallVel2Vesy).*sin(theta);
%Calculate d/ds(v dot s) in eq (40)
dvdotsds = oc.diffFT(vdots)/L;

%Compute the right hand side of equation (40)
rhs = -(dvdotsds + cur.*vdotn - ves.SPc*cur.*Esigma);
%The velocity components in eq (40) are nonlocal linear functionals of
%lambdaTilde. Solve the linear system for LambdaTilde in (39) using
%GMRES. Each iteration of GMRES requires a solution of Stokes equation.
%LambdaTilde is the lambda with a tilde in eq (39), in case it's not obv.
[lambTil,flag,relres,iter,resvec] = ...
      gmres(@(x) o.matvec40(x,ves,StokesMat),rhs,[],o.gmresTol,...
              o.gmresMaxIter,@(x) o.preconditioner(x));
% Calculate the Fourier derivative of lambdaTilde
dlamTilds = oc.diffFT(lambTil)/L;
% We can now compute the traction jump in first part of equation (39).
% This comes from applying the product rule and using Frenet-Seret.
tracJump = [(+lambTil.*cur.*sin(theta) - dlamTilds.*cos(theta));...
            (-lambTil.*cur.*cos(theta) - dlamTilds.*sin(theta))];
% Now we add the jump conditions in eq (39) to (33) which is in the 
% old variable tau
tau = tau + tracJump;
if o.confined
  [~,NearV2W] = getZone(ves,walls,2);
  kernel = @op.StokesSLPtar;
  SLP = @(X) op.StokesSLP(ves,StokesMat,X);
  %compute velocity due to the vesicle's traction and evaluate on the
  %outer wall using near singular integration
  vesVel2Wall = op.nearSingInt(ves,tau,SLP,...
                               NearV2W,kernel,kernel,walls,false,true);
  %build right hand side for the double-layer potential solver
  bgFlowWalls = o.bgFlow(walls.X, o.shearRate, o.farFieldFlow);
  rhsWalls = bgFlowWalls - vesVel2Wall;
  eta = o.etaSolver(rhsWalls);
end
%Add in the fdotn term for the semipermeability model
fdotn = tau(1:end/2).*sin(theta)-tau(end/2+1:end).*cos(theta);
%Compute u in equation (31) by evaluating the SLP for the u tilde in 
%equations (38) through (40) adding the contribution of the 
vel = op.StokesSLP(ves,StokesMat,tau);
uxvel = vel(1:N) + wallVel2Vesx + uinfx;
uyvel = vel(N+1:end) + wallVel2Vesy + uinfy;
%store the lagrange multiplier in ves.
ves.ten = lambTil;

end % usetself

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LHS = matvec40(o,Lambda,ves,StokesMat)   
% This function cooresponds to the matvec in equation 40 and returns the
% LHS of eq(40), (u \cdot s)_s + kappa * (u \cdot n) as LHS

theta = ves.theta; %shorthand for opening angle
cur = ves.cur; %shorthand for vesicle curvature
L = ves.L; %shorthand for vesicle length
N = ves.N; %shorthand for number of discretization points 
ulam = ves.viscIn/ves.viscOut; %define viscosity ratio
oc = curve(N); %shorthand for curve class
op = poten(N); %shorthand for poten class
% The the derivative of the rhs of eq(40)
drhsds = oc.diffFT(Lambda)/L;
% Compute the forces in equation (39)
tau = [+Lambda.*cur.*sin(theta) - drhsds.*cos(theta); ...
       -Lambda.*cur.*cos(theta) - drhsds.*sin(theta)];
% Form the velocity on the interface that solves (38) and (39), but
% without the log terms in the kernel
 krhs = StokesMat*tau;
% LogKerneltau1 and LogKerneltau2 are the log kernels integrated against
% the density functions tau1 and tau2.
LogKerneltau1 = op.IntegrateLogKernel(tau(1:N));
LogKerneltau2 = op.IntegrateLogKernel(tau(N+1:end));
% Calculate constants that multiply the weakly singular and regular
% parts of the integral operators
c1 = 1/(4*pi)*L/N;
c2 = -L/(8*pi);
% [utilde1 utilde2] is utilde in equations (39) and (40)
utilde1 = krhs(1:N)*c1 + LogKerneltau1*c2;
utilde2 = krhs(N+1:end)*c1 + LogKerneltau2*c2;
% Build left hand side of equation (40)
udotn = utilde1.*sin(theta) - utilde2.*cos(theta);
udots = utilde1.*cos(theta) + utilde2.*sin(theta);
dudotsds = oc.diffFT(udots)/L;
% LHS is (u \cdot s)_s + kappa * (u \cdot n) in eq (40)
LHS = (dudotsds + cur.*udotn + ves.SPc*cur.^2.*Lambda);

end %matvec40

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = preconditioner(o,rhs)
% preconditioner (equation (63) and (64)) for equation (40) which is the
% Schur complement for a Lagrange multiplier needed to enforce local
% inextensibility

N = numel(rhs);
x = fft(rhs)./[1 1:1:N/2 N/2-1:-1:1]';
x = real(ifft(x));

end % preconditioner

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ves,ux_old,uy_old,L,Ln,dcur0,N1,N2Hat,eta] = FirstSteps(...
        o,ves,walls,params,options,om)
% function [ves,ux_old,uy_old,L,Ln,dcur0,N1,N2Hat,cx,cy] = FirstSteps(...
%         o,ves,walls,params,options,om)
% Calculates the first time step [0,dt] using a first-order Euler method
% to find the vesicle position, curv, velocity, and density function.
% Returns ves, L, Ln and dcur0 to use for higher-order time stepping.     

time = 0; % current time
oc = curve(params.N); % define shorthand for curve class
L = ves.L; % shorthand vesicle length
[~,a_old,l_old] = oc.geomProp(ves.X);

% --------------- center of mass calculations
xnormal =  sin(ves.theta);
ynormal = -cos(ves.theta);
%cx0 = oc.centerOfMass(ves.X, ves.X(1:end/2),ves.L,xnormal);
%cy0 = oc.centerOfMass(ves.X, ves.X(end/2+1:end),ves.L,ynormal);
if o.confined
    % build right hand side for the double-layer potential solver
    bgFlowWalls = o.bgFlow(walls.X, o.shearRate, o.farFieldFlow);
    rhsWalls = bgFlowWalls;
    eta = o.etaSolver(rhsWalls);
else 
    eta = zeros(params.Nbd,1);
end

% Compute the x velocity, y velocity using the initialized concentration
% field. A step of Cahn-Hilliard is not taken until after this step.
%   Missing the u_s term in equation (33)
[ux,uy,fdotn,eta] = o.usetself(ves,eta);
om.initializeFiles(ves.X,ves.rcon,time,[ux;uy],ves.ten)
% put the x-y velocity into the normal and tangential velocity.
theta = ves.theta; %shorthand theta so we don't type ves. a million times
un = ux.*sin(theta) - uy.*cos(theta); % Normal Velocity
ut = ux.*cos(theta) + uy.*sin(theta); % Tangential Velocity
% Update length change over time using a first-order Euler method.  For
% subsequent time steps, will use a multistep method as described in
% equation (60)
%   Forward Euler for the length. The first element of dcur should be zero, 
%   so this is just checking for discretization and round-off errors i.e.
%   Ln = ves.L if params.dt*dcur(1) = 0.
dcur0 = oc.fdcur(ves,un);
Ln = L + params.dt*dcur0(1);
% Get the velocity of the vesicle by its velocity of the tangential
% angle, but without the stiff term. The nonlinear term N_1 defined in
% equation (54) but without the alpha derivative multiplied by the
% integral operator \mathcal{L}. It does include the semi-permeability
% terms using fdotn
N1 = oc.fthetaim(ves,un,ut,fdotn);
% Next, evolve the shape and the phase distribution in Fourier space.
% Fourier series of the derivative of the tangent angle. i.e. Fourier
% series of N_1
fsN1 = fft(N1);
% Fourier derivative of the tangent angle adjusted by a linear function
% so that we are taking the fft of a periodic function
fsTA = fft(theta - (2*pi*(0:ves.N-1)/ves.N)');
% Define the Fourier modes, Nk
Nk = 2*pi*[0 1:ves.N/2 ves.N/2-1:-1:1]';
% rsl is the stiff term in equation (55) that will be integrated
% implicitly using an integrating factor
rsl = ves.bendsti*(Nk/ves.L).^3/4 + ...
      ves.SPc*ves.bendsti*(Nk/ves.L).^4;
% rsln is the next step of rsl 'n' is for 'new' since we'll be using
% multistep
rsln = ves.bendsti*(Nk/Ln).^3/4 + ...
       ves.SPc*ves.bendsti*(Nk/Ln).^4;
% use trapezoid rule to approximate integral in equation (57). This next
% line is exactly equation (58) for the quadrature d1 is now exactly as
% in equation (57). Each Fourier mode has its own integating factor.
% This is the Forward Euler method that is analagous to equation (56)
ek = exp(-(params.dt*(rsl + rsln)/2));
% thetan is now the tangent angle of the new shape after taking a single
% step of Euler with the stiffest term treated implicitly and integrated
% with an integrating factor. Add back the linear function that makes
% theta a function that grows by 2*pi everytime you go completely around
% the shape.
thetan = real(ifft(ek.*(fsTA + params.dt*fsN1))) + ...
      2*pi*(0:ves.N-1)'/ves.N;
% -----  define the lipid species model for u  ----- 
if params.concentra > 0
  % Define the Fourier modes, but scaled by 2*pi. Note that these two
  % vectors will be nearly identical since L \approx Ln by
  % inextensibility form stiffest term that is treated implicitly, but
  % is also linear (and diagonal) in Fourier space
  rk = 2*pi*[0 1:ves.N/2 ves.N/2-1:-1:1]'; 
  % compute integrating factor components 
  rsl = params.epsch*(rk/ves.L).^4*params.consta;
  rsln = params.epsch*(rk/Ln).^4*params.consta;
  % compute the integrating factor for each Fourier mode using the
  % trapezoid rule. d1 is the integrating factor in equation (70)
  ek = exp(-(params.dt/params.nloop*(rsl + rsln)/2));
  % Take small time steps to move the lipid species from time 0 to time
  % dt
  for i=1:params.nloop
    % fncon is the non-linear term N_2 in equation (67)
    N2Hat = oc.frconim(ves,params.epsch,params.consta);
    % fcN2 are the fourier coefficients of N_2 as in equation (68)
    fsN2 = fft(N2Hat);
    % fcLS is the fourier coefficients of the lipid species
    % concentration as in equation (68)
    fsLS = fft(ves.rcon);
    % Compute vesicle concentration using first-order Euler method that
    % is analagous to equation (69)
    ves.rcon = real(ifft(ek.*(fsLS+params.dt/params.nloop*fsN2)));
  end 
  else
    N2Hat = [];
end
%              -----  update ves and area -----
% update ves.L
ves.L = Ln;
% update ves.theta
ves.theta = thetan;
% --------------- center of mass calculations -- Not using right now
  xnormal =  sin(ves.theta);
  ynormal = -cos(ves.theta);
  %cx = oc.centerOfMass(ves.X, ves.X(1:end/2),ves.L,xnormal);
  %cy = oc.centerOfMass(ves.X, ves.X(end/2+1:end),ves.L,ynormal);
% compute the average velocity
% avgux = sum(ux)*ves.L/ves.N;
% avguy = sum(uy)*ves.L/ves.N;
% Xprov = oc.recon(ves.N,0,0,ves.L,ves.theta);
 %cXprovx = oc.centerOfMass(Xprov, Xprov(1:end/2),ves.L,xnormal);
 %cXprovy = oc.centerOfMass(Xprov, Xprov(end/2+1:end),ves.L,ynormal);
 %cx = cx0 + params.dt*avgux;
 %cy = cy0 + params.dt*avguy;
 % ves.X(1:end/2) = Xprov(1:end/2) + cx - cXprovx;
 % ves.X(end/2+1:end) = Xprov(end/2+1:end) + cy - cXprovy;
%---------------------------------------------------------------------
% update the position with Forward Euler. At future time steps,
% second-order Adams-Bashforth will be used
ves.x0 = ves.x0 + params.dt*ux(1);
ves.y0 = ves.y0 + params.dt*uy(1);
if(strcmp(o.farFieldFlow, 'parabolic') || ...
  strcmp(o.farFieldFlow, 'ttube'))
  %pin the vesicle at x = 0
  ves.x0 = ves.x0 - mean(ves.X(1:end/2));
end
% set up variables for timestepping loop
ux_old = ux(1);
uy_old = uy(1);
% compute the new area and length of the vesicle
[~,a_new,l_new] = oc.geomProp(ves.X);
% compute the errors --- NOTE do this in monitor
ea = abs(a_new - a_old)./abs(a_old);
el = abs(l_new - l_old)./abs(l_old);
% plot and write the new data 
% om.plotData(ves.X,time,ea,el,[ux;uy],ves.rcon)
end % FirstSteps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ves = TimeStepLoop(o,ves,walls,params,om,ux_old,uy_old,L,Ln,...
                            dcur0,fntheta,N2Hat,eta_old)
% function ves = TimeStepLoop(o,ves,params,om,ux_old,uy_old,L,Ln,...
%                             dcur0,fntheta,N2Hat,cx0,cy0)
% Main time stepping routine which can be either Euler or multistep. Right
% now it is multistep, to change to Euler search for Euler and uncomment
% appropriate lines.

oc = curve(params.N);
nstep = round(params.T/params.dt); % total number of time steps
% Compute current area and length of the vesicle
[~,a_old,l_old] = oc.geomProp(ves.X);

% Entering time stepping loop
for ktime = 1:nstep
  time = ktime*params.dt;
  % compute the x- and y-components of the velocity. This is the routine
  % that calls GMRES which is used to solve equation (30) 
   [uxvel_loop,uyvel_loop,fdotn,eta] = o.usetself(ves,eta_old);
  % Save the velocity at the first discretization point
  uxvel_new = uxvel_loop(1);
  uyvel_new = uyvel_loop(1);
  % put the x-y velocity into the normal and tangential velocity.
  theta = ves.theta; % shorthand for theta
  % Compute the normal Velocity
  un = uxvel_loop.*sin(theta) - uyvel_loop.*cos(theta); 
  % Compute the tangential Velocity
  ut = uxvel_loop.*cos(theta) + uyvel_loop.*sin(theta); 
  % Update length change over time using a 2nd-order Adams Bashforth
  % method. Described in equation (60) 2nd-order Adams Bashforth for
  % the length. The first element of dcur should be zero, so this is
  % just checking for discretization and round-off errors. L = Ln = Lnn.
  dcur1 = oc.fdcur(ves,un);
  Lnn = Ln + params.dt*(3*dcur1(1)-dcur0(1))/2;
  % Update ves.L  
  ves.L = Lnn;
  % Get the velocity of the vesicle by its velocity of the tangential
  % angle without the stiff term. 
  fnthetan = oc.fthetaim(ves,un,ut,fdotn);
  %         -----  update the lipid species model for u  -----
  % Define the Fourier modes, but scaled by 2*pi. Note that these two
  % vectors will be nearly identical since L \approx Ln by
  % inextensibility form stiffest term that is treated implicitly, but
  % is also linear (and diagonal) in Fourier space
  rk = 2*pi*[0 1:ves.N/2 ves.N/2-1:-1:1]'; 
  % compute the nonlocal model for theta
  rsl = ves.bendsti*abs(rk/L).^3/4 + ...
        ves.SPc*ves.bendsti*(rk/ves.L).^4;
  rsln = ves.bendsti*abs(rk/Ln).^3/4 + ...
         ves.SPc*ves.bendsti*(rk/ves.L).^4;
  rslnn = ves.bendsti*abs(rk/Lnn).^3/4 + ...
          ves.SPc*ves.bendsti*(rk/ves.L).^4;
  % compute the local model for theta using the exponential time
  % integrators described in equations (58) and (59)	 
  d1 = exp(-(params.dt*(rsln+rslnn)/2));
  d2 = exp(-(params.dt*(rsl+rslnn)/2 + params.dt*rsln));
  % Compute the fourier series of fntheta
  fcfntheta = fft(fntheta);
  % Compute the Fourier series of fnthetan
  fcfnthetan = fft(fnthetan);
  % Compute the Fourier series of theta adjusted by a linear function so
  % that we are taking the fft of a periodic function
  fcthetan = fft(theta - 2*pi*(0:ves.N-1)'/ves.N);
  % Add back the linear function that makes theta a function that grows
  % by 2*pi everytime you go completely around the shape. thetann is now
  % the tangent angle of the new shape after taking a single step of a
  % 2nd order multistep method with the stiffest term treated implicitly
  % and integrated with the integrating factors d1 and d2 defined above.
  thetann = real(ifft(d1.*fcthetan + ...
        0.5*params.dt*(3*d1.*fcfnthetan - d2.*fcfntheta))) + ...
        2*pi*(0:ves.N-1)'/ves.N;
  
% ------    define the lipid species model for u  -----   
  if params.concentra > 0
    % Compute integrating factor components      
    rsl = params.epsch*(rk/ves.L).^4*params.consta;
    rsln = params.epsch*(rk/Ln).^4*params.consta;
    rslnn = params.epsch*(rk/Lnn).^4*params.consta;
    % Compute the integrating factors for each Fourier mode using the
    % trapezoid rule. These are dependent on the number of steps taken
    % to evolve the phase field surface.
    d1 = exp(-(params.dt/params.nloop*(rsln+rslnn)/2));
    d2 = exp(-(params.dt/params.nloop*((rsl+rslnn)/2+rsln)));
% ======    Evolve phase field on surface         ===
    for i = 1:params.nloop
      %evolve the phase field on the surface            
      N2Hatn = oc.frconim(ves,params.epsch,params.consta);
      %fcN2 are the fourier coefficients of N2
      fcN2 = fft(N2Hat);
      %fcN2n are the fourier coefficients of N2n
      fcN2n = fft(N2Hatn);
      %fcLSn is the fourier coefficients of the lipid species
      %concentration as in equation (68)
      fcLSn = fft(ves.rcon);
      %Compute rcon at current timstep, rconn
      ves.rcon = real(ifft(d1.*fcLSn + 0.5*params.dt/params.nloop*(3*d1.*...
              fcN2n - d2.*fcN2)));
      %update N2Hat for nloop
      N2Hat = N2Hatn;
    end    
  end
%                          ====================
  % update variables for time stepping loop
  dcur0 = dcur1;
  ves.theta = thetann;
  fntheta = fnthetan;   
%--------- center of mass calculations -- NOT using this section
  xnormal =  sin(ves.theta);
  ynormal = -cos(ves.theta);
  %cx = oc.centerOfMass(ves.X, ves.X(1:end/2),ves.L,xnormal);
  %cy = oc.centerOfMass(ves.X, ves.X(end/2+1:end),ves.L,ynormal);
  avgux = sum(uxvel_loop)*ves.L/ves.N;
  avguy = sum(uyvel_loop)*ves.L/ves.N;
  Xprov = oc.recon(ves.N,0,0,ves.L,ves.theta);    
  %cXprovx = oc.centerOfMass(Xprov, Xprov(1:end/2),ves.L,xnormal);
  %cXprovy = oc.centerOfMass(Xprov, Xprov(end/2+1:end),ves.L,ynormal);
  %cx = cx0 + params.dt*avgux;
  %cy = cy0 + params.dt*avguy;  
  %ves.X(1:end/2) = Xprov(1:end/2) + cx - cXprovx;
  %ves.X(end/2+1:end) = Xprov(end/2+1:end) + cy - cXprovy;
%----------------------------------------------------------------
  % update the single tracker point using Adams Bashforth
  ves.x0 = ves.x0 + 0.5*params.dt*(3*uxvel_new - ux_old);
  ves.y0 = ves.y0 + 0.5*params.dt*(3*uyvel_new - uy_old);
  % update the single tracker point using Forward Euler  
  %ves.x0 = ves.x0 + params.dt*uxvel_new;
  %ves.y0 = ves.y0 + params.dt*uyvel_new;
  % Update X
  ves.X = oc.recon(ves.N,ves.x0,ves.y0,ves.L,ves.theta);
  % CHECKPOINT: check the jacobian of the position
       %jac = oc.diffProp(ves.X);
       %figure(2); clf;
       %semilogy(abs(fftshift(fft(jac))))
       %pause(0.1)  
  % Update the curvature
  ves.cur = oc.acurv(ves.N,ves.theta,ves.L);
  
   if(strcmp(o.farFieldFlow, 'parabolic') || ...
           strcmp(o.farFieldFlow, 'ttube'))
       %pin the vesicle at x = 0
       ves.x0 = ves.x0 - mean(ves.X(1:end/2));
%       ves.X(1:end/2) = ves.X(1:end/2) - mean(ves.X(1:end/2));
   end
   
  % HACK: keep the vesicle centered at (0,0)
       %ves.X(1:end/2) = ves.X(1:end/2) - mean(ves.X(1:end/2));
       %ves.X(end/2+1:end) = ves.X(end/2+1:end) - mean(ves.X(end/2+1:end));
       
  % Print outputs.
  %if mod(ktime,o.saveRate) == 0
  cmod = ktime - o.saveRate*floor(ktime/o.saveRate);
  if cmod == 0
      terminate = om.outputInfo(ves.X,walls,ves.rcon,time,[uxvel_loop;uyvel_loop],ves.ten);
      if terminate
        msg = 'ERROR IN AREA IS TOO LARGE. STOPPING SIMULATION.';
        om.writeMessage(msg,'%s\n');
        break
      end
  end
  % Update ux_old, uy_old, and eta_old for timestepping loop
  ux_old = uxvel_new;
  uy_old = uyvel_new;
  eta_old = eta;
%   cx0 = cx;
%   cy0 = cy;
end

end % TimeStepLoop    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % methods

end % classdef
