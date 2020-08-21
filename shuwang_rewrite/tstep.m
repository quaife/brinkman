classdef tstep < handle
% This class implements the velocity evaulation schemes and the time
% stepping methods


properties

ves; % vesicle structure
shearRate; % background shear rate
bendsti; % maximum bending stiffness
bendratio; % ratio between max and min bending stiffness
kmatrix; % matrix for doing odd-even integration
viscIn;
viscOut;
gmresTol;
gmresMaxIter;

end % properties

methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = tstep(params,ves)
%o = tstep(prams,ves) is a constructor that initializes the class. 
%Take all elements of prams needed by the time stepper
o.ves = ves;
o.shearRate = params.shearRate;
o.bendsti = params.bendsti;
o.bendratio = params.bendratio;
o.viscIn = params.viscosityInside;
o.viscOut = params.viscosityOutside;
op = poten(params.N);
o.kmatrix = op.oddEvenMatrix; %construct the odd/even matrix
o.gmresTol = params.gmresTol; %GMRES tolerance
o.gmresMaxIter = params.gmresMaxIter; %maximum number of GMRES iterations

end % tstep: constructor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [uxvel, uyvel] = usetself(o)
%Usetself returns the x and y component of the velocity.

ves = o.ves; %shorthand for ves object
x0 = ves.x0; %shorthand for x component of tracking point
y0 = ves.y0; %shorthand for y component of tracking point
N = ves.N; %shorthand for Number of discratization points
L = ves.L; %shorthand for Length of vesicle
cur = ves.cur; %shorthand for curvature of vesicle
theta = ves.theta; %shorthand for opening angle of vesicle 
op = poten(N); %shorthand for poten class
oc = curve; %shorthand for curve class

%Compute Eu and Esigma, equations (13) and (14) 
[Eu,Esigma] = ves.variationsNonStiff;

%Compute the force in eq (33)
tau = [[+Esigma.*sin(theta) - Eu.*cos(theta)]; ...
       [-Esigma.*cos(theta) - Eu.*sin(theta)]];
%NOTE: missing us term??

%Construct Stokes matrix without the log singularity. ie. A3 only contains
%the rightmost kernel in equation (43)
StokesMat = op.StokesMatrixLogless(ves.X);

%form the velocity, k, on the interface corresponding to v^u in eq (33)
%If [[P^u n]]_sigma = tau, then k = stokesMatrix*tau
k = StokesMat*tau;

%ulam is the viscosity contrast
ulam = ves.viscIn/ves.viscOut;

%LogKernel1 and LogKernel2 are the log kernels in the left term in the
%right hand side of equation (43) integrated against the x and y
%components of the density functions tau that involves Eu and Esigma
LogKernel1 = op.IntegrateLogKernel(tau(1:N));
LogKernel2 = op.IntegrateLogKernel(tau(N+1:end));

%calculate constants that multiply the weakly singular and regular parts of 
%the integral operators LogKernel1 and LogKernel2
c1 = 1/(4*pi)*L/N;
c2 = -L/(8*pi);

%sigma1 and sigma2 are the solution of equation (33) (still unsure
%about the u_s term). Also, it includes the background shear flow which
%is not stated in equation (33), but instead in equations (6) and (7)
sigma1 = k(1:N)*c1 + LogKernel1*c2 + ves.X(N+1:end)*o.shearRate;
sigma2 = k(N+1:end)*c1 + LogKernel2*c2;

% Calculate v dot n in eq (40)
vdotn = sigma1.*sin(theta) - sigma2.*cos(theta);
% Calculate v dot s in eq (40)
vdots = sigma1.*cos(theta) + sigma2.*sin(theta);

%Calculate d/ds(v dot s) in eq (40)
IK = oc.modes(N);
dvdotsds = oc.diffFT(vdots,IK)/L;

%Compute the right hand side of equation (40)
rhs = -(dvdotsds + ves.cur.*vdotn);

%The velocity components in eq (40) are nonlocal linear functions of 
%lambdaTilde. Solve the linear system for LambdaTilde in (39) using GMRES.
%Each iteration of GMRES requires a solution of Stokes equation
%LambdaTilde is the lambda with a tilde in eq (39)
[lambTil,flag,relres,iter,resvec] = ...
      gmres(@(x) o.matvec40(x,StokesMat),rhs,[],o.gmresTol,...
              o.gmresMaxIter); %,@(x) o.preconditioner(x));

%calculate the Fourier derivative of lambdaTilde
dlamTilds = oc.diffFT(lambTil,IK)/L;

%We can now compute the traction jump in first part of equation (39).
%This comes from applying the product rule and using Frenet-Seret.
tracJump = [(-lambTil.*ves.cur.*sin(theta) + dlamTilds.*cos(theta));...
            (+lambTil.*ves.cur.*cos(theta) + dlamTilds.*sin(theta))];

%Adding the jump conditions in eq (39) to (33) which is in the variable tau
tau = -tau - tracJump;

%Compute u tilde in equations (38) through (40) without the weakly singular
%log kernel
k = StokesMat*tau;

%compute the weakly singluar log kernel part
force1 = op.IntegrateLogKernel(tau(1:N));
force2 = op.IntegrateLogKernel(tau(N+1:end));

%Calulating  u in eqatuion (31) by adding the results from the
%non-singular and weakly singular integral operators
uxvel = k(1:N)*c1 + force1*c2 + ves.X(N+1:end)*o.shearRate;
uyvel = k(N+1:end)*c1 + force2*c2;

end % usetself

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LHS = matvec40(o,rhs,StokesMat)   
%This function cooresponds to the matvec in equation 40 and returns
%the LHS of eq(40), (u \cdot s)_s + kappa * (u \cdot n)  as LHS

ves = o.ves; %shorthand for ves object
theta = ves.theta; %shorthand for opening angle
cur = ves.cur; %shorthand for vesicle curvature
L = ves.L; %shorthand for vesicle length
N = ves.N; %shorthand for number of discretization points 
ulam = ves.viscIn/ves.viscOut; %define viscosity ratio
oc = curve; %shorthand for curve class
op = poten(N); %shorthand for poten class
IK = oc.modes(N); %define the Fourier modes for differentiation

%The the derivative of the rhs of eq(40)
drhsds = oc.diffFT(rhs,IK)/L;

%Compute the forces
tau1 = -rhs.*cur.*sin(theta) + drhsds.*cos(theta);
tau2 = +rhs.*cur.*cos(theta) + drhsds.*sin(theta);

%form the velocity on the interface
krhs = StokesMat*[tau1;tau2];

%LogKerneltau1 and LogKerneltau2 are the log kernels integrated against the 
%density functions tau1 and tau2.
LogKerneltau1 = op.IntegrateLogKernel(tau1);
LogKerneltau2 = op.IntegrateLogKernel(tau2);

%calculate constants that multiply the weakly singular and regular parts of 
%the integral operators
c1 = 1/(4*pi)*L/N;
c2 = -L/(8*pi);

%[utilde1 utilde2] is utilde in equations (39) and (40)
utilde1 = krhs(1:N)*c1 + LogKerneltau1*c2;
utilde2 = krhs(N+1:end)*c1 + LogKerneltau2*c2;

% Build left hand side of equation (40)
udotn = utilde1.*sin(theta) - utilde2.*cos(theta);
udots = utilde1.*cos(theta) + utilde2.*sin(theta);
dudotsds = oc.diffFT(udots,IK)/L;

%LHS is (u \cdot s)_s + kappa * (u \cdot n) in eq (40)
LHS = (dudotsds + cur.*udotn);

end %matvec40

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = preconditioner(o,rhs)
% preconditioner (equation (63) and (64)) for equation (40) which is the
% Schur complement for a Lagrange multiplier needed to enforce local
% inextensibility
ves = o.ves;
N = ves.N;
x = fft(rhs,N)./[1 1:1:N/2 N/2-1:-1:1]';
x = real(ifft(x,N));
end %preconditioner

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ves,ux_old,uy_old,L,Ln,dcur0,N1,N2Hat] = FirstSteps(o,ves,params,options,om)
%refines the first time step [0,dt] and uses a first-order Euler method to 
%find the vesicle position, curv, velocity, and density function.  
%Returns ves, L, Ln and dcur0 to use for higher-order time stepping.     
    
time = 0; %current time
oc = curve; %define shorthand for curve class
L = ves.L; %shorthand vesicle length

%smooth the geometry by requiring that the opening tangent angle theta
%is band limited. Note that it will still have some small coefficients
%in the tails of the Fourier spectrum, but they will be much smaller
%than the original theta
ves.smoothGeom;
[~,a_old,l_old] = oc.geomProp(ves.X);

%compute the x velocity, y velocity using the initialized concentration
%field. A step of Cahn-Hilliard is not taken until after this step.
%However, not sure how the terms in equations (13) and (14) have been
%incorporated into the forces. Also missing the u_s term in equation
%(33)
[ux,uy] = o.usetself;

%Plot the position of the vesicle with the velocity vectors at the first 
%time step
%disp('step1')
% figure(1); clf; hold on;
% plot(ves.X(1:end/2),ves.X(end/2+1:end),'r')
% quiver(ves.X(1:end/2),ves.X(end/2+1:end),ux,uy)
% axis equal
% axis([-3 3 -3 3])
% pause(0.1)
% hold off
% pause

%put the x-y velocity into the normal and tangential velocity.
theta = ves.theta; %shorthand theta so we don't type ves. a million times
un = ux.*sin(theta) - uy.*cos(theta); %Normal Velocity
ut = ux.*cos(theta) + uy.*sin(theta); %Tangential Velocity

%Update length change over time using a first-order Euler method.
%For subsequent time steps, will use a multistep method as described in
%equation (60)
%   Forward Euler for the length. The first element of dcur should be zero, 
%   so this is just checking for discretization and round-off errors i.e.
%   Ln = ves.L if params.dt*dcur(1) = 0.
dcur0 = oc.fdcur(ves,un);
Ln = L + params.dt*dcur0(end);

%Get the velocity of the vesicle by its velocity of the tangential
%angle, but without the stiff term. The nonlinear term N_1
%defined in equation (54) but without the alpha derivative multiplied
%by the integral operator \mathcal{L}
N1 = oc.fthetaim(ves,un,ut);

%next evolve the shape and the phase distribution in Fourier space.
%Fourier series of the derivative of the tangent angle. i.e. Fourier
%series of N_1
fsN1 = fft(N1);

%Fourier derivative of the tangent angle adjusted by a linear function
%so that we are taking the fft of a periodic function
fsTA = fft(theta - (2*pi*(0:ves.N-1)/ves.N)');

%Define the Fourier modes, Nk
Nk = 2*pi*[0 1:ves.N/2 ves.N/2-1:-1:1]';

%rsl is the stiff term in equation (55) that will be integrated implicitly 
%using an integrating factor
rsl = ves.bendsti*(Nk/ves.L).^3/4;

%rsln is the next step of rsl 
%'n' is for 'new' since we'll be using multistep
rsln = ves.bendsti*(Nk/Ln).^3/4;

%use trapezoid rule to approximate integral in equation (57). This next
%line is exactly equation (58) for the quadrature
%d1 is now exactly as in equation (57). Each Fourier mode has its own
%integating factor. This is the Forward Euler method that is analagous 
%to equation(56)
ek = exp(-(params.dt*(rsl + rsln)/2));

%thetan is now the tangent angle of the new shape after taking a single
%step of Euler with the stiffest term treated implicitly and integrated
%with an integrating factor. Add back the linear function that makes theta a 
%function that grows by 2*pi everytime you go completely around the shape.
thetan = real(ifft(ek.*(fsTA + params.dt*fsN1)))+2*pi*(0:ves.N-1)'/ves.N;

%        -----  define the lipid species model for u  -----
%Define the Fourier modes, but scaled by 2*pi. Note that these two vectors 
%will be nearly identical since L \approx Ln by inextensibility
%form stiffest term that is treated implicitly, but is also linear (and
%diagonal) in Fourier space
rk = 2*pi*[0 1:ves.N/2 ves.N/2-1:-1:1]'; 

%compute integrating factor components 
rsl = params.epsch*(rk/ves.L).^4*params.consta;
rsln = params.epsch*(rk/Ln).^4*params.consta;

%compute the integrating factor for each Fourier mode using the
%trapezoid rule. d1 is the integrating factor in equation (70)
ek = exp(-(params.dt/params.nloop*(rsl + rsln)/2));

%Take small time steps to move the lipid species from time 0 to time dt
for i=1:params.nloop
  %fncon is the non-linear term N_2 in equation (67)
  N2Hat = oc.frconim(ves,params.epsch,params.consta);
  %fcN2 are the fourier coefficients of N_2 as in equation (68)
  fcN2 = fft(N2Hat);
  %fcLS is the fourier coefficients of the lipid species
  %concentration as in equation (68)
  fcLS = fft(ves.rcon);
  %Compute vesicle concentration using first-order Euler method that is 
  %analagous to equation (69)
  ves.rcon = real(ifft(ek.*(fcLS+params.dt/params.nloop*fcN2)));
end

%              -----  update ves and area -----
%update ves.L
ves.L = Ln;

%update ves.theta
ves.theta = thetan;

%update the position with Forward Euler. At future time steps,
%second-order Adams-Bashforth will be used
ves.x0 = ves.x0 + params.dt*ux(1);
ves.y0 = ves.y0 + params.dt*uy(1);

%Reconstruct ves.X with updated tracking point 
ves.X = oc.recon(ves.N,ves.x0,ves.y0,ves.L,ves.theta);    

%set up variables for timestepping loop
ux_old = ux(1);
uy_old = uy(1);

[~,a_new,l_new] = oc.geomProp(ves.X);
ea = abs(a_new - a_old)./abs(a_old);
el = abs(l_new - l_old)./abs(l_old);
om.plotData(ves.X,time,ea,el,[ux;uy])
om.initializeFiles(ves.X,ves.cur,time,[ux;uy])

end %FirstSteps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ves = TimeStepLoop(o,ves,params,om,ux_old,uy_old,L,Ln,...
                            dcur0,fntheta,N2Hat)

oc = curve;
nstep = round(params.T/params.dt); %total number of time steps
outpt = round(params.outpt/params.dt); %integer values for when output is 

[~,a_old,l_old] = oc.geomProp(ves.X);
%Set up empty array to store norms of x and y velocities
nn = [];

%Entering time stepping loop
for ktime = 1:nstep
  time = ktime*params.dt;
  
  %compute the x- and y-components of the velocity. This is the routine
  %that calls GMRES which is used to solve equation (30) 
  [uxvel_loop, uyvel_loop] = o.usetself;

  %Save the norm of x and y velocities
  nn = [nn;norm([uxvel_loop;uyvel_loop],inf)];
  
  %
  uxvel_new = uxvel_loop(1);
  uyvel_new = uyvel_loop(1);
  
  %put the x-y velocity into the normal and tangential velocity.
  theta = ves.theta; %shorthand for theta
  un = uxvel_loop.*sin(theta) - uyvel_loop.*cos(theta);%Normal Velocity
  ut = uxvel_loop.*cos(theta) + uyvel_loop.*sin(theta);%Tangential Velocity
%   
%   sum(un)
%   pause
  
  %Update length change over time using a 2nd-order Adams Bashforth method.
  %described in equation (60)
  %   2nd-order Adams Bashforth for the length. The first element of dcur 
  %   should be zero, so this is just checking for discretization and 
  %   round-off errors.
  dcur1 = oc.fdcur(ves,un);
  %Lnn = Ln + params.dt*(3*dcur1(end)-dcur0(end))/2;
  Lnn = Ln + params.dt*dcur1(end);
  %update ves.L  
  ves.L = Lnn;
  
  %Get the velocity of the vesicle by its velocity of the tangential angle
  %without the stiff term. 
  fnthetan = oc.fthetaim(ves,un,ut);
  
  %         -----  update the lipid species model for u  -----
  %Define the Fourier modes, but scaled by 2*pi. Note that these two 
  %vectors will be nearly identical since L \approx Ln by inextensibility
  %form stiffest term that is treated implicitly, but is also linear (and
  %diagonal) in Fourier space
  rk = 2*pi*[0 1:ves.N/2 ves.N/2-1:-1:1]'; 
  %compute the nonlocal model for theta
  rsl = ves.bendsti*(rk/L).^3/4;
  rsln = ves.bendsti*(rk/Ln).^3/4; 
  rslnn = ves.bendsti*(rk/Lnn).^3/4;
  %compute the local model for theta using some kind of exponential time 
  %integrators described in equations (56) and (57)	 
  d1 = exp(-(params.dt*(rsln+rslnn)/2));
  d2 = exp(-(params.dt*(rsl+rslnn)/2+params.dt*rsln));
  
  %Compute the fourier series of fntheta
  fcfntheta = fft(fntheta);
  %Compute the Fourier series of fnthetan
  fcfnthetan = fft(fnthetan);
  %Compute the Fourier series of theta adjusted by a linear function
  %so that we are taking the fft of a periodic function
  fcthetan = fft(theta -2*pi*(0:ves.N-1)'/ves.N);
  %add back the linear function that makes theta a function that grows
  %by 2*pi everytime you go completely around the shape. thetann is now the
  %tangent angle of the new shape after taking a single step of a 2nd order
  %multistep method with the stiffest term treated implicitly and 
  %integrated with the integrating factors d1 and d2 defined above.
  %thetann = real(ifft(d1.*fcthetan + 0.5*params.dt*(3*d1.*fcfnthetan- ...
  %          d2.*fcfntheta)))+2*pi*(0:ves.N-1)'/ves.N;
  thetann = real(ifft(d1.*(fcthetan + params.dt*fcfnthetan)))+2*pi*(0:ves.N-1)'/ves.N;
  %compute integrating factor components      
  rsl = params.epsch*(rk/ves.L).^4*params.consta;
  rsln = params.epsch*(rk/Ln).^4*params.consta;
  rslnn = params.epsch*(rk/Lnn).^4*params.consta;
  %compute the integrating factors for each Fourier mode using the
  %trapezoid rule. These are dependent on the number of steps taken to
  %evolve the phase field surface.
  d1 = exp(-(params.dt/params.nloop*(rsln+rslnn)/2));
  d2 = exp(-(params.dt/params.nloop*((rsl+rslnn)/2+rsln)));
  
%                 === Evolve phase field on surface ===
  if params.concentra > 0
    for i =1:params.nloop
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
      rconn = real(ifft(d1.*fcLSn + 0.5*params.dt/params.nloop*(3*d1.*...
              fcN2n - d2.*fcN2)));
      %update ves.rcon
      ves.rcon = rconn;
      %update fcN2 for nloop
      fcN2 = fcN2n;
    end    
  end
%                             ====================
  %update variables for time stepping loop
  dcur0 = dcur1;
  ves.theta = thetann;
  fntheta = fnthetan;   
  
  %update the single tracker point using Adams Bashforth
 %  ves.x0 = ves.x0 + 0.5*params.dt*(3*uxvel_new - ux_old);
 %  ves.y0 = ves.y0 + 0.5*params.dt*(3*uyvel_new - uy_old);
    ves.x0 = ves.x0 + params.dt*uxvel_new;
    ves.y0 = ves.y0 + params.dt*uyvel_new;

  %Update X
  ves.X = oc.recon(ves.N,ves.x0,ves.y0,ves.L,ves.theta);

  %HACK: keep the vesicle centered at (0,0)
  ves.X(1:end/2) = ves.X(1:end/2) - mean(ves.X(1:end/2));
  ves.X(end/2+1:end) = ves.X(end/2+1:end) - mean(ves.X(end/2+1:end));
  
  %Compute the errors in length and area
  [~,a_new,l_new] = oc.geomProp(ves.X);
  ea = abs(a_new - a_old)./abs(a_old);
  el = abs(l_new - l_old)./abs(l_old);
  
  %Print outputs
  om.outputInfo(ves.X,ves.cur,time,[uxvel_loop;uyvel_loop],ea,el)

end

%Update ux_old, and uy_old for timestepping loop
ux_old = uxvel_new;
uy_old = uyvel_new; 

end %TimeStepLoop    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % methods

end % classdef
