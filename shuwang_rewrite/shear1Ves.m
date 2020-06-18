
%%%%%%%%%%%%%%%%% Initialize parameters and options %%%%%%%%%%%%%%%%%%%%%%%
% TODO: SOME OF THESE ARE MORE OPTIONS THAN PARAMETERS
params.N = 64; % points on vesicle
params.dt = 1e-3; % time step size
params.T = 0.5; % time horizon
params.outpt = 1e-3; % ouptut frequency
params.concentra = 0.5; % constant, initial concentration of lipid species
params.oddeven = 0; % flag for initial lipid species profile?
params.shortax = 0.5; % short axis length
params.shearRate = 30; % shear rate
params.viscosityInside = 1.0;
params.viscosityOutside = 2.0;
params.bendsti = 1; % maximum bending stiffness
params.bendratio = 0.1; % ratio between max and min bending stiffness
params.consta = 100; % parameter 'a' in the Cahn-Hilliard energy
params.nloop = 20; 
% number of time steps of Cahn-Hilliard to be taken at each time step of
% the hydrodynamics
params.epsch = 5e-2; % small parameter  in the double-well potential 
params.gmresTol = 1e-10; %GMRES tolerance
params.gmresMaxIter = params.N; %maximum number of GMRES iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time = 0; %current time
nstep = round(params.T/params.dt); %total number of time steps
outpt = round(params.outpt/params.dt); %integer values for when output is 
                                       %saved
oc = curve;
%Form initial shape that has points equally spaced in arclength. alpha
%has the parameter values that give this parameterization
alpha = (0:params.N-1)'/params.N;
[alpha,X] = oc.initConfig(params.N,false,'ellipse',...
            'shortax',params.shortax,'parameter',alpha);
plot(X(1:end/2),X(end/2 +1:end))
pause(1)
%Define the initial concentration field
rcon = oc.initialConcentration(params.N,alpha,...
      params.concentra,params.oddeven);
%build object for the vesicle but without a band-limited opening angle
ves = capsules(X,rcon,params);
plot(ves.X(1:end/2),ves.X(end/2 +1:end))
pause(1)
%smooth the geometry by requiring that the opening tangent angle theta
%is band limited. Note that it will still have some small coefficients
%in the tails of the Fourier spectrum, but they will be much smaller
%than the original theta
ves.smoothGeom;
%Initialize time-stepping class
tt = tstep(params,ves);
%compute the x velocity, y velocity using the initialized concentration
%field. A step of Cahn-Hilliard is not taken until after this step.
%However, not sure how the terms in equations (13) and (14) have been
%incorporated into the forces. Also missing the u_s term in equation
%(33)
[un, ut] = tt.usetself;
%put the x-y velocity into the normal and tangential velocity.
theta = ves.theta;
unt  = un.*sin(theta) - ut.*cos(theta); %Tangential Velocity
utn = un.*cos(theta) + ut.*sin(theta); %Normal Velocity
%Update length change over time using a first-order Euler method.
%For subsequent time steps, will use a multistep method as described in
%equation (60)
%   Forward Euler for the length. The first element of dcur should be zero, 
%   so this is just checking for discretization and round-off errors i.e.
%   Ln = ves.L if params.dt*dcur(1) = 0.
dcur0 = oc.fdcur(ves,unt);
Ln0 = ves.L + params.dt*dcur0(1);
ves.L = Ln0;
%Get the velocity of the vesicle by its velocity of the tangential
%angle, but without the stiff term. fntheta is the nonlinear term N_1
%defined in equation (54) but without the alpha derivative multiplied
%by the integral operator \mathcal{L}
fntheta = oc.fthetaim(ves,unt,utn);
%next evolve the shape and the phase distribution in Fourier space.
%Fourier series of the derivative of the tangent angle. i.e. Fourier
%series of N_1
fsN1 = fft(theta);
%Fourier derivative of the tangent angle adjusted by a linear function
%so that we are taking the fft of a periodic function
fsTA = fft(theta - [2*pi*(0:ves.N-1)/ves.N]');
%Nk are the Fourier modes
Nk = 2*pi*[0 1:ves.N/2 ves.N/2-1:-1:1]';
%stiff term in equation (55) that will be integrated implicitly using
%an integrating factor
rsl = ves.bendsti*(Nk/ves.L).^3/4;
% 'n' is for 'new' since we'll be using multistep
rsln = ves.bendsti*(Nk/ves.L).^3/4;
%use trapezoid rule to approximate integral in equation (57). This next
%line is exactly equation (58) for the quadrature
%d1 is now exactly as in equation (57). Each Fourier mode has its own
%integating factor
d1 = exp(-(params.dt*(rsl + rsln)/2));
%This is the Forward Euler method that is analagous to equation(56) and
%add back the linear function that makes theta a function that grows
%by 2*pi everytime you go completely around the shape (ex. vesicle).
%thetan is now the tangent angle of the new shape after taking a single
%step of Euler with the stiffest term treated implicitly and integrated
%with an integrating factor
thetan = real(ifft(d1.*(fsTA + params.dt*fsN1)))+2*pi*(0:ves.N-1)'/ves.N;
ves.theta = thetan;
%lipid species model for u
rk = 2*pi*[0 1:ves.N/2 ves.N/2-1:-1:1]'; 
%Fourier modes, but scaled by 2*pi. Note that these two vectors will be
%nearly identical since L \approx Ln by inextensibility
%form stiffest term that is treated implicitly, but is also linear (and
%diagonal) in Fourier space
rsl = params.epsch*(rk/ves.L).^4*params.consta;
rsln = params.epsch*(rk/ves.L).^4*params.consta;
%compute the integrating factor for each Fourier mode using the
%trapezoid rule. d1 is the integrating factor in equation (70)
d1 = exp(-(params.dt/params.nloop*(rsl + rsln)/2));

%Take small time steps to move the lipid species from time 0 to time dt
for i=1:params.nloop
  %fncon is the non-linear term N_2 in equation (67)
  N2Hat = oc.frconim(ves,params.epsch,params.consta);
  %fcN2 are the fourier coefficients of N_2 as in equation (68)
  fcN2 = fft(N2Hat);
  %fcLS is the fourier coefficients of the lipid species
  %concentration as in equation (68)
  fcLS = fft(rcon);
  %First-order Euler method that is analagous to equation (69)
  ves.rcon = real(ifft(d1.*(fcLS+params.dt/params.nloop*fcN2)));
end

%New area
area = sum(sin(ves.theta).*X(1:end/2) - cos(ves.theta).*X(end/2+1:end))*...
       0.5*ves.L/ves.N;
%update the position with Forward Euler. At future time steps,
%second-order Adams-Bashforth will be used
ves.x0 = ves.x0 + params.dt*un(1);
ves.y0 = ves.y0 + params.dt*ut(1);
X = oc.recon(ves.N,ves.x0,ves.y0,ves.L,ves.theta);
plot(real(X(1:end/2)),real(X(end/2 +1:end)))
pause(1)

% From the second step, use multistep the evolve the dynamics.
for ktime = 1:nstep
  tcomp = ktime*params.dt;
  un0 = un(1);
  ut0 = ut(1);
  
  %compute the normal and tangential velocity. This is the routine that 
  %calls GMRES which is used to solve equation (30) in the Sohn et al JCP
  %paper (2010)
  [unloop, utloop] = tt.usetself;
  un1 = unloop(1);
  ut1 = utloop(1);
  %put the x-y velocity into the normal and tangential velocity.
  theta = ves.theta;
  unt  = unloop.*sin(theta) - utloop.*cos(theta); %Tangential Velocity
  utn = unloop.*cos(theta) + utloop.*sin(theta); %Normal Velocity
  %Update length change over time using a first-order Euler method.
  %For subsequent time steps, will use a multistep method as described in
  %equation (60)
  %   Forward Euler for the length. The first element of dcur should be zero, 
  %   so this is just checking for discretization and round-off errors i.e.
  %   Ln = ves.L if params.dt*dcur(1) = 0.
  dcur1 = oc.fdcur(ves,unt);
  Ln1 = Ln0 + params.dt*(3*dcur1(1)-dcur0(1))/2;
  ves.L = Ln1;
  %Get the velocity of the vesicle by its velocity of the tangential
  %angle, but without the stiff term. 
  fnthetan = oc.fthetaim(ves,unt,utn);
  %lipid species model for u
  %nonlocal model for theta
  rsl = ves.bendsti*(rk/ves.L).^3/4;
  rsln = ves.bendsti*(rk/Ln0).^3/4; 
  rslnn = ves.bendsti*(rk/Ln1).^3/4;
  %Local model	 
  d1 = exp(-(params.dt*(rsln+rslnn)/2));
  d2 = exp(-(params.dt*(rsl+rslnn)/2+params.dt*rsln));
  %using some kind of exponential time integrator (see equations (56)
  %and (57))
  %fourier series of fntheta
  fcfntheta = fft(fntheta);
  %Fourier series of fnthetan
  fcfnthetan = fft(fnthetan);
  %Fourier series of thetan adjusted by a linear function
  %so that we are taking the fft of a periodic function
  fcthetan = fft(thetan -2*pi*(0:ves.N-1)'/ves.N);
  %add back the linear function that makes thetan a function that grows
  %by 2*pi everytime you go completely around the shape (ex. vesicle).
  %thetann is now the tangent angle of the new shape after taking a single
  %step of Euler with the stiffest term treated implicitly and integrated
  %with an integrating factor
  ves.theta = real(ifft(d1.*fcthetan + 0.5*params.dt*(3*d1.*fcfnthetan- ...
            d2.*fcfntheta)))+2*pi*(0:ves.N-1)'/ves.N;
  %Do we need to keep this? there's no change between sl, sln, slnn
  
  rsl = params.epsch*(rk/ves.L).^4*params.consta;
  rsln = params.epsch*(rk/Ln0).^4*params.consta;
  rslnn = params.epsch*(rk/Ln1).^4*params.consta;
  
  d1 = exp(-(params.dt/params.nloop*(rsln+rslnn)/2));
  d2 = exp(-(params.dt/params.nloop*((rsl+rslnn)/2)+rsln));
  
  if params.concentra > 0
    for i =1:params.nloop
      %evolve the phase field on the surface            
      N2Hatn = oc.frconim(ves,params.epsch,params.consta);
      %fcN2 are the fourier coefficients of N2
      fcN2 = fft(N2Hat);
      %fcN2n are the fourier coefficients of N2n
      fcN2n = fft(N2Hatn);
      %fcLS is the fourier coefficients of the lipid species
      %concentration as in equation (68)
      fcLSn = fft(ves.rcon);  
      rconn = real(ifft(d1.*fcLSn + 0.5*params.dt/params.nloop*(3*d1.*...
              fcN2n - d2.*fcN2)));
      ves.rcon = rconn;
      fcN2 = fcN2n;
    end    
  end
  %update variables for loop
  dcur0 = dcur1;
  fntheta = fnthetan;    
  %update the single tracker point
  ves.x0 = ves.x0 + 0.5*params.dt*(3*un1 - un0);
  ves.y0 = ves.y0 + 0.5*params.dt*(3*ut1 - ut0);
  %Update X
  X = oc.recon(ves.N,ves.x0,ves.y0,ves.L,ves.theta);
  plot(real(X(1:end/2)),real(X(end/2 +1:end)))
  pause(1)
end
