
%%%%%%%%%%%%%%%%% Initialize parameters and options %%%%%%%%%%%%%%%%%%%%%%%
% TODO: SOME OF THESE ARE MORE OPTIONS THAN PARAMETERS
params.N = 2*96; % points on vesicle
params.dt = 1e-3*10; % time step size
params.T = 10; % time horizon
params.outpt = 1e-3; % ouptut frequency
params.concentra = 0.0; % constant, initial concentration of lipid species
params.oddeven = 0; % flag for initial lipid species profile?
params.shortax = 2.0; % short axis length
params.shearRate = 1; % shear rate
params.viscosityInside = 1.0;
params.viscosityOutside = 1.0;
params.bendsti = 1; % maximum bending stiffness
params.bendratio = 1; % ratio between max and min bending stiffness
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
%Define the initial concentration field
rcon = oc.initialConcentration(params.N,alpha,...
      params.concentra,params.oddeven);
%build object for the vesicle but without a band-limited opening angle
ves = capsules(X,rcon,params);
%clf;
% plot([ves.X(1:end/2);ves.X(1)],[ves.X(end/2 +1:end);ves.X(end/2+1)],'b-o')
% axis equal
% % hold on
% pause(0.1)
%semilogy(abs(fftshift(fft(ves.theta - 2*pi*alpha))))
%hold on

%smooth the geometry by requiring that the opening tangent angle theta
%is band limited. Note that it will still have some small coefficients
%in the tails of the Fourier spectrum, but they will be much smaller
%than the original theta
%clf; hold on
%plot(ves.X(1:end/2),ves.X(end/2+1:end))
ves.smoothGeom;
%plot(ves.X(1:end/2),ves.X(end/2+1:end),'r--')
%axis equal
%pause

%plot([ves.X(1:end/2);ves.X(1)],[ves.X(end/2 +1:end);ves.X(end/2+1)],'r--')
%axis equal
%semilogy(abs(fftshift(fft(ves.theta - 2*pi*alpha))),'ro')
%pause()

%plot(ves.theta)
%pause
%Initialize time-stepping class
tt = tstep(params,ves);

%compute the x velocity, y velocity using the initialized concentration
%field. A step of Cahn-Hilliard is not taken until after this step.
%However, not sure how the terms in equations (13) and (14) have been
%incorporated into the forces. Also missing the u_s term in equation
%(33)
[ux,uy] = tt.usetself;
%  disp('plotting ux and uy')
%  plot(ux)
%  hold on
%  plot(uy)
%  pause
% clf
% plot(ves.X(1:end/2),ves.X(end/2+1:end))
% pause
%plot([ves.X(1:end/2);ves.X(1)],[ves.X(end/2 +1:end);ves.X(end/2+1)],'b')
disp('step1')
figure(1); clf; hold on;
plot(ves.X(1:end/2),ves.X(end/2+1:end),'r')
quiver(ves.X(1:end/2),ves.X(end/2+1:end),ux,uy)
axis equal
axis([-3 3 -3 3])
pause(0.1)
hold off
% pause
% NOTE: THERE IS A BUG IN COMPUTING THE VELOCITY BEFORE THIS POINT

%put the x-y velocity into the normal and tangential velocity.
theta = ves.theta;
ut = ux.*sin(theta) - uy.*cos(theta); %Tangential Velocity
un = ux.*cos(theta) + uy.*sin(theta); %Normal Velocity

%Update length change over time using a first-order Euler method.
%For subsequent time steps, will use a multistep method as described in
%equation (60)
%   Forward Euler for the length. The first element of dcur should be zero, 
%   so this is just checking for discretization and round-off errors i.e.
%   Ln = ves.L if params.dt*dcur(1) = 0.
dcur0 = oc.fdcur(ves,ut);
Ln0 = ves.L + params.dt*dcur0(end);
ves.L = Ln0;
%Get the velocity of the vesicle by its velocity of the tangential
%angle, but without the stiff term. fntheta is the nonlinear term N_1
%defined in equation (54) but without the alpha derivative multiplied
%by the integral operator \mathcal{L}
fntheta = oc.fthetaim(ves,ut,un);

%next evolve the shape and the phase distribution in Fourier space.
%Fourier series of the derivative of the tangent angle. i.e. Fourier
%series of N_1
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

fsN1 = fft(fntheta);
% clf
% plot(fsN1)
% pause
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
%  disp('plotting rcon')
%  plot(ves.rcon)
%  pause%(0.1)
%New area
area = sum(sin(ves.theta).*ves.X(1:end/2) - ...
    cos(ves.theta).*ves.X(end/2+1:end))*0.5*ves.L/ves.N;
%update the position with Forward Euler. At future time steps,
%second-order Adams-Bashforth will be used
ves.x0 = ves.x0 + params.dt*ux(1);
ves.y0 = ves.y0 + params.dt*uy(1);
%============================================================
ves.X = oc.recon(ves.N,ves.x0,ves.y0,ves.L,ves.theta);
% plot([ves.X(1:end/2);ves.X(1)],[ves.X(end/2 +1:end);ves.X(end/2 +1)])
% axis equal
% pause

% From the second step, use multistep the evolve the dynamics.
ux_old = ux(1);
uy_old = uy(1);
nn = [];
for ktime = 1:nstep
  tcomp = ktime*params.dt;  
  %compute the x- and y-components of the velocity. This is the routine
  %that calls GMRES which is used to solve equation (30) in the Sohn et
  %al JCP paper (2010)
  [uxloop, uyloop] = tt.usetself;
%   clf;disp('step 2')
%   plot(uxloop)
%   hold on
%   %plot(uyloop)
%   pause
  disp('step')
  disp(ktime+1)
  figure(1); clf; hold on;
  plot(ves.X(1:end/2),ves.X(end/2+1:end),'r')
  quiver(ves.X(1:end/2),ves.X(end/2+1:end),uxloop,uyloop)
  nn = [nn;norm([uxloop;uyloop],inf)];
  disp(norm([uxloop;uyloop],inf))
  axis equal
  %axis([-3 3 -3 3])
  pause(0.1)
  hold off
  
  ux_new = uxloop(1);
  uy_new = uyloop(1);
  %put the x-y velocity into the normal and tangential velocity.
  ves.theta = thetan;
  theta = ves.theta;
  ut = uxloop.*sin(theta) - uyloop.*cos(theta); %Tangential Velocity
  un = uxloop.*cos(theta) + uyloop.*sin(theta); %Normal Velocity
  %Update length change over time using a 2nd-order Adams Bashforth method.
  %described in equation (60)
  %   2nd-order Adams Bashforth for the length. The first element of dcur should be zero, 
  %   so this is just checking for discretization and round-off errors i.e.
  %   Ln = ves.L if params.dt*dcur(1) = 0.
  dcur1 = oc.fdcur(ves,ut);
  Ln1 = Ln0 + params.dt*(3*dcur1(end)-dcur0(end))/2;
%  Multistep
  ves.L = Ln1;
  %Get the velocity of the vesicle by its velocity of the tangential
  %angle, but without the stiff term. 
  fnthetan = oc.fthetaim(ves,ut,un);
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
  %step of a 2nd order multistep method with the stiffest term treated
  %implicitly and integratedwith an integrating factor
  thetann = real(ifft(d1.*fcthetan + 0.5*params.dt*(3*d1.*fcfnthetan- ...
            d2.*fcfntheta)))+2*pi*(0:ves.N-1)'/ves.N;
  % Multistep

  %thetann = real(ifft(d1.*fcthetan + params.dt*d1.*fcfnthetan)) + ...
  %      2*pi*(0:ves.N-1)'/ves.N;
  
  rsl = params.epsch*(rk/ves.L).^4*params.consta;
  rsln = params.epsch*(rk/Ln0).^4*params.consta;
  rslnn = params.epsch*(rk/Ln1).^4*params.consta;
  
  d1 = exp(-(params.dt/params.nloop*(rsln+rslnn)/2));
  d2 = exp(-(params.dt/params.nloop*((rsl+rslnn)/2+rsln)));

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
  %pause
  %update variables for loop
  dcur0 = dcur1;
  thetan = thetann;
  fntheta = fnthetan;    
  %update the single tracker point
  ves.x0 = ves.x0 + 0.5*params.dt*(3*ux_new - ux_old);
  ves.y0 = ves.y0 + 0.5*params.dt*(3*uy_new - uy_old);
  % multistep
%  ves.x0 = ves.x0 + params.dt*ux1;
%  ves.y0 = ves.y0 + params.dt*uy1;
  % Euler
  %Update X
  ves.X = oc.recon(ves.N,ves.x0,ves.y0,ves.L,ves.theta);
  %plot(real(X(1:end/2)),real(X(end/2 +1:end)))
  %pause(1)
  ux_old = ux_new;
  uy_old = uy_new;
end
