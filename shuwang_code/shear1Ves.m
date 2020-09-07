%setup the intial data for the coefficient
shearate = 0;
shape = 3.0;
phi = 0;
scale = 1;

global consta eps_ch kmatrix velocity bendsti bendratio uinside uoutside m

initialdata = [96 1e-2 1/scale 1e-3 phi 0 shape  ...
    shearate 1 1 1 100 1.0 1.0];

ngrid = initialdata(1);
% the number of the grid points
dt = initialdata(2);
% time steps size
T = initialdata(3)*1;
% time horizon
outpt = initialdata(4);
% output data groups
concentra = initialdata(5);
% concentration of the phase
oddeven = initialdata(6);
% odd or even initial phase distribution
% oddeven symmetry: random when -1.  nonsymmetric but set by sin & cos
% when 0,  cos symmetry when 1,  sin symmetry when 2.

shortax = initialdata(7);
% shape parameter (short axis)
velocity = initialdata(8);
% velocity
bendsti = initialdata(9);
% bending stiffness
bendratio = initialdata(10);
% bending stiffness ratio
if concentra==0
  bendratio=1;
end
lambd = initialdata(11);
% viscosity ratio
consta = initialdata(12);
% the speed constant for phase decomposition
uinside = initialdata(13);
% viscosity inside
uoutside = initialdata(14);
% viscosity outside
eps_ch = 5e-2;
% constant for the double well potential of the phase distribution
nloop = 20;
 
ktime = 0;
m = ngrid;
% number of points on vesicle

tcomp=dt*ktime;
nstep = round(T/dt);
% number of time steps
outpt = round(outpt/dt);
% steps for saving to output

% set the initial condition.
[x,y,theta,rcon,sl] = initialsetup(shortax,ngrid,concentra,oddeven);
% here we keep the total arclength unchanged.
x0 = x(1);
y0 = y(1);
%%sl = 2.6442;
%%sl = 2*pi;
%sl = 4.844224110050042;
% For Shuwang: Why is sl defined here? It is actually an output
% of the call to initialsetup. It is just currently not asked
% for in the call above. NOTE: We requested sl to be an output in
% initialsetup
%x0 = x0/2;
%y0 = y0/2;
kmatrix = formkmatrix(ngrid);
% kmatrix is checkerboard pattern of 0s and 1s which is dot multiplied y
% the layer potential matricies so that the quadrature is odd/even
% trapezoid rule.

% compute the x velocity, y velocity, Lagrange multiplier, body forces,
% new shape, and center of mass using the initialized concentration
% field. A step of Cahn-Hilliard is not taken until after this step.
% However, not sure how the terms in equations (13) and (14) have been
% incorporated into the forces. Also missing the u_s term in equation
% (33)
[ux0,uy0,rlambdalnew,x,y,forc1,forc2,xcc,ycc] = ...
    usetself(x0,y0,sl,theta,rcon);
%  disp('plotting un and ut')
%  plot(ux0)
%  hold on
%  plot(uy0)
%  pause
% clf
% plot(x,y)
% pause
disp('step1')
figure(1); clf; hold on;
quiver(x(1:end-1),y(1:end-1),ux0 - 0*y(1:end-1),uy0)
axis equal
axis([-3 3 -3 3])
pause(0.01)
hold off
%pause
u1x = ux0(1);
u1y = uy0(1);
% put the x-y velocity into the normal and tangential velocity.
un  = ux0.*sin(theta) - uy0.*cos(theta); % ???Tangential Velocity???
utt = ux0.*cos(theta) + uy0.*sin(theta); % ???Normal Velocity???

% Update arc length change over time using a first-order Euler method.
% For subsequent time steps, will use a multistep method as described in
% equation (60)
fsl = forcsl(m,theta,un);
sln = sl + dt*fsl;
% Forward Euler for the arclength. fsl should be zero, so this is just
% checking for discretization and round-off errors

% Get the velocity of the vesicle by its velocity of the tangential
% angle, but without the stiff term. fntheta is the nonlinear term N_1
% defined in equation (54) but without the alpha derivative multiplied
% by the integral operator \mathcal{L}
fntheta = fthetaim(m,sl,theta,bendsti,un,utt);

% compute the non-stiff term for the velocity of the concentration
% (rcon) of the lipid species. This is what they call N_2.
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

fncon = frconim(m,sl,rcon,theta,bendsti,bendratio,eps_ch,consta);

% next evolve the shape and the phase distribution in Fourier space.
% Fourier series of the derivative of the tangent angle. i.e. Fourier
% series of N_1
temp1 = fft(fntheta); 
% clf
% plot(temp1)
% pause
% Fourier derivative of the tangent angle adjusted by a linear function
% so that we are taking the fft of a periodic function
temp2 = fft(theta - 2*pi*(0:m-1)/m);

% Nk are the Fourier modes
Nk = pi*2*[0 1:m/2 m/2-1:-1:1];
% stiff term in equation (55) that will be integrated implicitly using
% an integrating factor
rsl = bendsti*(Nk/sl).^3/4;
% 'n' is for 'new' since we'll be using multistep
rsln = bendsti*(Nk/sln).^3/4;
% use trapezoid rule to approximate integral in equation (57). This next
% line is exactly equation (58) for the quadrature
d1 = dt*(rsl + rsln)/2;
% d1 is now exactly as in equation (57). Each Fourier mode has its own
% integating factor
d1 = exp(-d1);
% This is the Forward Euler method that is analagous to equation(56)
temp4(1,1:m) = d1.*(temp2(1,1:m) + dt*temp1(1,1:m));
temp4 = real(ifft(temp4));

% add back on the linear function that makes theta a function that grows
% by 2*pi everytime you go completely around the shape (ex. vesicle).
% thetan is now the tangent angle of the new shape after taking a single
% step of Euler with the stiffest term treated implicitly and integrated
% with an integrating factor
thetan = temp4 + 2*pi*(0:m-1)/m;

% lipid species model for u
rk = 2*pi*[0 1:m/2 m/2-1:-1:1]; 
% Fourier modes scaled by 2*pi. Note that these two vectors will be
% nearly identical since sl is approximately sln by inextensibility.
rsl = eps_ch*(rk/sl).^4*consta;
rsln = eps_ch*(rk/sln).^4*consta;
% form stiffest term that is treated implicitly, but is also linear (and
% diagonal) in Fourier space

% compute the integrating factor for each Fourier mode using the
% trapezoid rule
d1 = dt/nloop*(rsl + rsln)/2;
% d1 is the integrating factor in equation (70)
d1 = exp(-d1);
% Take small time steps to move the lipid species from time 0 to time dt
for i=1:nloop
  % fncon is the non-linear term N_2 in equation (67)
  fncon = frconim(m,sl,rcon,theta,bendsti,bendratio,eps_ch,consta);
  % temp1c is the fourier coefficients of N_2 as in equation (68)
  temp1c = fft(fncon(1:m));
  % temp2c is the fourier coefficients of the lipid species
  % concentration as in equation (68)
  temp2c = fft(rcon(1:m));
  % First-order Euler method that is analagous to equation (69)
  temp4 = d1.*(temp2c + dt/nloop*temp1c);
  rconn = real(ifft(temp4));
  rcon = rconn;
end
%  disp('plotting rcon')
%  plot(rcon)
%  pause
% For Shuwang: Why do we have to take a time step size that is 1/20 the
% size of the time step size used for the hydrodynamics?

areasum = sum(sin(theta).*x(1:m)-cos(theta).*y(1,1:m))/2*sl/m;
% update the position with Forward Euler. At future time steps,
% second-order Adams-Bashforth will be used
x0 = x0 + dt*u1x;
y0 = y0 + dt*u1y;

% write results to variables that will store results at some time steps
% depending on parameter settings
tcomp = dt*ktime;
np = ktime/outpt+1;
lambda = rlambdalnew(1:m)';
xx = x(1:m)';
yy = y(1:m)'; 
AA(1) = areasum;
ss(1) = sl;
vx = ux0';
vy = uy0';

nn = [];
% From the second step, use multistep the evolve the dynamics.
for ktime = 1:nstep
  tic
  tcomp = ktime*dt;
  u1x0 = u1x;
  u1y0 = u1y;

  % compute the normal and tangential velocity. Inputs are the single
  % tracker point (x0,y0), total length sl, tangent angle thetan, and
  % concentration of lipid species rconn. This is the routine that calls
  % GMRES which is used to solve equation (30) in the Sohn et al JCP
  % paper (2010)
  % ux0 = x-velocity
  % uy0 = y-velocity
  % (x,y) = tracker point
  [ux0,uy0,rlambdalnew,x,y,forc1,forc2,xcc,ycc,area] = ...
      usetself(x0,y0,sl,thetan,rconn);
%   disp('here1')
%   norm(uy0)
%   pause
  disp('step')
  disp(ktime+1)
  figure(1); clf; hold on;
  plot(x,y,'r')
  quiver(x(1:end-1),y(1:end-1),ux0,uy0)
  nn = [nn;norm([ux0;uy0])];
  norm([ux0;uy0])
  axis equal
  axis([-3 3 -3 3])
  pause(0.1)
  hold off
  
  u1x=ux0(1);
  u1y=uy0(1);    

  un  = ux0.*sin(thetan) - uy0.*cos(thetan);
  utt = ux0.*cos(thetan) + uy0.*sin(thetan);
  
  fsln = forcsl(m,thetan,un);
  slnn = sln + dt*(3*fsln-fsl)/2;
  fnthetan = fthetaim(m,sln,thetan,bendsti,un,utt);
  % one 'n' after variable means at time n-1?
  % two 'n' after variable means at time n?

  temp1 = fft(fnthetan);
  temp2 = fft(thetan - 2*pi*(0:1:m-1)/m);
  temp3 = fft(fntheta);

  % nonlocal model for theta
  rsl = bendsti*(rk/sl).^3/4;
  rslnn = bendsti*(rk/slnn).^3/4;
  rsln = bendsti*(rk/sln).^3/4; 
  % Local model	 
  d1 = dt*(rsln+rslnn)/2;
  d2 = dt*(rsl+rslnn)/2+dt*rsln;
  d1 = exp(-d1);
  d2 = exp(-d2);
  % using some kind of exponential time integrator (see equations (56)
  % and (57))

  temp4 = d1.*temp2 + 1/2*dt*(3*d1.*temp1-d2.*temp3);
  temp4 = real(ifft(temp4));
  thetann = temp4 + 2*pi*(0:1:m-1)/m;
%   disp('here')
%   norm(fnthetan)
%   pause
  % this is what might need to change to do semi-permeability???

  rsl = eps_ch*(rk/sl).^4*consta;
  rsln = eps_ch*(rk/sln).^4*consta;
  rslnn = eps_ch*(rk/slnn).^4*consta;

  d1 = dt/nloop*(rsln+rslnn)/2;
  d2 = dt/nloop*(rsl+rslnn)/2+dt/nloop*rsln;
  d1 = exp(-d1);
  d2 = exp(-d2);
  
  if concentra>0
    disp('in here')
    for innerstep=1:nloop
      % evolve the phase field on the surface            
      fnncon=frconim(m,sl,rconn,thetan,bendsti,bendratio,...
          eps_ch,consta);        

      temp1c = fft(fnncon);
      temp2c = fft(rconn);
      temp3c = fft(fncon);       

      temp5 = d1.*temp2c + 0.5*(dt/nloop)*...
          (3*d1.*temp1c-d2.*temp3c);        
      temp5 = real(ifft(temp5));

      rconnn = temp5;
      rcon = rconn;
      rconn = rconnn;
      fncon = fnncon;
    end
  end
  %pause
  sl = sln;
  sln = slnn;
  fsl = fsln;
  theta = thetan;
  thetan = thetann;
  fntheta = fnthetan;
     
  x0 = x0 + 1/2*dt*(3*u1x - u1x0);
  y0 = y0 + 1/2*dt*(3*u1y - u1y0); 
  % update the single tracker point
  
  theta = real(theta);
  thetan = real(thetan);    
  rconn = real(rconn);
  rcon = real(rcon);
  x0 = real(x0);
  y0 = real(y0);
%  write results if the time is right
  [x,y] = recon(m,x0,y0,sl,theta);

%  clf;
%  plot(theta)
%  pause

  if mod(ktime,outpt)==0
    np = np+1;
    concen(1:m,np) = rconn';
    xx(1:m,np) = x(1,1:m)';
    yy(1:m,np) = y(1,1:m)'; 
    tt(np) = tcomp; 
    lambda(1:m,np) = rlambdalnew';
    the(1:m,np) = theta(1,1:m)';
    arclength(np) = abs(sl);
    AA(np) = area;
    vx(1:m,np) = ux0(1,1:m)';
    vy(1:m,np) = uy0(1,1:m)';
    ftx(1:m,np) = forc1(1,1:m)';
    fty(1:m,np) = forc2(1,1:m)';        
    inclinang(np) = inclinationAngle([x(1:m)';y(1:m)']);
  end

  toc
end

save(strcat('Results_phi',num2str(round(concentra*100)),...
    'shape',num2str(round(shortax*100)),'shear',...
    num2str(velocity),'.mat'),...
    'xx','yy','tt','lambda','concen','the','vx','vy',...
    'arclength','AA','inclinang','ftx','fty')

