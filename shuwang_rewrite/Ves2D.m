function ves = Ves2D(params,options)
% Ves2D does time stepping on the intial configuration X
% with parameters and options defined in prams and options.

oc = curve;

% Form initial shape that has points equally spaced in arclength. alpha
% has the parameter values that give this parameterization
upRate = 4;
h = (params.shortax - 1)^2/(params.shortax+1)^2;
% Lscale = 0.4626;
% L = Lscale*pi*(params.shortax + 1)*(1+h/4+h^2/64);
L = pi*(params.shortax + 1)*(1+h/4+h^2/64);
scale = 1;%0.6516;
angle = pi/2;
[alpha,X] = oc.initConfig(upRate*params.N,false,'ellipse',...
            'shortax',params.shortax, 'scale', scale, 'angle', angle);
% figure(3);clf;
% plot(alpha)
% pause
%Define the initial concentration field
rcon = oc.initialConcentration(upRate*params.N,alpha,...
      params.concentra,params.oddeven);
% figure(3)
% plot(rcon)
% pause
%build object for the vesicle but without a band-limited opening angle
ves = capsules(X,rcon,params);

% Downsample back to the desired resolution
ves.N = ves.N/upRate;
ves.X = ves.X(1:upRate:end);
ves.cur = ves.cur(1:upRate:end);
ves.theta = ves.theta(1:upRate:end);
ves.rcon = ves.rcon(1:upRate:end);
% semilogy(abs(fftshift(fft(ves.X(1:end/2)))))
% pause
% Smooth the geometry by requiring that the opening tangent angle theta
% is band limited. Note that it will still have some small coefficients
% in the tails of the Fourier spectrum, but they will be much smaller
% than the original theta
%ves.smoothGeom;
tt = tstep(params,ves); % set up tstep class
om = monitor(ves.X,params,options); % set up monitor class


% Take first step with first-order Euler
[ves,ux_old,uy_old,L,Ln,dcur0,fntheta,N2Hat,cx0,cy0] = ...
      tt.FirstSteps(ves,params,options,om);

% Begin time step routine
ves = tt.TimeStepLoop(ves,params,om,ux_old,uy_old,L,Ln,dcur0,...
                      fntheta,N2Hat,cx0,cy0);


end % Ves2D
