function ves = Ves2D(params,options)
% Ves2D does time stepping on the intial configuration X
% with parameters and options defined in prams and options.


% Form initial shape that has points equally spaced in arclength. alpha
% has the parameter values that give this parameterization
upRate = 4;
oc = curve(upRate*params.N);
h = (params.shortax - 1)^2/(params.shortax+1)^2;
scaleL = params.scaleL;
L = pi*(params.shortax + 1)*(1+h/4+h^2/64)*scaleL;

angle =  params.angle;
vesCenter = params.vesCenter;
geomCenter = params.geomCenter;
wallGeometry = params.wallGeometry;
vesGeometry = params.vesGeometry;
% [alpha,X] = oc.initConfig(upRate*params.N,false,...
%             'shortax',params.shortax, 'scale', scaleL, 'angle', angle, ...
%             'center', center, 'FFflow', FFflow);
[alpha,X] = oc.initConfig(upRate*params.N,false,...
            'shortax',params.shortax, 'scale', scaleL, 'angle', angle, ...
            'center', vesCenter, 'geometry', vesGeometry);
% [RA,~,L] = oc.geomProp(X)
% pause
% Define the initial concentration field
rcon = oc.initialConcentration(upRate*params.N,alpha,...
      params.concentra,params.oddeven);

% Build object for the vesicle but without a band-limited opening angle
ves = capsules(X,rcon,params);
% Downsample back to the desired resolution
ves.N = ves.N/upRate;
ves.X = ves.X(1:upRate:end);
ves.cur = ves.cur(1:upRate:end);
ves.theta = ves.theta(1:upRate:end);
ves.rcon = ves.rcon(1:upRate:end);
% Smooth the geometry by requiring that the opening tangent angle, theta,
% is band limited. Note that it will still have some small coefficients
% in the tails of the Fourier spectrum, but they will be much smaller
% than the original theta
ves.smoothGeom;

oc = curve(ves.N);

% Reconstruct the vesicle position using the new, smooth theta
ves.X = oc.recon(ves.N, ves.x0, ves.y0, ves.L, ves.theta);
if options.confined
  oc = curve(params.Nbd);
  [~,Xwalls] = oc.initConfig(params.Nbd,false,...
              'scale', 3, ...
              'center', geomCenter, 'geometry', wallGeometry);

  % Build object for the vesicle but without a band-limited opening angle
  walls = capsules(Xwalls,[],params);
else
  walls = [];
end
tt = tstep(params,options,ves,walls); % Shorthand for tstep class
om = monitor(ves.X,walls,params,options); % Shorthand for monitor class


% Take first step with first-order Euler to update for dX and dtheta
% [ves,ux_old,uy_old,L,Ln,dcur0,fntheta,N2Hat,cx0,cy0] = ...
%       tt.FirstSteps(ves,params,options,om);
[ves,ux_old,uy_old,L,Ln,dcur0,fntheta,N2Hat,eta_old] = ...
      tt.FirstSteps(ves,walls,params,options,om);
% Begin time step routine using multistep to update dX and dtheta 
% ves = tt.TimeStepLoop(ves,params,om,ux_old,uy_old,L,Ln,dcur0,...
%                       fntheta,N2Hat,cx0,cy0);
ves = tt.TimeStepLoop(ves,walls,params,om,ux_old,uy_old,L,Ln,dcur0,...
                      fntheta,N2Hat,eta_old);


end % Ves2D
