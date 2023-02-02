function ves = Ves2D_Continuation(params, options, fileName, restart, rwhere)

addpath output
filename = ['output/' fileName '.bin'];
[posx,posy,conc,~,~,~,~,~,~] = loadFile(filename);
if restart
    if rwhere == 1 %start from end
       % [~,closestIndex] = min(abs(max(squeeze(posx))-10));
	posx(:,:,end) = posx(:,:,end - 5000);
	posy(:,:,end) = posy(:,:,end - 5000);
	conc(:,:,end) = conc(:,:,end - 5000);
    end
end
X = [posx(:,:,end);posy(:,:,end)];
rcon = conc(:,:,end);

%Build object for the vesicle but without a band-limited opening angle
ves = capsules(X,rcon,params);

%Build the walls
if options.confined
  oc = curve(params.Nbd);
  [~,Xwalls] = oc.initConfig(params.Nbd,false,...
              'scale', 1, ...
              'center', params.geomCenter, 'geometry', params.wallGeometry);

  % Build object for the vesicle but without a band-limited opening angle
  walls = capsules(Xwalls,[],params);
else
  walls = [];
end

%Build classes
tt = tstep(params,options,ves,walls); % Shorthand for tstep class
om = monitor(ves.X,walls,params,options); % Shorthand for monitor class

% Take the next time step with first-order Euler to update for dX and
%  dtheta
[ves,ux_old,uy_old,L,Ln,dcur0,fntheta,N2Hat,eta_old] = ...
tt.FirstSteps(ves,walls,params,options,om);
% Begin time step routine using multistep to update dX and dtheta 
ves = tt.TimeStepLoop(ves,walls,params,om,ux_old,uy_old,L,Ln,dcur0,...
                      fntheta,N2Hat,eta_old);
end
