function Xfinal  = Ves2D_hack(X,Xwalls,prams,options,Xtra,pressTar)
% Ves2D does time stepping on the intial configuration X
% with parameters and options defined in prams and options.
% Also can pass a set of initial tracer locations (Xtra) and 
% target points where one wants the pressure and stress (pressTar)

global matvecs ;  % number of matvecs
global derivs  ;  % number of times we compute differential
                  % operators for preconditioning
global fmms       % number of fmm calls

if nargin == 4
  Xtra = [];       % Xtra is positions of tracers for visualization
  pressTar = [];   % points to compute pressure for postprocessing
elseif nargin == 5
  pressTar = [];   
end

matvecs = 0; 
% counter for the total number of time steps
derivs = 0; 
% counter for the total number of times derivative operators are
% computed
fmms = 0; 
% counter for the total number of FMM calls
              
nreject = 0; % count for the number of rejected time steps
naccept = 0; % count for the number of accepted time steps
countGMRES = 0; % count for the number of total GMRES iterations

N = prams.N; % Number of points per vesicle
nv = prams.nv; % Number of vesicles
Nbd = prams.Nbd; % number of points per solid wall
nvbd = prams.nvbd; % number of solid wall components

om = monitor(X,Xwalls,options,prams);
% Create class for doing input/output

if options.profile
  profile off; profile on -timer real;
%  profile off; profile on
end
% Turn on profiler

tstart = tic;
% Start a timer
c = clock;
% get current time in case tic and toc give funny anwers
message = ['Initial Month is ' num2str(c(2),'%d')];
om.writeMessage(message,'%s\n')
message = ['Initial Day is ' num2str(c(3),'%d')];;
om.writeMessage(message,'%s\n')
message = ['Initial Hour is ' num2str(c(4),'%d')];;
om.writeMessage(message,'%s\n')
message = ['Initial Minute is ' num2str(c(5),'%d')];;
om.writeMessage(message,'%s\n')
message = ['Initial Second is ' num2str(round(c(6)),'%d')];;
om.writeMessage(message,'%s\n')
message = ' ';
om.writeMessage(message,'%s\n')

oc = curve;
tt = tstep(options,prams);
% build an object of class tstep with required options and parameters
Xstore = zeros(2*N,nv,options.order);
sigStore = zeros(N,nv,options.order);
uStore = zeros(2*N,nv,options.order);
etaStore = zeros(2*prams.Nbd,prams.nvbd,options.order);
RSstore = zeros(3,prams.nvbd,options.order);
% need to store options.order previous time steps to do
% higher-order time stepping
% RSstore is the rotlets and stokeslets stored as
% [stokeslet1;stokeslet2;rotlet]
Xstore(:,:,1) = X;
% initial configuration from user

if ~options.confined
  uStore(:,:,1) = tt.farField(X);
  etaStore = [];
  RSstore = [];
  walls = [];
  wallsCoarse = [];
else
  [walls,wallsCoarse] = tt.initialConfined(prams,Xwalls); 
end
% initial velocity on the vesicles is set to the background velocity.
% If flow is unbounded, there is no density function eta.  If bounded,
% compute a structure for the solid walls 
%
[Xstore,sigStore,uStore,etaStore,RSstore,Xtra] = ...
    tt.firstSteps(options,prams,...
    Xstore(:,:,end),sigStore(:,:,end),uStore(:,:,end),...
    walls,wallsCoarse,om,Xtra,pressTar);
% For higher-order methods (only 2nd order for now), need to initially
% take smaller time steps so that we have two initial conditions

time = (options.order-1)*tt.dt;
% initial time.  firstSteps took the first time steps so that there is
% enough data to use the time stepping order that is desired

accept = true;
% Main time stepping loop
while time < prams.T - 1e-10
  if time+tt.dt > prams.T
    tt.dt = prams.T - time;
  end
  % make sure we land exactly on the time horizon
  dt = tt.dt;
  tt.currentTime = time;
  time = time + dt; % current time
  
  tTstep = tic;
  [X,sigma,u,eta,RS,iter,accept,dtScale,res,iflag] = ...
      tt.timeStepGL(Xstore,sigStore,uStore,...
          etaStore,RSstore,prams.kappa,...
          prams.viscCont,walls,wallsCoarse,om,time,accept);
  countGMRES = countGMRES + iter;
  tTstep = toc(tTstep);
    
  if options.profile
    fprintf('Time to correct area and length     %5.1e\n',toc);
  end

  if ~accept
    time = time - dt;
  end
  % go back to old time

  if 1
  for k = 1:prams.nv
    z = X(1:end/2,k) + 1i*X(end/2+1:end,k);
    zh = fftshift(fft(z));
    zh(1) = 0;
    z = ifft(ifftshift(zh));
    X(1:end/2,k) = real(z);
    X(end/2+1:end,k) = imag(z);
  end
  % remove Nyquist Fourier mode

  xmean1 = mean(X(1:end/2,1));
  xmean2 = mean(X(1:end/2,2));
  xmid = xmean1 + xmean2;
  X(1:end/2,1) = X(1:end/2,1) - xmid/2;
  X(1:end/2,2) = X(1:end/2,2) - xmid/2;
  % shift horizontally so that the y axis is centered between the
  % vesicles midpoints

  ymean1 = mean(X(end/2+1:end,1));
  ymean2 = mean(X(end/2+1:end,2));
  ymid = ymean1 + ymean2;
  X(end/2+1:end,1) = X(end/2+1:end,1) - ymid/2;
  X(end/2+1:end,2) = X(end/2+1:end,2) - ymid/2;
  % shift vertically so that the x axis is centered between the
  % vesicles midpoints
  end

  if accept
    vesicle = capsules(X,sigma,u,prams.kappa,prams.viscCont);
    [shearStress,normalStress] = ...
        vesicle.computeShearStress(options,prams);
    % compute the shear and normal stress along the vesicles

    terminate = om.outputInfo(X,sigma,u,eta,RS,...
        Xwalls,Xtra,time,iter,dtScale,res,iflag);
    % check if we have violated the error in area or length also plot
    % and save the current solution to dat file.  Print information to
    % log file and console
    message = ['A time step takes ' num2str(tTstep,'%2.2e') ' seconds'];
    om.writeMessage(message,'%s\n')
    
    if terminate
      Xfinal = X;
      totTime = toc(tstart);
      message = ['Final Month is ' num2str(c(2),'%d')];
      om.writeMessage(message,'%s\n')
      message = ['Final Day is ' num2str(c(3),'%d')];
      om.writeMessage(message,'%s\n')
      message = ['Final Hour is ' num2str(c(4),'%d')];
      om.writeMessage(message,'%s\n')
      message = ['Final Minute is ' num2str(c(5),'%d')];
      om.writeMessage(message,'%s\n')
      message = ['Final Second is ' num2str(round(c(6)),'%d')];
      om.writeMessage(message,'%s\n')
      message = ' ';
      om.writeMessage(message,'%s\n')
      om.summary(matvecs,derivs,fmms,...
          naccept,nreject,countGMRES,totTime);
      % write a final summary with number of matvecs, accepts/rejects,
      % and the CPU time
      return;
    end
    % if error in area or length is too large, stop simulation
    naccept = naccept + 1;
  else
    nreject = nreject + 1;
  end % if accept
  % save data if solution was accepted, compute pressure and stress

  for k = 1:options.order-1
    Xstore(:,:,k) = Xstore(:,:,k+1);
    sigStore(:,:,k) = sigStore(:,:,k+1);
    uStore(:,:,k) = uStore(:,:,k+1);
    etaStore(:,:,k) = etaStore(:,:,k+1);
    RSstore(:,:,k) = RSstore(:,:,k+1);
  end
  % Save the time steps still required if using second
  % order time stepping
  Xstore(:,:,options.order) = X;
  sigStore(:,:,options.order) = sigma;
  uStore(:,:,options.order) = u;
  etaStore(:,:,options.order) = eta;
  RSstore(:,:,options.order) = RS;
  % update the positions, tension, and velocity field of
  % the vesicles, and the density function and rotlet and
  % stokeslets
end
% end of main 


totTime = toc(tstart);
% Total time doing time stepping
c = clock;
% get current time in case tic and toc give funny anwers
message = ['Final Month is ' num2str(c(2),'%d')];
om.writeMessage(message,'%s\n')
message = ['Final Day is ' num2str(c(3),'%d')];
om.writeMessage(message,'%s\n')
message = ['Final Hour is ' num2str(c(4),'%d')];
om.writeMessage(message,'%s\n')
message = ['Final Minute is ' num2str(c(5),'%d')];
om.writeMessage(message,'%s\n')
message = ['Final Second is ' num2str(round(c(6)),'%d')];
om.writeMessage(message,'%s\n')
message = ' ';
om.writeMessage(message,'%s\n')

if options.profile
  profile off;
  p = profile('info');
  filename = [options.logFile(1:end-4) 'Profile'];
  save([filename '.mat'],'p');
%  profview
  profsave(profile('info'),filename);
end
% Save profile

om.summary(matvecs,derivs,fmms,...
    naccept,nreject,countGMRES,totTime);
% write a final summary with number of matvecs, accepts/rejects,
% and the CPU time

Xfinal = X;
% final configuation

message = 'SIMULATION SUCCESSFULLY COMPLETED';
om.writeMessage(message,'%s\n')
om.writeStars

