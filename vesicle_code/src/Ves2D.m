function Xfinal = Ves2D(X,Xwalls,prams,options,pressTar)
% Ves2D does time stepping on the intial configuration X
% with parameters and options defined in prams and options.
% Also can pass a set of initial target points where one wants 
% the pressure and stress (pressTar)

global matvecs ;  % number of matvecs
global derivs  ;  % number of times we compute differential
                  % operators for preconditioning
global fmms       % number of fmm calls


% TODO: THIS SHOULDN'T GO HERE
pressDrop = 200;
%pressDrop = 5;


om = monitor(X,Xwalls,options,prams);
% Create class for doing input/output

if nargin == 4
  pressTar = [];   % points to compute pressure for postprocessing
else
  om.initializePressure(pressTar); 
end
Ntar = 21;
ytar = linspace(-3,3,Ntar)';
dy = ytar(2) - ytar(1);
pressTar = [[0.5*ones(Ntar,1);ytar] [17.5*ones(Ntar,1);ytar]]; 

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

tt = tstep(options,prams);
% build an object of class tstep with required options and parameters

Xstore = zeros(2*N,nv);
sigStore = zeros(N,nv);
uStore = zeros(2*N,nv);
etaStore = zeros(2*prams.Nbd,prams.nvbd);
RSstore = zeros(3,prams.nvbd);
% palce to store previous time step to do time stepping
% RSstore is the rotlets and stokeslets stored as
% [stokeslet1;stokeslet2;rotlet]
Xstore = X;
% initial configuration from user

if ~options.confined
  uStore = tt.farField(X);
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
[Xstore,sigStore,uStore,etaStore,RSstore] = ...
    tt.firstSteps(options,prams,Xstore,sigStore,uStore,...
    walls,wallsCoarse,om);
% Was for higher-order multistep methods which have been fazed out. Now
% it just computes the initial density function and tension that are
% needed for SDC


% compute pressure and stress
if nargin == 5
  op = poten(walls.N);
  [~,pressDLPtar] = op.exactPressDL(walls,etaStore,[],pressTar,1);
  pressDLPtar = pressDLPtar(1:end/2);

  vesicle = capsules(Xstore,sigStore,uStore,...
    prams.kappa,prams.viscCont,options.semipermeable,...
    prams.fluxCoeff,options.fluxShape);
  tracJump = vesicle.tracJump(Xstore,sigStore);
  % compute traction

  op = poten(vesicle.N);
  [~,pressSLPtar] = op.exactPressSL(vesicle,tracJump,[],pressTar,1);
  pressSLPtar = pressSLPtar(1:end/2);

  om.writePressure(pressDLPtar + pressSLPtar);
  % write the pressure contributions to file
end

time = 0;
% initial time

accept = true;
% Main time stepping loop
nstep = 0;
while time < prams.T - 1e-10
%%Hacking for time-varying periodic flow 02/21/2020
%  tt.farField = @(X) tt.bgFlow(X,options.farField,...
%	options.farFieldSpeed*(1+sin(time))/2);
%%Hacking for time-varying periodic flow 02/21/2020
  if time+tt.dt > prams.T
    tt.dt = prams.T - time;
  end
  % make sure we land exactly on the time horizon
  dt = tt.dt;
  tt.currentTime = time;
  time = time + dt; % current time
  
  tTstep = tic;

% find the protein locations in parameter space with respect to the
% memebranes tracker points

  vesicle = capsules(Xstore,sigStore,uStore,...
      prams.kappa,prams.viscCont,options.semipermeable,...
      prams.fluxCoeff,options.fluxShape);

  % TODO: THIS SHOULD BE A FUNCTION IN TSTEP
  if 1
    [ssig,eeta] = vesicle.computeSigAndEta(tt,walls);

    [~,pressDLPtar] = tt.opWall.exactPressDL(walls,eeta,[],pressTar,1);

    tracJump = vesicle.tracJump(Xstore,ssig);
    % compute traction
    figure(2); clf;
    plot(ssig)
    figure(1)

    [~,pressSLPtar] = tt.op.exactPressSL(vesicle,tracJump,[],pressTar,1);

    pressL = pressSLPtar(1:end/2,1) + pressDLPtar(1:end/2,1);
    pressR = pressSLPtar(1:end/2,2) + pressDLPtar(1:end/2,2);
    avePressL = dy/(ytar(end) - ytar(1))* ...
        (0.5*pressL(1) + sum(pressL(2:end-1)) + 0.5*pressL(end));
    avePressR = dy/(ytar(end) - ytar(1))* ...
        (0.5*pressR(1) + sum(pressR(2:end-1)) + 0.5*pressR(end));

    dpress = avePressR - avePressL;
    
    options.farFieldSpeed = -pressDrop/dpress;
    tt.farField = @(X) tt.bgFlow(X,...
        options.farField,'k',options.farFieldSpeed);
    [walls,wallsCoarse] = tt.initialConfined(prams,Xwalls); 
%    dpress
%    [dpress pressDrop options.farFieldSpeed max(walls.u)]
%    clf;
%    plot(walls.u)
%    pause
    % change velocity field speed so that it maintains a constant
    % pressure drop
  end

    
  [X,sigma,u,eta,RS,iter,accept,dtScale,res,iflag] = ...
      tt.timeStepGL(vesicle,etaStore,RSstore,...
          walls,wallsCoarse,om,time,accept);
  countGMRES = countGMRES + iter;
  tTstep = toc(tTstep);

  if options.profile
    fprintf('Time to correct area and length     %5.1e\n',toc);
  end

  if ~accept
    time = time - dt;
  end
  % go back to old time

  if 0
  xmid = mean(X(1:end/2,1));
  ymid = mean(X(end/2+1:end,1));
  X(1:end/2) = X(1:end/2) - xmid;
  X(end/2+1:end) = X(end/2+1:end) - ymid;
  end
  % Shift single vesicle as in a fluid trap

  if 0
  for k = 1:prams.nv
    z = X(1:end/2,k) + 1i*X(end/2+1:end,k);
    zh = fftshift(fft(z));
    zh(1) = 0;
    z = ifft(ifftshift(zh));
    X(1:end/2,k) = real(z);
    X(end/2+1:end,k) = imag(z);
  end
  end
  % remove Nyquist Fourier mode

  if 0
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
  % shift vertically so that the x axis is centered between the vesicles
  end
  % Shift vesicle doublet as in a fluid trap

  if accept
    nstep = nstep + 1;
    if (mod(nstep,prams.saveRate) == 1 || prams.saveRate == 1)
      om.saveData = true;
    else
      om.saveData = false;
    end
    terminate = om.outputInfo(X,sigma,u,eta,RS,...
        Xwalls,time,iter,dtScale,res,iflag);
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
  % save data if solution was accepted

  if options.xshift && max(X(1:end/2)) > options.xshiftLoc
    X(1:end/2) = X(1:end/2) - options.xshiftVec;
  end
  % Shift vesicle back xshiftVec units

  Xstore = X;
  sigStore = sigma;
  uStore = u;
  etaStore = eta;
  RSstore = RS;
  % update the positions, tension, and velocity field of the vesicles,
  % and the density function and rotlet and stokeslets

  % compute pressure and stress
  if (accept && nargin == 5)
    op = poten(walls.N);
    [~,pressDLPtar] = op.exactPressDL(walls,etaStore,[],pressTar,1);
    pressDLPtar = pressDLPtar(1:end/2);

    vesicle = capsules(Xstore,sigStore,uStore,...
      prams.kappa,prams.viscCont,options.semipermeable,...
      prams.fluxCoeff,options.fluxShape);
    tracJump = vesicle.tracJump(Xstore,sigStore);
    % compute traction

    op = poten(vesicle.N);
    [~,pressSLPtar] = op.exactPressSL(vesicle,tracJump,[],pressTar,1);
    pressSLPtar = pressSLPtar(1:end/2);

    om.writePressure(pressDLPtar + pressSLPtar);
    % write the pressure contributions to file
  end

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


