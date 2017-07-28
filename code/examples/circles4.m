clear all
addpath ../src

load radii4.dat;
load centers4.dat;
radii = radii4;
centers = centers4;

prams.Nouter = 512;
% number of pois on outer solid wall
prams.Ninner = 64;
% number of points per circle exclusion
prams.nv = numel(radii);
% number of exclusions
prams.gmresTol = 1e-8;
% gmres tolerance
prams.maxIter = min(2*(prams.Nouter + prams.nv*prams.Ninner),5000);
% maximum number of gmres iterations

% Different options
options.bieSolve = true; 
options.computeEuler = true;
options.axis = [-1 1 -0.8 0.8];
options.dataFile = 'output/circles4Data.bin';
options.logFile = 'output/circles4.log';
options.farField = 'pipe';
options.fmm = true;
options.saveData = true;
options.usePlot = true;
options.verbose = true;

oc = curve;
Xouter = oc.initConfig(prams.Nouter,'square', ...
         'outerDimensions',[-3 3 -1 1]);
% outer most boundary

Xinner = oc.initConfig(prams.Ninner,'circles', ...
          'nv',prams.nv, ...
          'center',centers, ...
          'radii',radii);
% inner geometry pore(s)

om = monitor(options,prams);
if options.usePlot
  om.plotGeometry(Xinner,Xouter);
  drawnow()
end
%theta = linspace(0,2*pi,1000);
%hold on;
%plot(0.7*exp(1i*theta),'r--')
%pause
%% plot the figure if desired


if options.bieSolve
  stokesSolver(Xinner,Xouter,options,prams);
end
% solve density function and write to .bin files.  If this
% calculation has already been done, Eulerian grid points can be
% computed without solving for the density function again

if options.computeEuler
%  prams.xmin = min(Xouter(1:end/2)) + 0.1;
%  prams.xmax = max(Xouter(1:end/2)) - 0.1;
%  prams.ymin = min(Xouter(end/2+1:end)) + 0.05;
%  prams.ymax = max(Xouter(end/2+1:end)) - 0.05;
  prams.xmin = -1;
  prams.xmax = 1;
  prams.ymin = -0.9;
  prams.ymax = 0.9;
  prams.nx = 500;
  prams.ny = 500;

  [eX,eY,u,v] = computeEulerVelocity(options,prams,...
      radii,centers);

  clf; hold on
  plot(Xouter(1:end/2),Xouter(end/2+1:end))
  quiver(eX,eY,u,v)
  fill(Xinner(1:end/2,:),Xinner(end/2+1:end,:),'k')
  axis equal
end


