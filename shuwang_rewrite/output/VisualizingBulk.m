addpath ..
load('~/projects/brinkman/vesicle_code/docs/confinedMulti/makefigs/AshleyFiles/Stenosis_RAp4_MCp5.mat');

N = length(posx1(:,:,1));
oc = curve(N);
op = poten(N);

Nbd = 1024;
geomCenter = [0;0];
wallGeometry = 'longchoke';
ocw = curve(Nbd);
[~,XWalls] = ocw.initConfig(Nbd,false,'scale', 1, ...
             'center', geomCenter, 'geometry', wallGeometry);

options.confined = true;
params.shearRate = 0.25;
params.farFieldFlow = 'longchoke';
params.bendsti = 1; %maximum bending stiffness
params.bendratio = 0.1; %ratio between max and min bending 
                                %stiffness
params.viscosityInside = 1;
params.viscosityOutside = 1;
params.SPcoeff = 0; %Semi-permeability coefficient
params.gmresTol = 1e-8; %GMRES tolerance
params.gmresMaxIter = 20;
params.saveRate = 1;
walls = capsules(XWalls,[]);
tol = 1e-4;
%tol = 1e-1;

% care about vesicle when it's at the third, fourth, fifth, and sixth of
% these values
meancenters = linspace(-22,22,8);

centers = squeeze(mean(posx1));
[~,k] = min(abs(centers - meancenters(6)));

xx1 = posx1(:,:,k);
yy1 = posy1(:,:,k);
concc1 = conc1(:,:,k);

ves = capsules([xx1;yy1],concc1,params); 
tt = tstep(params,options,ves,walls);
rhsWalls = tt.bgFlow(walls.X, tt.shearRate, tt.farFieldFlow);
eta = tt.etaSolver(rhsWalls);

for i = 1:length(k)
  xx1 = posx1(:,:,k(i));
  yy1 = posy1(:,:,k(i));
  concc1 = conc1(:,:,k(i));
  
  ves = capsules([xx1;yy1],concc1,params); 
  tt = tstep(params,options,ves,walls);

  [xtar,ytar] = meshgrid(min(xx1)-1:.01:max(xx1)+1,-.69:.01:.69);
  xtar = xtar(:); ytar = ytar(:);
  targets.N = numel(xtar);
  targets.X = [xtar;ytar];

  % set large value for relative norm and set counter to 0
  relnorm = 1;
  counter = 0;
  % Picard iterate until convergence in the density function on outer
  % wall is found
  while relnorm > tol
    disp([counter relnorm])
    eta_old = eta;
    [uxvel,uyvel,~,eta,uxtar,uytar,pressure] = tt.usetself(ves,eta,targets);
    relnorm = norm(eta-eta_old)/norm(eta_old);
    counter = counter + 1;
  end
  
  % uxtar = uxtar - mean(squeeze(xvel1(:,:,k)));
  figure(2); clf; hold on
  quiver(xtar,ytar,uxtar-mean(uxvel),uytar)
  %quiver(xtar,ytar,uxtar,uytar)
  plot(walls.X(1:end/2),walls.X(end/2+1:end),'k','linewidth',3)
  plot(ves.X(1:end/2),ves.X(end/2+1:end),'r','linewidth',3)
 % title(k(i))
  axis equal
  %axis([min(xx1)-1, max(xx1)+1,min(yy1)-1, max(yy1)+1 ])
  xlim([mean(xx1)-2, mean(xx1)+2])
  ylim([-.7 .7])
  set(gca,'xtick',[])
  set(gca,'ytick',[])
  set(gca,'xcolor','white')
  set(gca,'ycolor','white')
  title(time1(k(i)))
  %pause
end

ny = 139;
nx = numel(xtar)/ny;
posx = posx1(:,:,k);
posy = posy1(:,:,k);
conc = conc1(:,:,k);
ten = ten1(:,:,k);
xtar = reshape(xtar,ny,nx);
ytar = reshape(ytar,ny,nx);
uxtar = reshape(uxtar,ny,nx);
uytar = reshape(uytar,ny,nx);
wallsx = walls.X(1:end/2);
wallsy = walls.X(end/2+1:end);

filename = 'Stenosis_RAp4_MCp5_pos4.mat';
save(filename,'posx','posy','conc','ten','eta','xtar','ytar','uxtar','uytar',...
      'wallsx','wallsy');

