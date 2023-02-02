addpath ..
%file1 = 'Chi2p5_shax5p7_scL0p311_Conc0p5_Beta0_n1024_nbd1024_dt5en7_bmax1_bmin0p1_eps0p04_a100_contracting_left.bin';
%[posx1,posy1,conc1,ea1,el1,time1,xvel1,yvel1,ten1] = loadFile(file1);
% posx1= posx3;
% posy1 = posy3;
% time1 = time3;
% conc1 = conc3;



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
% a = find(time1>=28);
% 
% b = find(time1>=53);

istart = 35000;
irate = 1;%5000;
iend = 35000;%76018;%numel(time1);
count = 1;


xx1 = posx1(:,:,istart);
yy1 = posy1(:,:,istart);
concc1 = conc1(:,:,istart);
ves = capsules([xx1;yy1],concc1,params); 

tt = tstep(params,options,ves,walls);
rhsWalls = tt.bgFlow(walls.X, tt.shearRate, tt.farFieldFlow);
eta = tt.etaSolver(rhsWalls);

ytar = 0;%linspace(-0.5,0.5,100)';
xtar = 26*ones(size(ytar));
ytar = [ytar;ytar];
xtar = [-xtar;xtar];
targets.N = numel(xtar);
targets.X = [xtar;ytar];

rhsWalls = tt.bgFlow(walls.X, tt.shearRate, tt.farFieldFlow);
etawithout = tt.etaSolver(rhsWalls);
[~,~,~,~,~,~,pressureVoid] = tt.usetself(ves,etawithout,targets);
PressDropNoVes = pressureVoid(end/2+1:end) - pressureVoid(1:end/2);


for k = istart:irate:iend
xx1 = posx1(:,:,k);
yy1 = posy1(:,:,k);
concc1 = conc1(:,:,k);

ves = capsules([xx1;yy1],concc1,params); 
tt = tstep(params,options,ves,walls);

relnorm = 1;
counter = 0;
while relnorm > tol
  eta_old = eta;
  [~,~,~,eta,uxtar(:,k),uytar(:,k),pressure] = tt.usetself(ves,eta,targets);
  relnorm = norm(eta-eta_old)/norm(eta_old);
  counter = counter + 1;
end

PressDropVes(count) = pressure(2) - pressure(1);
pressDrop(count) = mean(PressDropVes-PressDropNoVes);
count = count + 1;

end
hold on
%plot(mean(squeeze(posx1(:,:,istart:irate:iend)))/20,-pressDrop*1e-1*0.000145,'r')
% plot(mean(squeeze(posx1(:,:,istart:irate:iend))),-pressDrop*1e-1,'b','linewidth',3)