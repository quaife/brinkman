function [] = MakeSnaps(file,sp,cb,ha)
%if ~exist('posx1')
%  load Stenosis_RAp6_MCp5.mat
%end
load(file);

%addpath ../../../../../shuwang_rewrite
addpath ~/projects/brinkman/shuwang_rewrite

bmin = .1;
bmax = 1;
a = 100;

points = linspace(-22,22,8);
%points = [-25,12.9,13.5,14.9271,15.6729,15.871,16.5];
%points = [0.05,0.15,0.2,0.28,0.31,0.34,0.4,0.5];

a = find(mean(squeeze(posx1))>=points(1));
b = find(mean(squeeze(posx1))>=points(2));
c = find(mean(squeeze(posx1))>=points(3));
d = find(mean(squeeze(posx1))>=points(4));
e = find(mean(squeeze(posx1))>=points(5));
f = find(mean(squeeze(posx1))>=points(6));
g = find(mean(squeeze(posx1))>=points(7));
h = find(mean(squeeze(posx1))>=points(8));
k = [a(1),b(1),c(1),c(1),c(1),c(1),d(1),e(1),f(1),g(1),h(1)];

Nbd = 2048;
geomCenter = [0;0];
wallGeometry = 'longchoke';
oc = curve(Nbd);
[~,Xwalls] = oc.initConfig(Nbd,false,...
             'scale', 1, ...
             'center', geomCenter, 'geometry', wallGeometry);
[xwalls,ywalls] = oc.getXY(Xwalls);
xwalls = [xwalls(:);xwalls(1:4)];
ywalls = [ywalls(:);ywalls(1:4)];

N = length(posx1(:,:,1));
oc = curve(N);

%subplot(3,1,sp);
%figure(1); clf
plot(xwalls,ywalls, 'k', 'linewidth', 2)
axis equal
hold on

color1 = [0 0 1];
%color2 = [0 1 1];
color2 = [0.5 1 0];
color3 = [1 0 0];
t = linspace(0,1,256)';
map1 = (1-t)*color1 + t*color2;
map2 = (1-t(2:end))*color2 + t(2:end)*color3;
map = [map1;map2];
colormap(map);

for i = 1:length(k)
  xx1 = posx1(:,:,k(i));
  yy1 = posy1(:,:,k(i));
  tt1 = conc1(:,:,k(i));
  xxx1 = [xx1(:,:);xx1(1:4,:)];
  yyy1 = [yy1(:,:);yy1(1:4,:)];
  ttt1 = [tt1(:,:);tt1(1:4,:)];

  if strcmp(file,'AshleyFiles/Stenosis_RAp6_SC.mat') || ...
     strcmp(file,'AshleyFiles/Stenosis_RAp4_SC.mat')
%    rbn = bmax*ones(N,1);
%    rbn1 = [rbn(:,:);rbn(1:4,:)];
%    h = cline(xxx1,yyy1,rbn1);
%    set(h,'linewidth',3);
    plot(xxx1,yyy1,'color',map(end,:),'linewidth',3);
  elseif strcmp(file,'AshleyFiles/Stenosis_RAp6_SCp55.mat') || ...
         strcmp(file,'AshleyFiles/Stenosis_RAp4_SCp55.mat')
%    rbn = (bmax+bmin)/2*ones(N,1);
%    rbn1 = [rbn(:,:);rbn(1:4,:)];
%    h = cline(xxx1,yyy1,rbn1);
%    set(h,'linewidth',3);
    plot(xxx1,yyy1,'color',map(256,:),'linewidth',3);
  else
    rbn = (bmin - bmax)/2*tanh(3*(tt1 - 0.5)) + (bmax + bmin)/2;
    rbn1 = [rbn(:,:);rbn(1:4,:)];
    h = cline(xxx1,yyy1,rbn1);
    set(h,'linewidth',3);
  end
  plot3(xxx1(1),yyy1(1),100,'k.','markersize',20)
end
xlim([-25 25])
ylim([-2.5 2.5])

set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'xcolor','white')
set(gca,'ycolor','white')
set(gca,'box','off')
set(gca, 'visible', 'off')

h2 = colorbar;
set(h2, 'limits', ([0.1 1]))
%set(h2,'xtick',(0.2:0.2:1));
set(h2,'xtick',[0.1 1]);
set(h2, 'FontSize',16)
if ~cb
  set(h2,'visible','off')
end
%set(s,'Position',hsize);
%f = figure(1);
%set(f,'Units','Inches');
%pos = get(f,'Position');
%set(f,'PaperPositionMode','Auto','PaperUnits','Inches',...
%    'PaperSize',[pos(3), pos(4)])
%print(f,'-dpdf','SPMCp5right.pdf','-r0')
