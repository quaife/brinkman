addpath ../../../../../shuwang_rewrite

% clc;clear;clf;
% addpath ..
% addpath D:\reserach\Runs for paper\stenosis
% file1 = "4a_newBendingModel.bin";
% [posx1,posy1,conc1,ea1,el1,time1,xvel1,yvel1,ten1] = loadFile(file1);
% file2 = "4a_newBendingModel_cont.bin";
% [posx2,posy2,conc2,ea2,el2,time2,xvel2,yvel2,ten2] = loadFile(file2);
% 
% posx1 = posx1(:,:,1:10:end);
% posy1 = posy1(:,:,1:10:end);
% conc1 = conc1(:,:,1:10:end);
% ten1 = ten1(:,:,1:10:end);
% time1 = time1(1:10:end);
% ea1 = ea1(1:10:end);
% el1 = el1(1:10:end);

% posx2 = posx2(:,:,1:10:end);
% posy2 = posy2(:,:,1:10:end);
% conc2 = conc2(:,:,1:10:end);
% ten2 = ten2(:,:,1:10:end);
% time2 = time2(1:10:end);
% ea2 = ea2(1:10:end);
% el2 = el2(1:10:end);
% 
%Runs that restarted bring back one
% posx2 = posx2(:,:,2:end); 
% posy2 = posy2(:,:,2:end);
% conc2 = conc2(:,:,2:end);
% ea2 = ea2(2:end);
% el2 = el2(2:end);
% time2 = time2(2:end);
% xvel2 = xvel2(:,:,2:end);
% yvel2 = yvel2(:,:,2:end);
% ten2 = ten2(:,:,2:end);
% 
% [sx11,sx12,sx13] = size(posx1);
% [sx21,sx22,sx23] = size(posx2);
% posx(1:sx11, sx12, 1:sx13) = posx1;
% posx(1:sx11, sx12, sx13+1:sx23+sx13) = posx2;
% 
% [sy11,sy12,sy13] = size(posy1);
% [sy21,sy22,sy23] = size(posy2);
% posy(1:sy11, sy12, 1:sy13) = posy1;
% posy(1:sy11, sy12, sy13+1:sy23+sy13) = posy2;
% 
% [sc11,sc12,sc13] = size(conc1);
% [sc21,sc22,sc23] = size(conc2);
% conc(1:sc11, sc12, 1:sc13) = conc1;
% conc(1:sc11, sc12, sc13+1:sc23+sc13) = conc2;
% 
%  time2 = time2 + time1(end);
%  time = [time1;time2];

% posx1 = posx; posy1 = posy; conc1=conc; time1 = time;

% pause
bmin = .1;
bmax = 1;
a = 100;

points = [6,12.9,13.5,14.9271,15.6729,15.871,16.5];
%points = [0.05,0.15,0.2,0.28,0.31,0.34,0.4,0.5];

a = find(mean(squeeze(posx1))>=points(1));
b = find(mean(squeeze(posx1))>=points(2));
c = find(mean(squeeze(posx1))>=points(3));
% d = find(mean(squeeze(posx1))>=points(4));
% e = find(mean(squeeze(posx1))>=points(5));
% f = find(mean(squeeze(posx1))>=points(6));
% g = find(mean(squeeze(posx1))>=points(7));
% a = find(time1>=points(1));
% b = find(time1>=points(2));
% c = find(time1>=points(3));
% d = find(time1>=points(4));
% e = find(time1>=points(5));
% f = find(time1>=points(6));
% f = find(time1>=0.34);
% g = find(time1>=0.4);
% h = find(time1>=0.5);
k = [a(1),b(1),c(1),c(1),c(1),c(1)];%,d(1),e(1),f(1),g(1)];
%k = [a(1),b(1),c(1),d(1),e(1),f(1),g(1),h(1)];

Nbd = 2048;
geomCenter = [0;0];
wallGeometry = 'longchoke';
oc = curve(Nbd);
[~,Xwalls] = oc.initConfig(Nbd,false,...
             'scale', 1, ...
             'center', geomCenter, 'geometry', wallGeometry);
[xwalls,ywalls] = oc.getXY(Xwalls);
xwalls = [xwalls(:);xwalls(1)];
ywalls = [ywalls(:);ywalls(1)];

N = length(posx1(:,:,1));
oc = curve(N);
[~, ~, Len] = oc.geomProp([posx1(:,:,1);posy1(:,:,1)]);

eps = 0.04;
figure(4); clf
plot(xwalls,ywalls, 'k', 'linewidth', 3)
axis equal
hold on
squishx = mean(squeeze(posx1));

colormap(turbo);

for i = 1:length(k)
    xx1 = posx1(:,:,k(i));
    yy1 = posy1(:,:,k(i));
    tt1 = conc1(:,:,k(i));
    xxx1 = [xx1(:,:);xx1(1,:)];
    yyy1 = [yy1(:,:);yy1(1,:)];
    ttt1 = [tt1(:,:);tt1(1,:)];
    
    N = numel(conc1(:,:,k(i)));
    rbn = (bmin - bmax)/2*tanh(3*(tt1 - 0.5)) + (bmax + bmin)/2;
    rbn1 = [rbn(:,:);rbn(1,:)];
    rbn
    pause

    s = subplot(1,length(k),i);
    h = cline(xxx1,yyy1,rbn1);
    %plot(xxx1,yyy1,'color', [0.4796 0.0158 0.0106],'linewidth',3) %at top
    %plot(xxx1,yyy1,'color', [0.6332 0.9919 0.2394],'linewidth',3) %at middle
    hold on
    plot(xwalls,ywalls,'k','linewidth',2)
    axis equal
    xlim([mean(xxx1)-2.5, mean(xxx1)+2.5])
    ylim([-2.5 2.5])
%     xlim([-0.5, 0.5])
%     ylim([-0.4,0.4])
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    set(gca,'xcolor','white')
    set(gca,'ycolor','white')
    set(gca,'box','off')
    set(gca, 'visible', 'off')
   % set(findall(gca, 'type', 'text'), 'visible', 'on')
   % h3 = title(points(i));
   % set(h3,'FontSize',20)
   % hold on
%     plot3(posx1(1,:,k(i)),posy1(1,:,k(i)),100,'k.','markersize',20)
%     view(2)
    hsize = s.Position;
end
% h2 = colorbar;
% set(h2, 'limits', ([0 1]))
% set(h2, 'FontSize',16)
% set(s,'Position',hsize);
% pause
f = figure(4);
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'-dpdf','SPMCp5right.pdf','-r0')
