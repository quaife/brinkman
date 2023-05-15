addpath ..
% addpath 'D:\reserach\Runsforpaper\stenosis'

set(gcf,'Position',get(0,'Screensize'));
set(0,'DefaultAxesFontSize',22)
% options.savefig = false;
% clf;

irate = 500; % controls the speed of the visualization

bmax = 1;
bmin = 0.1;

%file = '4a_newBendingModel.bin';
% file = 'Chi400_shax3p45_scL0p49_Conc0p5_Beta0_n1024_nbd1024_dt1en6_bmax1_bmin0p001_eps0p04_longchoke_oldBending.bin';
%[posx1,posy1,conc1,ea1,el1,time1,xvel1,yvel1] = loadFile(file);
  
istart = 1;
iend = numel(time1); 
ntime = numel(time1);
count = 1;

if 1
    Nbd = 1024;
    geomCenter = [0;0];
    wallGeometry = 'longchoke';
    oc = curve(Nbd);
    [~,Xwalls] = oc.initConfig(Nbd,false,...
                 'scale', 1, ...
                 'center', geomCenter, 'geometry', wallGeometry);
    [xwalls,ywalls] = oc.getXY(Xwalls);
else
    xwalls = 0;
    ywalls = 0;
end

   %colormap(winter)
   colormap(turbo);
  % ax = [-1 1 -1 1];
for k = istart:irate:iend

    xx1 = posx1(:,:,k);
    yy1 = posy1(:,:,k);
    tt =  conc1(:,:,k);
    xxx1 = [xx1(:,:);xx1(1,:)];
    yyy1 = [yy1(:,:);yy1(1,:)];
    ttt1 = [tt(:,:);tt(1,:)];
    N = length(ttt1);
    rbn = (bmin - bmax)/2*tanh(3*(ttt1 - 0.5)) + (bmax + bmin)/2;
    
    %h = cline(xxx1,yyy1,rbn);
    %plot(xxx1,yyy1,'color', [0.4796 0.0158 0.0106],'linewidth',3) %at top
    plot(xxx1,yyy1,'color', [0.6332 0.9919 0.2394],'linewidth',3) %at middle

    hold on 
    xwalls = [xwalls(:,:);xwalls(1,:)];
    ywalls = [ywalls(:,:);ywalls(1,:)];
    plot(xwalls, ywalls, 'k', 'linewidth', 3)
    plot3(posx1(1,:,k),posy1(1,:,k),100,'k.','markersize',15)
    view(2)
    %title('15%% floppy, $\Chi = 0.5 \mu$ m/s','interpreter','latex')
    axis equal
    xlim([mean(xxx1)-2.5, mean(xxx1)+2.5])
    ylim([-2.5 2.5])
    %axis(ax)
    %xlim([-0.5 0.5])
    %ylim ([-0.5 0.5])
    %axis([min(xx1)-1, max(xx1)+1,min(yy1)-1, max(yy1)+1 ])
    %set(h,'Units','Inches');
    %axis([min(xx1)-1, max(xx1)+1,min(yy1)-1, max(yy1)+1 ])  
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    set(gca,'xcolor','white')
    set(gca,'ycolor','white')
    axis off
    %set(gca,'visible','off')
   

%     q = colorbar;
%     set(q,'Limits',[0,1],'Ticks',[0, 0.2, 0.4, 0.6, 0.8,1])
    
    %pos = get(h,'Position');
    %set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    filename = ['ver3a_' num2str(count,'%04.f') '.pdf'];
  
    print(gcf, filename,'-dpdf','-fillpage')
    pause(0.01);clf;
    count = count + 1;
end
