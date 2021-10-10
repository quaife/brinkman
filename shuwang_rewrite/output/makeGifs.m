clc;clear;close all;
addpath ..

name = 'NewBendingMod_constrained_1en6';
name2 = 'Chi400_shax3p45_scL0p49_Conc0p3_Beta0_n1024_nbd1024_dt1en6_bmax1_bmin0p001_eps0p04_longchoke';

ptitle ={'\chi = 400, u = 0.3','Confined, New bending'};
ptitle2 ={'\chi = 400, u = 0.3','Confined, Old bending'};
file1 = [name '.bin'];
file2 = [name2 '.bin'];
ax = [-10 10 -5 5];

[posx,posy,conc,ea,el,time,xvel1,yvel1,ten] = loadFile(file1);
[posx2,posy2,conc2,ea2,el2,time2,xvel12,yvel12,ten2] = loadFile(file2);

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

% plot(mean(squeeze(posy)))
% pause

%plot(time,squeeze(mean(posy)))
 
 irate = 100; 
 istart = 1;
 iend = numel(time2);
 %iend = find(time == 4);
 ntime = iend;
 ntime = numel(time2);
 N = length(posx(:,:,1));
 time = time2(istart:iend);
 oc = curve(N);
 h = figure;
 set(gca,'LineWidth',3)
 set(0,'DefaultAxesFontSize',22)
 %set(gcf,'position',[500,0,900,1000])
 set(gcf,'Position',get(0,'Screensize'))
% ax = [-12 12 -4 4];
 ax = [-4 4 -4 4]; 

  for k = istart:irate:iend
       clf; 
       xx1 = posx(:,:,k);
       yy1 = posy(:,:,k);
       xx2 = posx2(:,:,k);
       yy2 = posy2(:,:,k);
       
       xxx1 = [xx1(:,:);xx1(1,:)];
       yyy1 = [yy1(:,:);yy1(1,:)];
       
       xxx2 = [xx2(:,:);xx2(1,:)];
       yyy2 = [yy2(:,:);yy2(1,:)];
       
       N = numel(conc(:,:,k));
       rbn = 1 * (ones(N,1) - conc(:,:,k)) + 0.1*conc(:,:,k);
       ttt1 = [rbn(:,:);rbn(1,:)];
       
       N = numel(conc2(:,:,k));
       rbn2 = 1 * (ones(N,1) - conc2(:,:,k)) + 0.1*conc2(:,:,k);
       ttt2 = [rbn2(:,:);rbn2(1,:)];
       
       subplot(2,1,1)
       plot(xxx1(1,:),yyy1(1,:),'k.','markersize',30);
       hold on
       cline(xxx1,yyy1,ttt1);
       set(gca,'linewidth',20)
       plot(xwalls,ywalls,'r','linewidth',2)
       axis equal
       axis([min(xxx1)-1, max(xxx1)+1,min(yyy1)-1, max(yyy1)+1 ])
       set(gca,'xtick',[])
       set(gca,'ytick',[])
       set(gca,'xcolor','white')
       set(gca,'ycolor','white')
       titleStr = ['t = ' num2str(time(k),'%4.2e') ...
         ' eA = ' num2str(ea(k),'%4.2e') ...
         ' eL = ' num2str(el(k),'%4.2e')];
       title([ptitle, titleStr]);
       hold off
       subplot(2,1,2)
       plot(xxx2(1,:),yyy2(1,:),'k.','markersize',30);
       hold on 
       cline(xxx2,yyy2,ttt2);
       set(gca,'linewidth',20)
       plot(xwalls,ywalls,'r','linewidth',2)
       axis equal
       axis([min(xxx2)-1, max(xxx2)+1,min(yyy2)-1, max(yyy2)+1 ])
       set(gca,'xtick',[])
       set(gca,'ytick',[])
       set(gca,'xcolor','white')
       set(gca,'ycolor','white')
       titleStr = ['t = ' num2str(time(k),'%4.2e') ...
         ' eA = ' num2str(ea2(k),'%4.2e') ...
         ' eL = ' num2str(el2(k),'%4.2e')];
       title([ptitle2, titleStr]);
       hold off
       
       % Capture the plot as an image 
       frame = getframe(h); 
       im = frame2im(frame); 
       [imind,cm] = rgb2ind(im,256); 
       % Write to the GIF File 
       if k == 1 
           imwrite(imind,cm,[name '.gif'],'gif', 'Loopcount',inf); 
       else 
           imwrite(imind,cm,[name '.gif'],'gif','WriteMode','append'); 
       end 

  end
 
