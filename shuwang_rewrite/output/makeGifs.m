clc;clear;close all;
addpath ..

name = 'confined_examples';
name1 = 'confined_larger_800_conc0_T10';
name2 = 'unconfined_larger_800_conc0_T10';
ptitle ={'\chi = 800, u = 0.3','confined'};
file1 = [name1 '.bin'];
file2 = [name2 '.bin'];
ax = [-2 2 -2 2];

[posx,posy,conc,ea,el,time,xvel1,yvel1,ten] = loadFile(file1);
[posx2,posy2,conc2,ea2,el2,time2,xvel12,yvel12,ten2] = loadFile(file2);

if 1
    Nbd = 512;
    geomCenter = [0;0];
    wallGeometry = 'tube';
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
 
 irate = 500; 
 istart = 1;
 iend = numel(time);
 %iend = find(time == 4);
 ntime = iend;
 ntime = numel(time);
 N = length(posx(:,:,1));
 time = time(istart:iend);
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
       
       vec1 = [xx1(:,:);xx1(1,:)];
       vec2 = [yy1(:,:);yy1(1,:)];
       
       N = numel(conc(:,:,k));
       rbn = 1 * (ones(N,1) - conc(:,:,k)) + 0.1*conc(:,:,k);
       vec4 = [rbn(:,:);rbn(1,:)];
       
       xx2 = posx2(:,:,k);
       yy2 = posy2(:,:,k);
       
       vec5 = [xx2(:,:);xx2(1,:)];
       vec6 = [yy2(:,:);yy2(1,:)];
       
       N = numel(conc2(:,:,k));
       rbn2 = 1 * (ones(N,1) - conc2(:,:,k)) + 0.1*conc2(:,:,k);
       vec7 = [rbn2(:,:);rbn2(1,:)];
       
       
       
       %cline(vec1,vec2,vec4);
       plot(vec1,vec2, 'r')
       set(gca,'linewidth',20)
       %colorbar
       hold on
       plot(vec5,vec6, 'r--')
       %cline(vec5,vec6,vec7);
       plot(vec1(1,:),vec2(1,:),'k.','markersize',30);
       plot(vec5(1,:),vec6(1,:),'k.','markersize',30);
      % plot(xwalls,ywalls,'r','linewidth',2)
       plot([ax(1) ax(2)],[0 0],'k--')
       axis equal
       axis(ax)
      
       set(gca,'xtick',[])
       set(gca,'ytick',[])
       set(gca,'xcolor','white')
       set(gca,'ycolor','white')
       titleStr = ['t = ' num2str(time(k),'%4.2e') ...
         ' eA = ' num2str(ea(k),'%4.2e') ...
         ' eL = ' num2str(el(k),'%4.2e')];
       %title([ptitle, titleStr]);
       title(['confined examples ', titleStr]);
       legend('confined','unconfined')
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
 
