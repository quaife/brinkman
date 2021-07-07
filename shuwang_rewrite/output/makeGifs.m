clc;clear;close all;
addpath ..

name = 'Confined_Chi400_Scale0p35_shortax2p7_conc0p3';
%name2 = 'Confined_Chi400_Scale0p35_shortax2p7_conc0p3';

ptitle ={'\chi = 400, u = 0.3','Confined'};
file1 = [name '.bin'];
%file2 = [name2 '.bin'];
ax = [-10 10 -5 5];

[posx,posy,conc,ea,el,time,xvel1,yvel1,ten] = loadFile(file1);
%[posx2,posy2,conc2,ea2,el2,time2,xvel12,yvel12,ten2] = loadFile(file2);

if 1
    Nbd = 4*128;
    geomCenter = [0;0];
    wallGeometry = 'choke';
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
 
 irate = 10; 
 istart = 2;
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
       %xx2 = posx2(:,:,k);
       %yy2 = posy2(:,:,k);
       
       vec1 = [xx1(:,:);xx1(1,:)];
       vec2 = [yy1(:,:);yy1(1,:)];
       
       %vec3 = [xx1(:,:);xx1(1,:)];
       %vec4 = [yy1(:,:);yy1(1,:)];
       
       N = numel(conc(:,:,k));
       rbn = 1 * (ones(N,1) - conc(:,:,k)) + 0.1*conc(:,:,k);
       vec5 = [rbn(:,:);rbn(1,:)];
       
%        N = numel(conc2(:,:,k));
%        rbn2 = 1 * (ones(N,1) - conc2(:,:,k)) + 0.1*conc2(:,:,k);
%        vec6 = [rbn2(:,:);rbn2(1,:)];
%        
       plot(vec1, vec2, 'r')
       hold on
      % plot(vec3, vec4, 'b--')
       plot(vec1(1,:),vec2(1,:),'k.','markersize',30);
       %plot(vec3(1,:),vec4(1,:),'k.','markersize',30);
       cline(vec1,vec2,vec5);
       set(gca,'linewidth',20)
       %colorbar
       plot(xwalls,ywalls,'r','linewidth',2)
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
       title([ptitle, titleStr]);
       %title(['confined examples ', titleStr]);
       %legend('unconfined','confined')
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
 
