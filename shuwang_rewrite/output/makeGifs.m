clc;clear;close all;
addpath ..

name = 'Chi200_RA0p85_Conc0p3_Beta0_y0p1_eps0p04_n20';
ptitle ={'\chi = 200, RA_0 = 0.85, u = 0.3,', '\beta = 0, \epsilon = 0.04'};
file = [name '.bin'];
ax = [-2 2 -2 2];

[posx,posy,conc,ea,el,time,xvel1,yvel1,ten] = loadFile(file);
% plot(mean(squeeze(posy)))
% pause

%plot(time,squeeze(mean(posy)))
 
 irate = 100; 
 istart = 1;
 iend = numel(time);
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
 ax = [-2 2 -2 2];
 
  for k = istart:irate:iend
       clf; 
       xx1 = posx(:,:,k);
       yy1 = posy(:,:,k);
       
       vec1 = [xx1(:,:);xx1(1,:)];
       vec2 = [yy1(:,:);yy1(1,:)];
       
       N = numel(conc(:,:,k));
       rbn = 1 * (ones(N,1) - conc(:,:,k)) + 0.1*conc(:,:,k);
       vec4 = [rbn(:,:);rbn(1,:)];
       
       
       cline(vec1,vec2,vec4);
       set(gca,'linewidth',3)
       colorbar
       hold on
       plot(vec1(1,:),vec2(1,:),'k.','markersize',30);
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
 
