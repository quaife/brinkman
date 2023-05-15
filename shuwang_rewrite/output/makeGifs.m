clc;clear;close all;
addpath ..

name =  'Figure 6a multicomponent';
name1 = 'Chi2p85_shax3p085_scL0p5393_Conc0p5_Beta5en4_n1024_nbd1024_dt5en6_bmax1_bmin0p1_eps0p04_a100_contracting_left';
%name2 = 'Chip25_shax5p7_scL0p311_Conc0p5_Beta0_n1024_nbd1024_dt1en4_bmax1_bmin0p1_eps0p04_a100_longchoke_rand';


ptitle ={'RA = 0.65, \chi = 5.7\mu/s, u = 0.3, \beta = 1e-4','Confined'};
% ptitle2 ={'\chi = 400, u = 0.5','Confined'};
file1 = [name1 '.bin'];
% file2 = [name2 '.bin'];
ax = [-25 25 -5 5];

[posx,posy,conc,ea,el,time,xvel1,yvel1,ten] = loadFile(file1);
%[posx2,posy2,conc2,ea2,el2,time2,xvel12,yvel12,ten2] = loadFile(file2);

% posx = posx(:,:,1:2:end);
% posy = posy(:,:,1:2:end);
% conc = conc(:,:,1:2:end);
% 
% posx2 = posx2(:,:,1:5:end);
% posy2 = posy2(:,:,1:5:end);
% conc2 = conc2(:,:,1:5:end);
% ea2 = ea2(1:10:end);
% el2 = el2(1:10:end);

bmin1 = 0.1;
bmax1 = 1;
% bmin2 = 0.1;

if 1
    Nbd = 1024;
    geomCenter = [0;0];
    wallGeometry = 'contracting';
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
 
 irate = 50; 
 istart = 1;
 iend = numel(time);
 %iend = find(time == 4);

 N = length(posx(:,:,1));
 %time = time2(istart:iend);
 oc = curve(N);
 h = figure;
 set(gca,'LineWidth',3)
 set(0,'DefaultAxesFontSize',22)
 %set(gcf,'position',[500,0,900,1000])
 set(gcf,'Position',get(0,'Screensize'))
% ax = [-12 12 -4 4];
%  ax = [-4 4 -4 4]; 

  for k = istart:irate:iend
       clf; 
       xx1 = posx(:,:,k);
       yy1 = posy(:,:,k);
%        if k >= numel(time2)
%          xx2 = posx2(:,:,end);
%          yy2 = posy2(:,:,end);
%          rbn2 = (bmin2 - 1)/2*tanh(3*(conc2(:,:,end) - 0.5)) + (1 + bmin2)/2;
%        else
%          xx2 = posx2(:,:,k);
%          yy2 = posy2(:,:,k);
%          rbn2 = (bmin2 - 1)/2*tanh(3*(conc2(:,:,k) - 0.5)) + (1 + bmin2)/2;
%        end
%        
       xxx1 = [xx1(:,:);xx1(1,:)];
       yyy1 = [yy1(:,:);yy1(1,:)];
%        
%        xxx2 = [xx2(:,:);xx2(1,:)];
%        yyy2 = [yy2(:,:);yy2(1,:)];
       
       N = numel(conc(:,:,k));
       rbn = (bmin1 - bmax1)/2*tanh(3*(conc(:,:,k) - 0.5)) + (bmax1 + bmin1)/2;
       ttt1 = [rbn(:,:);rbn(1,:)];
       
       
%        ttt2 = [rbn2(:,:);rbn2(1,:)];
%        
%        subplot(2,1,1)
       plot(xxx1(1,:),yyy1(1,:),'k.','markersize',30);
       hold on
       cline(xxx1,yyy1,ttt1);
       colorbar
       plot(xwalls,ywalls,'r','linewidth',1)
       axis equal
       axis([min(xxx1)-1, max(xxx1)+1,min(yyy1)-1, max(yyy1)+1 ])
       %axis([14.75, 15,-0.5, -0.3])
       set(gca,'xtick',[])
       set(gca,'ytick',[])
       set(gca,'xcolor','white')
       set(gca,'ycolor','white')
       titleStr = ['t = ' num2str(time(k),'%4.2e') ...
         ' eA = ' num2str(ea(k),'%4.2e') ...
         ' eL = ' num2str(el(k),'%4.2e')];
       title([ptitle, titleStr]);
%        hold off
       
%        subplot(2,1,2)
%        plot(xxx2(1,:),yyy2(1,:),'k.','markersize',30);
%        hold on 
%        cline(xxx2,yyy2,ttt2);
%        colorbar
%        set(gca,'linewidth',20)
%        plot(xwalls,ywalls,'r','linewidth',2)
%        axis equal
%        axis([min(xxx2)-1, max(xxx2)+1,min(yyy2)-1, max(yyy2)+1 ])
%        set(gca,'xtick',[])
%        set(gca,'ytick',[])
%        set(gca,'xcolor','white')
%        set(gca,'ycolor','white')
%        if k >= numel(time2)
%          titleStr = ['t = ' num2str(time(end),'%4.2e') ...
%            ' eA = ' num2str(ea2(end),'%4.2e') ...
%            ' eL = ' num2str(el2(end),'%4.2e')];
%          title([ptitle2, titleStr]);
%        else
%          titleStr = ['t = ' num2str(time(k),'%4.2e') ...
%            ' eA = ' num2str(ea2(k),'%4.2e') ...
%            ' eL = ' num2str(el2(k),'%4.2e')];
%          title([ptitle2, titleStr]);
%        end
%        hold off
%        
%        % Capture the plot as an image 
       frame = getframe(h); 
       im = frame2im(frame); 
       [imind,cm] = rgb2ind(im,256); 
       % Write to the GIF File 
       if k == istart 
           imwrite(imind,cm,[name '.gif'],'gif', 'Loopcount',inf); 
       else 
           imwrite(imind,cm,[name '.gif'],'gif','WriteMode','append'); 
       end 

  end
 
