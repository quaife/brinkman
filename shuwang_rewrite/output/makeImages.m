%clc;clear;close all
clf
addpath ..
%set(gcf,'Position',get(0,'Screensize'));
set(0,'DefaultAxesFontSize',22)
options.savefig = false;

%file = 'extensional_RAp8_Conc0p3_Chi10_beta0_n500_dt0en4.bin';
%file = 'parabolic_Conc0p3_Chi600_beta0_n500_dt0en4_h2.bin';

%file = 'Parabolic_RA0p95_Conc0_Chi600_beta0.bin';
%file = 'Parabolic_RA0p85_Conc0_Chi600_beta0.bin';
%file = 'Parabolic_RA0p75_Conc0_Chi600_beta0.bin';
%file = 'Parabolic_RA0p6_Conc0_Chi600_beta0.bin';
%file = 'Parabolic_RA0p95_Conc0p3_Chi600_beta0.bin';
file = 'test2.bin';
ax = [-2 2 -2 2];

[posx,posy,conc,ea,el,time,xvel1,yvel1,ten] = loadFile(file);

%plot(time,squeeze(mean(posy)))
 
 irate = 1000; 
 istart = 1;
 iend = numel(time);
 ntime = iend;
 ntime = numel(time);
 
 time = time(istart:iend);
 oc = curve;
 
  for k = istart:irate:iend
        clf; 
        xx1 = posx(:,:,k);
        yy1 = posy(:,:,k);
        tt = conc(:,:,k);
 %       [~,area(k),~] = oc.geomProp([xx1;yy1]); 
 % 
        vec1 = [xx1(:,:);xx1(1,:)];
        vec2 = [yy1(:,:);yy1(1,:)];
       vec3 = [tt(:,:);tt(1,:)];
       % --------------FIRST SUBPLOT: POSITION --------------------------
       subplot(2,2,1)
       plot(vec1,vec2,'r','linewidth',3);
       hold on
       plot(vec1(1,:),vec2(1,:),'k.','markersize',30)
       axis equal
       axis(ax)
      
       set(gca,'xtick',[])
       set(gca,'ytick',[])
       set(gca,'xcolor','white')
       set(gca,'ycolor','white')
       titleStr = ['t = ' num2str(time(k),'%4.2e') ...
         ' eA = ' num2str(ea(k),'%4.2e') ...
         ' eL = ' num2str(el(k),'%4.2e')];
       title(titleStr);
       
       % --------------SECOND SUBPLOT CHOICE 1: Concentration -----------
 %       subplot(2,2,2)
 %       h = cline(vec1,vec2,vec3);
 %       set(h,'linewidth',3)
 %       colorbar
 %       axis equal
 %      
 %       set(gca,'xtick',[])
 %       set(gca,'ytick',[])
 %       set(gca,'xcolor','white')
 %       set(gca,'ycolor','white')
 %       title('Concentration');
 %       hold off
 %       
        % ------------SECOND SUBPLOT COICE TWO: Bending Modulus ---------
       subplot(2,2,2)
       N = length(conc(:,:,k));
       rbn = 1 * (ones(N,1) - conc(:,:,k)) + 0.1*conc(:,:,k);
       vec4 = [rbn(:,:);rbn(1,:)];
       h = cline(vec1,vec2,vec4);
       set(h,'linewidth',3)
       colorbar
       hold on
       plot(vec1(1,:),vec2(1,:),'k.','markersize',30)
       
       axis equal
       axis(ax)
      
       set(gca,'xtick',[])
       set(gca,'ytick',[])
       set(gca,'xcolor','white')
       set(gca,'ycolor','white')
       title('Bending Modulus');
       hold off
       
       % --------------THIRD SUBPLOT: Tension ---------------------------
       subplot(2,2,3)
       plot(ten(:,:,k),'linewidth',3)
       xlim([0 N])
       title('\Lambda^{LL}');
       
       % --------------FOURTH SUBPLOT: Curvature squared ----------------
       subplot(2,2,4)
       [~,~,cur] = oc.computeOpeningAngle(N,[posx(:,:,k);posy(:,:,k)]); 
       plot(cur.^2,'linewidth',3)
       xlim([0 N])
       title('\kappa^2');
       
      pause(0.1);
  end
 
