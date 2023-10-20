%clc;clear;close all

addpath ..

% set(gcf,'Position',get(0,'Screensize'));
% set(0,'DefaultAxesFontSize',22)
% options.savefig = false;

%file = 'extensional_RAp8_Conc0p3_Chi10_beta0_n500_dt0en4.bin';
%file = 'parabolic_Conc0p3_Chi600_beta0_n500_dt0en4_h2.bin';

%file = 'Chi800_RA0p5_Conc0p3_Beta0_y0p1_eps0p04_n20.bin';
% name = 'Chi200_RA0p95_Conc0p3_Beta0_y0p1_eps0p04_n20.jpg';
ax = [-2 2 -2 2];

file = 'blah2.bin';
%file = 'relaxation1VesA.bin';
%file = 'longChoke_Chi400_Scale0p49_shortax3p45_conc0p3.bin';
%file = 'longChoke_Chi400_Scale0p59_shortax2p75_conc0p3.bin';
%file = 'longChoke_Chi400_Scale0p71_shortax2p20_conc0p3.bin';
%file = 'longChoke_Chi400_Scale0p85_shortax1p70_conc0p3.bin';
%file = 'longChoke.bin';
%file1 = 'Unconfined_larger_800.bin';
[posx,posy,conc,ea,el,time,xvel1,yvel1,ten] = loadFile(file);
%[posx1,posy1,conc1,ea1,el1,time1,xvel11,yvel11,ten1] = loadFile(file1);
%ten = -ten;

plot(posx(:,:,end),posy(:,:,end))
hold on
%plot(posx1(:,:,end),posy1(:,:,end), 'r--')
% title('Confined vs unconfined for larger ves, x = 800, conc = 0') 
% legend('confined','unconfined')

if 0
    Nbd = 128;
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
% plot(posx(:,:,1),posy(:,:,1))
% hold on 
% plot(xwalls,ywalls)
% pause
% N = numel(posx(:,:,1));
% oc = curve(N);
% CoM = [];
% for k = 1:length(time)
%     C = oc.centerOfMass([posx(:,:,k);posy(:,:,k)]);
%     CoM = [CoM;C];
% end
% plot(time,CoM(:,2),'linewidth',3)
% semilogy(time,abs(CoM(:,2)))
% title('Center of mass: $\chi$ = 400, $RA_0$ = 0.95, $u$ = 0.3','interpreter','latex')
% xlabel('time')
% ylabel('y')
% saveas(gcf,name)
% vec1 = [posx(:,:,end);posx(1,:,end)];
% vec2 = [posy(:,:,end);posy(1,:,end)];
% N = numel(conc(:,:,end));
% rbn = 1 * (ones(N,1) - conc(:,:,end)) + 0.1*conc(:,:,end);
% vec4 = [rbn(:,:);rbn(1,:)];
% h = cline(vec1,vec2,vec4);
% set(h,'linewidth',3)
% axis equal
% 
% set(gca,'xtick',[])
% set(gca,'ytick',[])
% set(gca,'xcolor','white')
% set(gca,'ycolor','white')
% hold off

% 
% % plot(mean(squeeze(posy)))
% % pause
% 
% %plot(time,squeeze(mean(posy)))


 irate = 1000; 
 istart = 1;
 iend = numel(time);
 ntime = iend;
 ntime = numel(time);
 N = length(posx(:,:,1));
 time = time(istart:iend);
 oc = curve(N);
 
  for k = istart:irate:iend
       clf; 
       xx1 = posx(:,:,k);
       yy1 = posy(:,:,k);
       tt = conc(:,:,k);
       [ra(k),area(k),~] = oc.geomProp([xx1;yy1]); 
 % 
       vec1 = [xx1(:,:);xx1(1,:)];
       vec2 = [yy1(:,:);yy1(1,:)];
       vec3 = [tt(:,:);tt(1,:)];
       % --------------FIRST SUBPLOT: POSITION --------------------------
       subplot(2,2,1)
       plot(vec1,vec2,'r','linewidth',3);
       hold on
       plot(vec1(1,:),vec2(1,:),'k.','markersize',20)
       plot([ax(1) ax(2)],[0 0],'k--','linewidth',2)
%       plot(xwalls,ywalls,'k','linewidth',2)
       axis equal
       axis(ax)
      
       set(gca,'xtick',[])
       set(gca,'ytick',[])
       set(gca,'xcolor','white')
       set(gca,'ycolor','white')
       titleStr = ['t = ' num2str(time(k),'%4.2e') ...
         ' eA = ' num2str(ea(k),'%4.2e') ...
         ' eL = ' num2str(el(k),'%4.2e')];
       title(titleStr,'fontsize',20);
       
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
       N = numel(conc(:,:,k));
       rbn = 0.9/2*tanh(3*(conc(:,:,k) - 0.5)) + 1.1/2;
%       rbn = 1 * (ones(N,1) - conc(:,:,k)) + 0.1*conc(:,:,k);
       vec4 = [rbn(:,:);rbn(1,:)];
       h = cline(vec1,vec2,vec4);
       set(h,'linewidth',3)
       h2 = colorbar;
       set(h2,'fontsize',16)
       hold on
       plot(vec1(1,:),vec2(1,:),'k.','markersize',20)
%       plot(xwalls,ywalls,'k','linewidth',2)
    
       axis equal
       axis(ax)
      
       set(gca,'xtick',[])
       set(gca,'ytick',[])
       set(gca,'xcolor','white')
       set(gca,'ycolor','white')
       title('Bending Modulus','fontsize',20);
       hold off
       
       % --------------THIRD SUBPLOT: Tension ---------------------------
       subplot(2,2,3)
       plot(conc(:,:,k),'linewidth',3)
       xlim([0 N])
       title('Concentration','fontsize',20);
%       title('\Lambda^{LL}');
      set(gca,'fontsize',16)
       
       subplot(2,2,4)
       plot(rbn,'linewidth',3)
       xlim([0 N])
       title('Bending Modulus','fontsize',20);
%       title('\Lambda^{LL}');
      set(gca,'fontsize',16)
       
       % --------------FOURTH SUBPLOT: Curvature squared ----------------
%        subplot(2,2,4)
% %       [~,~,cur] = oc.computeOpeningAngle(N,[posx(:,:,k);posy(:,:,k)]); 
% pause
%        [~,~,cur] = oc.diffProp([xx1;yy1]); 
%        plot(cur.^2,'linewidth',3)
%        xlim([0 N])
%        title('\kappa^2');
%        
       pause(0.1);
  end
 
