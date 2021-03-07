addpath ..
set(gcf,'Position',get(0,'Screensize'));
set(0,'DefaultAxesFontSize',22)
options.savefig = false;

%file = 'Chi10_ra065_beta0p1_conc0p5.bin';
%file = 'Chi10_ra065_beta0_conc0p5.bin';
%file = 'Chi20_ra065_beta10_conc0p3.bin';
%file = 'Chi0_ra065_beta10_conc0p5.bin';
%file = 'Chi4_ra07564_beta0_conc0p48_MoreEL.bin';
%file = 'Chi4_ra07564_beta0_conc0_MoreEL.bin';
file = 'Chi0_ra07564_beta0_conc0p3_dtp25em5.bin';
[posx,posy,conc,ea,el,time,xvel1,yvel1,ten] = loadFile(file);

irate = 1; 
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
%plot(time, area)