addpath ../../src
set(0,'DefaultAxesFontSize',22)
options.savefig = false;

irate = 5; % controls the speed of the visualization

if 0
  file = 'parabolic1VesData.bin';
end
if 0
%  file = 'extensional1VesCleanData.bin'; irate = 2;
  file = 'extensional1VesSemiData.bin'; irate = 2;
%  file = 'extensional1VesCircleCleanData.bin'; irate = 2;
%  file = 'extensional1VesCircleSemiData.bin'; irate = 1;

%  file = 'extensional1VesData.bin'; 


%file = '~/projects/brinkman/vesicle_code/docs/adhesion/makefigs/extensional_adR4em1adS7em1Chi1em1_ra070/extensional2VesData.bin';
%  file = '~/projects/brinkman/vesicle_code/results/extensional2Ves/adR4em1adS7em1Chi1em2_ra070/extensional2VesData.bin';
%  file = '~/projects/brinkman/vesicle_code/results/yuan_runs/june_28_2019/extensional2VesData75be.bin';
  ax = 1*[-3 3 -3 3];
  options.confined = false;
end
if 0
  file = 'extensionalManyVesData.bin';
  ax = [-12 12 -2 2];
  options.confined = false;
end
if 0
  file = '1Vesicle_Shear2_beta1_RApt65.bin';
  ax = (2*[-3 3 -3 3]);
  options.confined = false;
end

if 1
  file = 'relaxation1VesData.bin';
  ax = [-3 3 -3.5 3.5];
  options.confined = false;
end

if 0
%  file = '~/projects/brinkman/vesicle_code/results/May092019/relaxation1Ves/kappa1em1_beta1ep0_ra0p95/relaxation1VesData.bin';
  file = '~/projects/brinkman/vesicle_code/results/relaxationManyVes/Segment1/relaxationManyVesData.bin';
  ax = [-6 6 -6 6];
  options.confined = false;
  options.savefig = false;
  count = 1;
end
if 0
  file = 'relaxation2VesData.bin';
%  file = '~/presentations/2018/lifeSciences2018/results/relaxation/RA65_Range8_Strength2/relaxation2VesData.bin';
  ax = [-4 4 -3 3];
  options.confined = false;
end
if 0
  file = 'relaxation4VesData.bin';
  ax = [-1 1 -1 1];
  options.confined = false;
end
if 0
  file = 'nshearVes03n.Data.bin';
  ax = [-2 2 -2 2];
  options.confined = false;
end
if 0
%  file = 'shear1VesData.bin';
  file = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi8p0em1_ra065_beta1p0em1/shear1VesData.bin';
  ax = [-5 5 -5 5];
  options.confined = false;
end
if 0
%  file = '~/projects/brinkman/vesicle_code/results/shear2Ves/adR1em1adS1e0Chi5em1_ra090/shear2VesData.bin';
  file = '~/projects/brinkman/vesicle_code/results/shear2Ves/adR1em1adS3em1Chi2p5em1_ra090/shear2VesData.bin';
  ax = [-3 3 -3 3];
  options.confined = false;
end
if 0
  file = 'choke1VesData.bin';
  ax = [-10.5 10.5 -3.5 3.5];
  options.confined = true;
end

if 0
  file = 'shear1VesData_ASP_FC0FS1_RApt6_H3_Rpt6.bin';
  options.confined = false;
  ax = [-3 3 -2 2];
end
if 0
  file = 'shear1VesData_ASP_FCpt1FS1_RApt6_H3_Rpt6.bin';
  options.confined = false;
  ax = [-3 3 -2 2];
end
if 0
  file = 'shear1VesData_ASP_FCpt1FS1_RApt6_H5_Rpt6.bin';
  options.confined = false;
  ax = [-3 3 -2 2];
end


[posx,posy,ten,wallx,wally,ea,el,time,n,nv] = loadFile(file);
%[posx,posy,ten,~,~,wallx,wally,ea,el,time,n,nv] = loadFileOld(file);
% load positions, tension, stresses, errors, time, number of points, and
% number of vesicles
istart = 1;
iend = numel(time);
ntime = numel(time);

oc = curve;
[~,~,L] = oc.geomProp([posx(:,1,1);posy(:,1,1)]);
area = zeros(size(time));
ra = zeros(size(time));
%incAng = zeros(size(time));
trac = zeros(2*n,nv,ntime);
cur = zeros(n,nv,ntime);
sa = zeros(n,nv,ntime);
for k = 1:numel(time)
%for k = numel(time)-11:numel(time)
  [ra(k),area(k)] = oc.geomProp([posx(:,1,k);posy(:,1,k)]);
%  incAng(k) = InclinationAngle(posx(:,1,k),posy(:,1,k));

  vesicle = capsules([posx(:,:,k);posy(:,:,k)],[],[],1,1);
%  trac(:,:,k) = vesicle.tracJump([posx(:,:,k);posy(:,:,k)],0*ten(:,:,k));
  [sa(:,:,k),~,cur(:,:,k)] = oc.diffProp([posx(:,1,k);posy(:,1,k)]);
  ten(:,:,k) = ten(:,:,k) + 1.5*cur(:,:,k).^2;
end


%min_ten = floor(min(min(min(ten))));
%max_ten = floor(max(max(max(ten))));
min_ten = -0;
max_ten = +0.2;

figure(1); clf
for k = istart:irate:iend
%  xx = interpft(posx(:,:,k),256); yy = interpft(posy(:,:,k),256);  
  xx = posx(:,:,k);
  yy = posy(:,:,k);
  tt = ten(:,:,k);
  vec1 = [xx(:,:);xx(1,:)];
  vec2 = [yy(:,:);yy(1,:)];
  vec3 = [tt(:,:);tt(1,:)];
  if 1
    clf; hold on;
%    plot(vec1,vec2,'r','linewidth',3)
%%    plot(vec1(1,:),vec2(1,:),'b.','markersize',20)
    for j = 1:1
      h = cline(vec1(:,j),vec2(:,j),vec3(:,j));
      set(h,'linewidth',3)
    end
    colorbar
%    pause
    if options.confined
      vec1 = [wallx(:,:);wallx(1,:)];
      vec2 = [wally(:,:);wally(1,:)];
      plot(vec1,vec2,'k','linewidth',3)
    end
    hold off
    axis equal
    axis(ax)
    titleStr = ['t = ' num2str(time(k),'%4.2e') ...
      ' eA = ' num2str(ea(k),'%4.2e') ...
      ' eL = ' num2str(el(k),'%4.2e')];
    title(titleStr)
  end
%  if 0
%  tt = interpft(ten(:,:,k),96);
%  ss = interpft(shearStress(:,:,k),96);
%  ns = interpft(normalStress(:,:,k),96);
%  clf;
%  subplot(1,3,1);hold on
%  vec3 = [tt(:,:);tt(1,:)];
%  for j = 1:nv
%    h = cline(vec1(:,j),vec2(:,j),vec3(:,j));
%    set(h,'LineWidth',4);
%  end
%  axis equal
%  axis(ax);
%  set(gca,'xtick',[]);
%  set(gca,'ytick',[]);
%  set(gca,'xcolor','w');
%  set(gca,'ycolor','w');
%  title('Tension')
%  colorbar
  caxis([min_ten max_ten])
%
%  subplot(1,3,2);hold on
%  vec3 = [ss(:,:);ss(1,:)];
%  for j = 1:nv
%    h = cline(vec1(:,j),vec2(:,j),vec3(:,j));
%    set(h,'LineWidth',4);
%  end
%  axis equal
%  axis(ax);
%  set(gca,'xtick',[]);
%  set(gca,'ytick',[]);
%  set(gca,'xcolor','w');
%  set(gca,'ycolor','w');
%  title('Shear Stress')
%  colorbar
%  caxis([min_ss max_ss])
%
%  subplot(1,3,3);hold on
%  vec3 = [ns(:,:);ns(1,:)];
%  for j = 1:nv
%    h = cline(vec1(:,j),vec2(:,j),vec3(:,j));
%    set(h,'LineWidth',4);
%  end
%  axis equal
%  axis(ax);
%  set(gca,'xtick',[]);
%  set(gca,'ytick',[]);
%  set(gca,'xcolor','w');
%  set(gca,'ycolor','w');
%  title('Normal Stress')
%  colorbar
%  caxis([min_ns max_ns])
%  
%  titleStr = ['t = ' num2str(time(k),'%4.2e') ...
%      ' eA = ' num2str(ea(k),'%4.2e') ...
%      ' eL = ' num2str(el(k),'%4.2e')];
%  suptitle(titleStr)
%  end
  if options.savefig
    filename = ['./frames/image', sprintf('%04d',count),'.pdf'];
    count = count+1;
    figure(1);
    print(gcf,'-dpdf','-r300',filename);
  end
  pause(0.01)
%  pause
end



