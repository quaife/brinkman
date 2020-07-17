addpath ../../src
set(0,'DefaultAxesFontSize',22)
options.savefig = false;
options.pressure = false;

irate = 1; % controls the speed of the visualization

if 0
  file = 'parabolic1Ves2Data.bin';
%  file = '~/projects/brinkman/vesicle_code/results/parabolic_offcenter/pflowR10u1e2B1em4Data.bin';
%  file = '~/projects/brinkman/vesicle_code/results/parabolic_offcenter/pflowR10u1e0p5B1em4bData.bin';
%  file = '~/projects/brinkman/vesicle_code/results/parabolic_offcenter/pflowR10u1em0p0B1em4bData.bin';
%  file = '~/projects/brinkman/vesicle_code/results/parabolic_offcenter/pflowR10u1em0p5B1em4bData.bin';
%  file = '~/projects/brinkman/vesicle_code/results/parabolic_offcenter/pflowR10u1em1p0B1em4bData.bin';
%  file = '~/projects/brinkman/vesicle_code/results/parabolic_offcenter/pflowR10u1ep0B1em4bData.bin';
%  file = '~/projects/brinkman/vesicle_code/results/parabolic_offcenter/pflowR10u1ep1p0B1em4dData.bin';
  file = '~/projects/brinkman/vesicle_code/results/parabolic_offcenter/pflowR10u1ep1p0B0em4dData.bin';

  ax = [-3 3 -3 3];
  irate = 10;
  options.confined = false;
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

if 0
  file = 'relaxation1VesData.bin';
%  file = '~/projects/brinkman/vesicle_code/results/Apr142020/starBeta1em3/relaxation1VesData.bin';
%    file = '~/projects/brinkman/vesicle_code/results/Apr232020/ellipseBeta1em5/relaxation1VesData.bin';
%  ax = [-3 3 -3.5 3.5];
  ax = 2*[-1 1 -1 1];
  options.confined = false;
  options.savefig = false;
  beta = 1e-5;
  count = 1;
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
%  file = 'shear1VesDData_Part4.bin';
  file = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e1_ra065_beta1p0em5/shear1VesData_Part5.bin';
  ax = [-5 5 -5 5];
  options.confined = false;
  beta = 0.2;
end
if 0
%  file = '~/projects/brinkman/vesicle_code/results/shear2Ves/adR1em1adS1e0Chi5em1_ra090/shear2VesData.bin';
  file = '~/projects/brinkman/vesicle_code/results/shear2Ves/adR1em1adS3em1Chi2p5em1_ra090/shear2VesData_Part2.bin';
  ax = [-3 3 -3 3];
  options.confined = false;
end

if 1
  file = 'chokeMulti1VesData.bin';
  ax = [-1.3 11.3 -0.1 5.6];
  options.confined = true;
  options.pressure = false;
  options.savefig = false;
  count = 1;
  irate = 1;
end

if 0
  file = 'choke1VesData.bin';
%file = '~/projects/brinkman/vesicle_code/results/choke1VesLong/beta0Scale1p44_kappa1e0_farfield1e0/choke1VesData.bin'
%file = '~/projects/brinkman/vesicle_code/results/choke1VesLong/beta0Scale1p44_kappa1e1_farfield1e0/choke1VesData.bin';
%file = '~/projects/brinkman/vesicle_code/results/choke1VesLong/beta0Scale1p44_kappa1e2_farfield1e0/choke1VesData.bin';
%file = '~/projects/brinkman/vesicle_code/results/choke1VesLong/beta0Scale1p44_kappa1e0_farfield5e2/choke1VesData.bin';
%file = '~/projects/brinkman/vesicle_code/results/choke1VesLong/beta0Scale1p44_kappa1e0_farfield1e2/choke1VesData.bin';
file = '~/projects/brinkman/vesicle_code/results/choke1VesLong/beta0Scale1p44_kappa1e0_farfield5e2_offcenter/choke1VesData.bin';
%file = '~/projects/brinkman/vesicle_code/results/choke1VesLong/beta1em3Scale1p44_kappa1e0_farfield5e2_offcenter/choke1VesData.bin';
%file = '~/projects/brinkman/vesicle_code/results/choke1VesLong/beta1em4Scale1p44_kappa1e0_farfield5e2_offcenter/choke1VesData.bin';
%file = '~/projects/brinkman/vesicle_code/results/choke1VesLong/beta1em5Scale1p44_kappa1e0_farfield5e2_offcenter/choke1VesData.bin';
  ax = [-50 50 -12.5 12.5];
  options.confined = true;
  options.pressure = true;
  options.savefig = false;
  count = 1;
  irate = 10;
end
if 0
%  file = 'slit1VesData.bin';
%  file = '~/projects/brinkman/vesicle_code/results/slit/Scale0p55_Beta0/slit1VesData.bin';
  file = '~/projects/brinkman/vesicle_code/results/slit/Scale0p55_Beta1em3/slit1VesData.bin';
%  file = '~/projects/brinkman/vesicle_code/results/slit/Scale0p55_Beta1em4/slit1VesData.bin';
%  file = '~/projects/brinkman/vesicle_code/results/slit/Scale0p55_Beta1em5/slit1VesData.bin';
  ax = [-5 5 -2.5 2.5];
  options.confined = true;
  options.pressure = true;
  options.savefig = false;
  count = 1;
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
% load positions, tension, errors, time, number of points, and
% number of vesicles
if options.pressure
  fileName = [file(1:end-8) 'Pressure.bin'];
  [pressx,pressy,press] = loadPressure(fileName);
  dpress = diff(press);
end

istart = 1;
iend = numel(time);
ntime = iend;
ntime = numel(time);
%time = time(istart:iend);

oc = curve;
[~,~,L] = oc.geomProp([posx(:,1,1);posy(:,1,1)]);
area = zeros(ntime,1);
length = zeros(ntime,1);
ra = zeros(ntime,1);
incAng = zeros(ntime,1);
trac = zeros(2*n,nv,ntime);
normal = zeros(2*n,nv,ntime);
flux = zeros(n,ntime);
cur = zeros(n,nv,ntime);
jac = zeros(n,nv,ntime);
sa = zeros(n,nv,ntime);
velBen = zeros(2*n,ntime);
velTen = zeros(2*n,ntime);
velSLP = zeros(2*n,ntime);
velFlux = zeros(2*n,ntime);
cx = zeros(ntime,1); % x coordinate of center of mass
cy = zeros(ntime,1); % y coordinate of center of mass

op = poten(n);
for k = istart:1:iend
  [ra(k),area(k),length(k)] = oc.geomProp([posx(:,1,k);posy(:,1,k)]);
%  incAng(k) = InclinationAngle(posx(:,1,k),posy(:,1,k));

  if 1
    vesicle = capsules([posx(:,:,k);posy(:,:,k)],[],[],1,1);
    normal(1:end/2,:,k) = +vesicle.xt(end/2+1:end,:);
    normal(end/2+1:end,:,k) = -vesicle.xt(1:end/2,:);
%    trac(:,:,k) = vesicle.tracJump([posx(:,:,k);posy(:,:,k)],ten(:,:,k));
    [jac(:,:,k),~,cur(:,:,k)] = oc.diffProp([posx(:,1,k);posy(:,1,k)]);
%    ten(:,:,k) = ten(:,:,k) + 1.5*cur(:,:,k).^2;
    cx(k) = 1/(2*area(k))*sum(...
        posx(:,:,k).^2.*normal(1:end/2,:,k).*jac(:,:,k))*2*pi/n;
    cy(k) = 1/(2*area(k))*sum(...
        posy(:,:,k).^2.*normal(end/2+1:end,:,k).*jac(:,:,k))*2*pi/n;
  end
%
%  G = op.stokesSLmatrix(vesicle);
%  velSLP(:,k) = G*trac(:,:,k);
%
%  flux(:,k) = -(trac(1:end/2,:,k).*normal(1:end/2,:,k) + ...
%          trac(end/2+1:end,:,k).*normal(end/2+1:end,:,k));
%  velFlux(:,k) = [flux(:,k).*normal(1:end/2,:,k); ...
%                  flux(:,k).*normal(end/2+1:end,:,k)];
end

if 0
dArea = zeros(ntime,1);
dArea2 = zeros(ntime,1);
for k = istart:1:iend
  dArea(k) = beta*2*pi/n*sum(-flux(:,k).*jac(:,k));
  dArea2(k) = -beta*2*pi/n*sum((-cur(:,k).^3/2 + ...
      cur(:,k).*ten(:,k)).*jac(:,k));
  % not sure why the negative sign is necessary based on the writeup.
  % Probably an IBP error
end
end


%min_ten = floor(min(min(min(ten))));
%max_ten = floor(max(max(max(ten))));
min_ten = -15;
max_ten = +15;
min_flux = -15;
max_flux = +2;

figure(1); clf
for k = istart:irate:iend
%  xx = interpft(posx(:,:,k),256); yy = interpft(posy(:,:,k),256);  
  xx = posx(:,:,k) - 0*cx(k);
  yy = posy(:,:,k);
  tt = ten(:,:,k);
  vec1 = [xx(:,:);xx(1,:)];
  vec2 = [yy(:,:);yy(1,:)];
  vec3 = [tt(:,:);tt(1,:)];
  vec4 = [flux(:,k);flux(1,k)];
  if 1
    clf; hold on;
    plot(vec1,vec2,'r-','linewidth',3)
%    plot(0,cy(k),'k.','markersize',10);
%    plot([-5 5],[0 0],'k--')
%    plot(vec1(1,:),vec2(1,:),'b.','markersize',20)
%    for j = 1:1
%      subplot(1,2,1)
%      h = cline(vec1(:,j),vec2(:,j),vec3(:,j));
%      set(h,'linewidth',3)
%    end
%    colorbar
%    caxis([min_ten max_ten])
%    axis equal
%    axis(ax)
%    for j = 1:1
%      subplot(1,2,2)
%      h = cline(vec1(:,j),vec2(:,j),vec4(:,j));
%      set(h,'linewidth',3)
%    end
%    colorbar
%    caxis([min_flux max_flux])
    if options.confined
      vec1 = [wallx(:,:);wallx(1,:)];
      vec2 = [wally(:,:);wally(1,:)];
      plot(vec1,vec2,'k','linewidth',3)
    end
    hold off
    axis equal
    axis(ax)
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    set(gca,'xcolor','white')
    set(gca,'ycolor','white')
    titleStr = ['t = ' num2str(time(k),'%4.2e') ...
      ' eA = ' num2str(ea(k),'%4.2e') ...
      ' eL = ' num2str(el(k),'%4.2e')];
    title(titleStr)
    pause(0.01)
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
    ti = get(gca,'TightInset');
    set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
    set(gca,'units','centimeters')
    pos = get(gca,'Position');
    ti = get(gca,'TightInset');
    set(gcf, 'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);

    print(gcf,'-dpdf','-r300',filename);
  end
  pause(0.01)
%  pause
end





