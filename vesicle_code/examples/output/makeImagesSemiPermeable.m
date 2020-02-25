addpath ../../src
set(0,'DefaultAxesFontSize',22)
options.savefig = false;

irate = 10000; % controls the speed of the visualization
%file = 'relaxation1VesData.bin';
%file = 'shear1VesData.bin';
%file = 'extensional1VesData.bin';
%file = '~/projects/brinkman/vesicle_code/results/semipermeable/shear1VesData_fluxCoeff1_fluxShape3_fluxWidth4.bin';
%file = '~/projects/brinkman/vesicle_code/examples/output/shear1VesData.bin';
file = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e1_ra065_beta2p0e0/shear1VesData.bin';
chi = 5.0; beta = 1.0;

[posx,posy,sigma,wallx,wally,ea,el,time,n,nv] = loadFile(file);
% load positions, "tension", stresses, errors, time, number of points, and
% number of vesicles

op = poten(n);

xmin = min(min(posx));
xmax = max(max(posx));
ymin = min(min(posy));
ymax = max(max(posy));
ax = [xmin xmax ymin ymax];

istart = 1;
iend = numel(time);
ntime = numel(time);

area = zeros(ntime,1);
ra = zeros(ntime,1);
oc = curve;
[~,~,L] = oc.geomProp([posx(:,1,1);posy(:,1,1)]);
tension = zeros(n,ntime);
trac = zeros(2*n,ntime);
normal = zeros(2*n,ntime);
vel = zeros(2*n,ntime);
cur = zeros(n,ntime);
jac = zeros(n,ntime);
flux = zeros(n,ntime);
normal_trac = zeros(n,ntime);
IA = zeros(ntime,1); % inclination angle
for k = 1:iend
  [ra(k),area(k)] = oc.geomProp([posx(:,:,k);posy(:,:,k)]);
  ves = capsules([posx(:,:,k);posy(:,:,k)],sigma(:,k),[],1,1);
  trac(:,k) = ves.tracJump([posx(:,:,k);posy(:,:,k)],sigma(:,k));

  [jac(:,k),tang,cur(:,k)] = oc.diffProp([posx(:,1,k);posy(:,1,k)]);
  normal(:,k) = [tang(end/2+1:end);-tang(1:end/2)];
  tension(:,k) = sigma(:,1,k) + 1.5*cur(:,k).^2;
end

for k = istart:irate:iend

  %calculate the moment of inertia using the formula in section 4.1 of
  %Rahimian, Veerapaneni, Biros 2010
  center = [mean(posx), mean(posy)];
  r = [posx(:,:,k),posy(:,:,k)];
  rdotn = posx(:,:,k).*normal(1:end/2,k) + ...
          posy(:,:,k).*normal(end/2+1:end,k);
  J(1,1) = 2*pi/n*sum(rdotn.*(r(:,2).^2).*jac(:,k));
  J(1,2) = 2*pi/n*sum(rdotn.*(-r(:,1).*r(:,2)).*jac(:,k));
  J(2,1) = J(1,2);
  J(2,2) = 2*pi/n*sum(rdotn.*(r(:,1).^2).*jac(:,k));

  %Use the moment of inertia tensor to calculate the inclination angle
  [V,D] = eig(J);
  [~,I] = min(diag(D));
  V = V(:,I);
  if V(1) || V(2) < 0 
    z = -V(1) + 1i*(-V(2));
  else
    z = V(1) + 1i*(V(2));
  end
  %In radians
  IA(k) = angle(z);

  normal_trac(:,k) = ...
       trac(1:end/2,k) .* normal(1:end/2,k) + ...
       trac(end/2+1:end,k) .* normal(end/2+1:end,k);
  vec3 = [tt;tt(1)];
  G = op.stokesSLmatrix(ves);
  vel(:,k) = chi*[posy(:,:,k);zeros(n,1)] + G*trac(:,k);
  flux(:,k) = vel(1:n,k) .* normal(1:n,k) + ...
              vel(n+1:2*n,k) .* normal(n+1:2*n,k);
end
min_ten = floor(min(min(min(sigma))));
max_ten = floor(max(max(max(sigma))));

figure(1); clf
for k = istart:irate:iend
  clf
  subplot(3,3,1)
  xx = posx(:,:,k);
  yy = posy(:,:,k);
  vec1 = [xx(:,:);xx(1,:)];
  vec2 = [yy(:,:);yy(1,:)];
  plot(vec1,vec2,'r','linewidth',3)
  hold on;
  plot(vec1(1,:),vec2(1,:),'b.','markersize',20)
  hold off
  axis equal
  axis(ax)
  titleStr = ['t = ' num2str(time(k),'%4.2e') ...
    ' eL = ' num2str(el(k),'%4.2e')];
  title(titleStr)

  subplot(3,3,2);hold on
  tt = tension(:,k);
  vec3 = [tt;tt(1)];
  for j = 1:nv
    h = cline(vec1(:,j),vec2(:,j),vec3(:,j));
    set(h,'LineWidth',4);
  end
  hold on
  plot(vec1(1,:),vec2(1,:),'b.','markersize',20)
  hold off
  axis equal
  axis(ax);
  set(gca,'xtick',[]);
  set(gca,'ytick',[]);
  set(gca,'xcolor','w');
  set(gca,'ycolor','w');
  title('Tension')
  colorbar


  subplot(3,3,3);hold on
  tt = flux(:,k);
  vec3 = [tt;tt(1)];
  for j = 1:nv
    h = cline(vec1(:,j),vec2(:,j),vec3(:,j));
    set(h,'LineWidth',4);
  end
  hold on
  plot(vec1(1,:),vec2(1,:),'b.','markersize',20)
  hold off
  axis equal
  axis(ax);
  set(gca,'xtick',[]);
  set(gca,'ytick',[]);
  set(gca,'xcolor','w');
  set(gca,'ycolor','w');
  title('Normal Flux')
  colorbar


  subplot(3,3,4);hold on
  tt = trac(1:end/2,k);
  vec3 = [tt;tt(1)];
  for j = 1:nv
    h = cline(vec1(:,j),vec2(:,j),vec3(:,j));
    set(h,'LineWidth',4);
  end
  hold on
  plot(vec1(1,:),vec2(1,:),'b.','markersize',20)
  hold off
  axis equal
  axis(ax);
  set(gca,'xtick',[]);
  set(gca,'ytick',[]);
  set(gca,'xcolor','w');
  set(gca,'ycolor','w');
  title('Traction x')
  colorbar

  subplot(3,3,5);hold on
  tt = trac(end/2+1:end,k);
  vec3 = [tt;tt(1)];
  for j = 1:nv
    h = cline(vec1(:,j),vec2(:,j),vec3(:,j));
    set(h,'LineWidth',4);
  end
  hold on
  plot(vec1(1,:),vec2(1,:),'b.','markersize',20)
  hold off
  axis equal
  axis(ax);
  set(gca,'xtick',[]);
  set(gca,'ytick',[]);
  set(gca,'xcolor','w');
  set(gca,'ycolor','w');
  title('Traction y')
  colorbar

  subplot(3,3,6);hold on
  vec3 = [normal_trac(:,k);normal_trac(1,k)];
  for j = 1:nv
    h = cline(vec1(:,j),vec2(:,j),vec3(:,j));
    set(h,'LineWidth',4);
  end
  hold on
  plot(vec1(1,:),vec2(1,:),'b.','markersize',20)
  hold off
  axis equal
  axis(ax);
  set(gca,'xtick',[]);
  set(gca,'ytick',[]);
  set(gca,'xcolor','w');
  set(gca,'ycolor','w');
  title('Normal Traction')
  colorbar

  subplot(3,3,7);
  plot(time(1:k),area(1:k))
  xlim([0 time(end)])
  ylim([min(area) max(area)])
  title('Area')
  xlabel('time')

  subplot(3,3,8)
  plot(time(1:k),squeeze(mean(posx(:,:,1:k))))
  xlim([0 time(end)])
  ylim([min(min(mean(posx))) max(max(mean(posx)))])
  title('Mean x Position')
  xlabel('time')

  subplot(3,3,9)
  plot(time(1:k),squeeze(mean(posy(:,:,1:k))))
  xlim([0 time(end)])
  ylim([min(min(mean(posy))) max(max(mean(posy)))])
  title('Mean y Position')
  xlabel('time')

  pause(0.01)

end

%figure(2); clf; hold on;
%plot(time(1:end-1),diff(area)./diff(time));
%
%dArea = zeros(ntime,1);
%%dArea2 = zeros(ntime,1);
%for k = 1:ntime
%  dArea(k) = beta*2*pi/n*sum(normal_trac(:,k).*jac(:,k));
%%  dArea2(k) = -beta*2*pi/n*sum((cur(:,k).^3 + ...
%%      cur(:,k).*sigma(:,k)).*jac(:,k));
%  % not sure why the negative sign is necessary based on the writeup.
%  % Probably an IBP error
%end
%plot(time,dArea,'r--')
%%plot(time,dArea2,'k--')


