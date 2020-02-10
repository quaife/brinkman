addpath ../../src
set(0,'DefaultAxesFontSize',22)
options.savefig = false;

irate = 1; % controls the speed of the visualization
%file = 'relaxation1VesData.bin';
%file = 'shear1VesData.bin';
file = 'extensional1VesData.bin';
%file = '~/projects/brinkman/vesicle_code/results/semipermeable/shear1VesData_fluxCoeff1_fluxShape3_fluxWidth4.bin';
%file = '~/projects/brinkman/vesicle_code/examples/output/shear1VesData.bin';

[posx,posy,sigma,wallx,wally,ea,el,time,n,nv] = loadFile(file);
% load positions, "tension", stresses, errors, time, number of points, and
% number of vesicles

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
cur = zeros(n,ntime);
jac = zeros(n,ntime);
for k = 1:numel(time)
  ves = capsules([posx(:,:,k);posy(:,:,k)],sigma(:,k),[],1,1);
  trac(:,k) = ves.tracJump([posx(:,:,k);posy(:,:,k)],sigma(:,k));

  [jac(:,k),tang,cur(:,k)] = oc.diffProp([posx(:,1,k);posy(:,1,k)]);
  normal(:,k) = [tang(end/2+1:end);-tang(1:end/2)];
  tension(:,k) = sigma(:,1,k) + 1.5*L*cur(:,k).^2;
  [ra(k),area(k)] = oc.geomProp([posx(:,:,k);posy(:,:,k)]);
end
min_ten = floor(min(min(min(sigma))));
max_ten = floor(max(max(max(sigma))));

figure(1); clf
for k = istart:irate:iend
  clf
  subplot(3,3,1:3)
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

  subplot(3,3,4:6);hold on
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
%  caxis([min_ten max_ten])


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
%dArea1 = zeros(ntime,1);
%dArea2 = zeros(ntime,1);
%for k = 1:ntime
%  dArea1(k) = 2*pi/n*sum((...
%      trac(1:n,k).*normal(1:n,k) + ...
%      trac(n+1:2*n,k).*normal(n+1:2*n,k)).*jac(:,k));
%
%  dArea2(k) = -2*pi/n*sum((cur(:,k).^3 + ...
%      cur(:,k).*sigma(:,k)).*jac(:,k));
%  % not sure why the negative sign is necessary based on the writeup.
%  % Probably a IBP error
%end
%plot(time,dArea1,'r--')
%plot(time,dArea2,'k--')


