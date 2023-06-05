if 0
load('Stenosis_RAp4_SC_bulk_targets.mat');
load('Stenosis_RAp4_SC_bulk.mat');
load('Stenosis_RAp4_SC_bulk_pressure.mat');
load('walls.dat');
end

xx = xtar; yy = ytar;
uu = reshape(uxtar,139,509);
vv = reshape(uytar,139,509);
xwalls = XWalls(1:end/2);
ywalls = XWalls(end/2+1:end);

%uu = uu - mean(uxvel);

posx = squeeze(posx1(:,:,45001));
posy = squeeze(posy1(:,:,45001));

xmin = 2.0;
xmax = 2.5;
% index of points in a given x-range
s = find((xx(1,:) > xmin) & (xx(1,:) < xmax));

figure(1); clf; hold on
plot(posx,posy,'r','linewidth',2);
plot(xwalls,ywalls,'k','linewidth',2);
step = 1;
quiver(xx(1:step:end,1:step:end),yy(1:step:end,1:step:end),...
       uu(1:step:end,1:step:end),vv(1:step:end,1:step:end),'b')

%axis([min(min(xx)) max(max(xx)) min(min(yy)) max(max(yy))])
axis equal;
axis([xmin xmax -0.8 -0.3])

posx_neg = posx(posy < -0.4);
posy_neg = posy(posy < -0.4);
%k = 1;
figure(2); clf; hold on; figure(1);
for k = 1:1
%for k = 1:numel(s)
  [~,kmin] = min(abs(posx_neg - xx(1,s(k))));
  s2 = find(yy(:,1) < posy_neg(kmin));
  plot(xx(s2,s(k)),yy(s2,s(k)),'g.')

  figure(2);
%  plot(yy(s2,s(k)),uu(s2,s(k)))
  yyy = yy(s2,s(k)) - yy(s2(1),s(k));
  h = yyy(end);

  chi = 1/h;
%  vel = chi*uu(s2,s(k));
  vel = uu(s2,s(k));
%  vel = vel/vel(end);

%  plot(yyy/yyy(end),chi*uu(s2,s(k)))
  plot(yyy/h,(h*vel)/(h*vel(end)))
  figure(1);
end

z = linspace(0,1,100);
figure(2); hold on;
%plot(z,z.^2,'k--','linewidth',2)
%plot(z,z,'k--','linewidth',2)
plot(z,-z.*(z-2),'k--','linewidth',2)


modes = (-512:511)';
modes(abs(modes) > 50) = 0;
sig = squeeze(ten1(:,:,45001));
sigh = fftshift(fft(sig));
Dsigh = 1i*modes.*sigh;
Dsig = ifft(ifftshift(Dsigh));
figure(6); clf;
cline(posx,posy,Dsig)
hold on;
plot3(posx(1),posy(1),1000,'k.','markersize',20)
%plot3(posx(30),posy(30),1000,'r.','markersize',20)
axis equal
colorbar

figure(5);
plot(sig);

