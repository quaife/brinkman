function [htop,hbot,maxtop,maxbot,ptop,pbot,indPos,indNeg,spos,sneg] = ...
        BLvelocities3(icol,fileName)
%icol = 279;
iplot = ~true;

%load bulkVel_RAp4SCatT45.mat;
%fileName = 'AshleyFiles/Stenosis_RAp4_SCp55_pos4.mat';
load(fileName)

xx = xtar;
yy = ytar;
uu = uxtar;
vv = uytar;

% distance the vesicle points are 'ahead' or 'behind' the cross section
z = posx - xx(1,icol);
z = [z;z(1)];
% crosses between index s and index s+1
s = find(z(2:end)./z(1:end-1) < 0);

[~,spos] = max(posy(s));
[~,sneg] = min(posy(s));
if numel(s) == 0
  indPos = (1:size(xtar,1))';
  indNeg = (1:size(xtar,1))';
else
  spos = s(spos);
  sneg = s(sneg);

  indPos = find(yy(:,icol) > posy(spos));
  indNeg = find(yy(:,icol) < posy(sneg));
end

if iplot
  figure(1); clf; hold on
  plot(posx,posy,'r','linewidth',2);
  axis equal;
%  plot(xx(:,icol),yy(:,icol),'k.');
  plot(xx(indNeg,icol),yy(indNeg,icol),'k.');
  plot(xx(indPos,icol),yy(indPos,icol),'k.');
  plot(xx(1,:),-0.7*ones(size(xx(1,:))),'k','linewidth',2);
  plot(xx(1,:),+0.7*ones(size(xx(1,:))),'k','linewidth',2);
  quiver(xx(indNeg,icol),yy(indNeg,icol),uu(indNeg,icol),vv(indNeg,icol),'b')
  quiver(xx(indPos,icol),yy(indPos,icol),uu(indPos,icol),vv(indPos,icol),'b')
  axis([xx(1,icol) - 0.7 xx(1,icol) + 0.7 -0.71 0.71])
end




if iplot
  figure(2); clf; hold on;
  plot(yy(:,icol),uu(:,icol))
  plot(yy(indNeg,icol),uu(indNeg,icol),'k','linewidth',2)
  plot(yy(indPos,icol),uu(indPos,icol),'k','linewidth',2)
end

zbot = yy(indNeg,icol) + 0.7;
pbot = polyfit(zbot,uu(indNeg,icol),2);
if iplot
  figure(3); clf; hold on
  plot(zbot,uu(indNeg,icol),'linewidth',2)
  xlabel('Distance from bottom wall')
  ylabel('x component of velocity');
  plot(zbot,polyval(pbot,zbot),'k--','linewidth',2)
  h = legend('Velocity','Quadratic Fit');
  set(h,'location','northwest')
end

ztop = 0.7 - yy(indPos,icol);
ptop = polyfit(ztop,uu(indPos,icol),2);
if iplot
  figure(4); clf; hold on
  plot(ztop,uu(indPos,icol),'linewidth',2)
  xlabel('Distance from top wall')
  ylabel('x component of velocity');
  plot(ztop,polyval(ptop,ztop),'k--','linewidth',2)
  h = legend('Velocity','Quadratic Fit');
  set(h,'location','northwest')
end

htop = max(ztop);
hbot = max(zbot);
maxtop = max(uu(indPos,icol));
maxbot = max(uu(indNeg,icol));

