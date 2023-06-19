fileName = 'AshleyFiles/Stenosis_RAp4_MCp5_pos4.mat';
load(fileName);

N = numel(ten);
modes = [(0:N/2-1) (-N/2:-1)]';

addpath ~/projects/brinkman/shuwang_rewrite
params.bendsti = 1;
params.bendratio = 1e-1;
params.viscosityInside = 1;
params.viscosityOutside = 1;
params.SPcoeff = 0;
ves = capsules([posx;posy],conc,params);
oc = curve(N);

% subtract off the phase pieces from the force which are part of the
% tension in the code, but not in the paper
ten = ten + oc.CHforce(ves,0.04,100);

RHSflux = oc.fluxj(ves,0.04,100);
rcon_s = oc.diffFT(ves.rcon)/ves.L;

dten = oc.diffFT(ten)/ves.L - rcon_s.*RHSflux;


%figure(5); clf;
%plot(dten)
%pause




istart = 1;
iend = size(xtar,2);
%iend = 465;
ratiotop = [];
ratiobot = [];
lineartop = [];
linearbot = [];
quadtop = [];
quadbot = [];
%K = [145 190 226 237 255 279 300];
K = istart:iend;
htop = [];
hbot = [];
indtop = cell(1,numel(K));
indbot = cell(1,numel(K));
spos = cell(1,numel(K));
sneg = cell(1,numel(K));
for icol = K
%  figure(5); hold on
  [htop(icol),hbot(icol),maxtop,maxbot,ptop,pbot,...
        indtop{icol},indbot{icol},spos{icol},sneg{icol}] = BLvelocities3(icol,fileName);
%  pause
  lineartop = [lineartop ptop(2)];
  linearbot = [linearbot pbot(2)];
  quadtop = [quadtop ptop(1)];
  quadbot = [quadbot pbot(1)];
  ratiotop = [ratiotop ptop(1)/ptop(2)];
  ratiobot = [ratiobot pbot(1)/pbot(2)];
%  plot(htop,maxtop,'ro','markersize',5)
%  plot(hbot,maxbot,'ko','markersize',5)
%  plot(maxtop,htop,'ro','markersize',5)
%  plot(maxbot,hbot,'ko','markersize',5)
%  pause(1e-1)
end

load(fileName);
figure(1); clf; hold on;
%plot(xtar(:,istart:iend),ytar(:,istart:iend),'b.');
plot(wallsx,wallsy,'k','linewidth',2);
quiver(xtar(:,K),ytar(:,K),uxtar(:,K),uytar(:,K),'b')
%plot(posx,posy,'r')
%h = fill(posx,posy,'r');
%set(h,'linestyle','none')
h = cline(posx,posy,conc);
set(h,'linewidth',2)
colorbar
axis equal;
axis([min(posx) - 0.1 max(posx) + 0.1 -0.71 0.71])

figure(2); clf; 
subplot(3,1,1);
plot(xtar(1,istart:iend),lineartop)
title('Top Boundary Layer','fontsize',20)
ylabel('Linear term')
subplot(3,1,2);
plot(xtar(1,istart:iend),quadtop)
ylabel('Quadratic term')
subplot(3,1,3);
plot(xtar(1,istart:iend),ratiotop)
ylabel('Quadratic/Linear term')
xlabel('x')

figure(3); clf; 
subplot(3,1,1);
plot(xtar(1,istart:iend),linearbot)
title('Bottom Boundary Layer','fontsize',20)
ylabel('Linear term')
subplot(3,1,2);
plot(xtar(1,istart:iend),quadbot)
ylabel('Quadratic term')
subplot(3,1,3);
plot(xtar(1,istart:iend),ratiobot)
ylabel('Quadratic/Linear term')
xlabel('x')

figure(5); clf; hold on
surf(xtar,ytar,log10(uxtar.^2 + uytar.^2))
view(2); shading interp; axis equal;
colorbar
h = fill3(posx,posy,100*ones(size(posx)),'r');
set(h,'linestyle','none')
axis equal;
axis([min(posx) - 0.1 max(posx) + 0.1 -0.71 0.71])

figure(6); clf; hold on
h = cline(posx,posy,ten);
set(h,'linewidth',2)
plot(wallsx,wallsy,'k','linewidth',2);
axis equal;
axis([min(posx) - 0.1 max(posx) + 0.1 -0.71 0.71])
colorbar
pause

if 1
% find indicies where quadratic term changes sign and when its
% derivative changes sign (max and mins)
s1 = find(quadtop(2:end).*quadtop(1:end-1) < 0);
dquadtop = smooth(diff(quadtop))';
% points where derivative changes sign (ie. max or min)
s2 = find(dquadtop(2:end).*dquadtop(1:end-1) < 0);

s = [s1 s2];
for j = 1:numel(s)
  disp(s(j))
  figure(1); clf; hold on;
  %plot(xtar(:,istart:iend),ytar(:,istart:iend),'b.');
  plot(wallsx,wallsy,'k','linewidth',2);
  quiver(xtar(:,K),ytar(:,K),uxtar(:,K),uytar(:,K),'b')
%  h = fill(posx,posy,'r');
%  set(h,'linestyle','none')
  h = cline(posx,posy,conc);
  colorbar
  set(h,'linewidth',2)
  plot(xtar(indtop{s(j)},s(j)),ytar(indtop{s(j)},s(j)),'k--','linewidth',2);
  axis equal;
  axis([min(posx) - 0.1 max(posx) + 0.1 -0.71 0.71])

  figure(2); clf;
  subplot(3,1,1); hold on
  plot(xtar(1,istart:iend),lineartop)
  plot(xtar(1,s(j)),lineartop(s(j)),'k.','markersize',20)
  subplot(3,1,2); hold on
  plot(xtar(1,istart:iend),quadtop)
  plot(xtar(1,s(j)),quadtop(s(j)),'k.','markersize',20)
  subplot(3,1,3); hold on
  plot(xtar(1,istart:iend),ratiotop)
  plot(xtar(1,s(j)),ratiotop(s(j)),'k.','markersize',20)
  figure(4);
  plot(+0.7-ytar(indtop{s(j)},s(j)),uxtar(indtop{s(j)},s(j)))
  xlabel('Distance from top wall')
  ylabel('Velocity')
  figure(6); clf; hold on
  plot(dten);
  plot(spos{s(j)},dten(spos{s(j)}),'k.','markersize',20);

%  h = cline(posx,posy,ten);
%  set(h,'linewidth',2)
%  plot(wallsx,wallsy,'k','linewidth',2);
%  plot(xtar(indtop{s(j)},s(j)),ytar(indtop{s(j)},s(j)),'k--','linewidth',2);
%  axis equal;
%  axis([min(posx) - 0.1 max(posx) + 0.1 -0.71 0.71])
%  colorbar
  pause
end

end


% find indicies where quadratic term changes sign and when its
% derivative changes sign (max and mins)
s1 = find(quadbot(2:end).*quadbot(1:end-1) < 0);
dquadbot = smooth(diff(quadbot))';
% points where derivative changes sign (ie. max or min)
s2 = find(dquadbot(2:end).*dquadbot(1:end-1) < 0);

s = [s1 s2];
for j = 1:numel(s)
  s(j)
  figure(1); clf; hold on;
  %plot(xtar(:,istart:iend),ytar(:,istart:iend),'b.');
  plot(wallsx,wallsy,'k','linewidth',2);
  quiver(xtar(:,K),ytar(:,K),uxtar(:,K),uytar(:,K),'b')
%  h = fill(posx,posy,'r');
%  set(h,'linestyle','none')
  h = cline(posx,posy,conc);
  set(h,'linewidth',2)
  colorbar
  axis equal;
  axis([min(posx) - 0.1 max(posx) + 0.1 -0.71 0.71])

  plot(xtar(indbot{s(j)},s(j)),ytar(indbot{s(j)},s(j)),'k--','linewidth',2);
  figure(3); clf
  subplot(3,1,1); hold on
  plot(xtar(1,istart:iend),linearbot)
  plot(xtar(1,s(j)),linearbot(s(j)),'k.','markersize',20)
  subplot(3,1,2); hold on
  plot(xtar(1,istart:iend),quadbot)
  plot(xtar(1,s(j)),quadbot(s(j)),'k.','markersize',20)
  subplot(3,1,3); hold on
  plot(xtar(1,istart:iend),ratiobot)
  plot(xtar(1,s(j)),ratiobot(s(j)),'k.','markersize',20)
  figure(4);
  plot(+0.7+ytar(indbot{s(j)},s(j)),uxtar(indbot{s(j)},s(j)))
  xlabel('Distance from bottom wall')
  ylabel('Velocity')

  figure(6); clf; hold on
  plot(dten);
  plot(sneg{s(j)},dten(sneg{s(j)}),'k.','markersize',20);
%  h = cline(posx,posy,ten);
%  set(h,'linewidth',2)
%  plot(wallsx,wallsy,'k','linewidth',2);
%  plot(xtar(indbot{s(j)},s(j)),ytar(indbot{s(j)},s(j)),'k--','linewidth',2);
%  axis equal;
%  axis([min(posx) - 0.1 max(posx) + 0.1 -0.71 0.71])
%  colorbar

  pause
end

