addpath ..
% file1 = 'Chip25_shax5p7_scL0p311_Conc0p5_Beta0_n1024_nbd1024_dt1en4_bmax1_bmin0p1_eps0p04_a100_longchoke_rand.bin';
% [posx1,posy1,conc1,ea1,el1,time1,xvel1,yvel1,ten1] = loadFile(file1);


irate = 1; 
istart = 1;
iend = numel(time1);

N = length(posx1(:,:,1));
oc = curve(N);

if 0
    Nbd = 1024;
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


count = 1;
Tangent1 = [];
Tangent2 = [];
Tangent3 = [];
Tangent4 = [];

for k = istart:irate:iend
  %==================== Phase & Reference Point Velocities
  [~, tang, ~] = oc.diffProp([posx1(:,:,k);posy1(:,:,k)]);
  %dummy = tang(1:end/2).*xvel1(:,:,k)+tang(end/2+1:end).*(yvel1(:,:,k);
  %Tangent(count) = dummy(1);
  Tangent1(count) = tang(1)*(xvel1(1,:,k)-mean(squeeze(xvel1(:,:,k)))) + tang(N+1)*yvel1(1,:,k);
  Tangent2(count) = tang(N/4)*(xvel1(N/4,:,k)-mean(squeeze(xvel1(:,:,k)))) + tang(N/4+N+1)*yvel1(N/4,:,k);
  Tangent3(count) = tang(N/2)*(xvel1(N/2,:,k)-mean(squeeze(xvel1(:,:,k)))) + tang(N/2+N+1)*yvel1(N/2,:,k);
  Tangent4(count) = tang(3*N/4)*(xvel1(3*N/4,:,k)-mean(squeeze(xvel1(:,:,k)))) + tang(3*N/4+N+1)*yvel1(3*N/4,:,k);
  %Tangent2(count) = tang(N/2+1).*(xvel(N/2+1,:,k)-mean(squeeze(xvel1(:,:,ind)))) + tang(N+1+N/2).*yvel(N/2+1,:,k); 
%   [~,indmax] = max(conc1(:,:,k));
%   indleft = indmax - round(ubar*N/2);
  %indright = indmax + round(ubar*N/2);
%   [~, tangphase, ~] = oc.diffProp([posx1(indright,:,k);posy1(indright,:,k)]);
%   TangentPhase(count) = tangphase(1).*xvel(indright,:,k) + tangphase(2).*yvel(indright,:,k);
%   if(indleft < 0)
%     maxleft = 1 - round(ubar*N/2);
%     diff = indleft - maxleft;
%     [~, tangphase, ~] = oc.diffProp([posx1(1:diff,:,k);posy1(1:diff,:,k)]);
%     TangentPhase(count) = mean(tangphase(1:end/2).*xvel(indleft:indright,:,k)+tangphase(end/2+1:end).*yvel(indleft:indrigh,:,k));
%   else
%     [~, tangphase, ~] = oc.diffProp([posx1(indleft:indright,:,k);posy1(indleft:indright,:,k)]);
%     TangentPhase(count) = mean(tangphase(1:end/2).*xvel(indleft:indright,:,k)+tangphase(end/2+1:end).*yvel(indleft:indrigh,:,k));
%   end

  count = count + 1;
end
figure(1)
plot(mean(squeeze(posx1(:,:,istart:irate:iend))),Tangent1, 'r', 'linewidth', 3)
hold on
plot(mean(squeeze(posx1(:,:,istart:irate:iend))),Tangent2, 'b', 'linewidth', 3)
plot(mean(squeeze(posx1(:,:,istart:irate:iend))),Tangent3, 'g', 'linewidth', 3)
plot(mean(squeeze(posx1(:,:,istart:irate:iend))),Tangent4, 'm', 'linewidth', 3)
xlabel('Mean x position','fontsize',16)
ylabel('Tank treading velocity','fontsize',16)
title('$\alpha = 0.6$','Interpreter','latex','fontsize',24)
xline(11,'k--')
legend('$1$','$\frac{N}{4}$','$\frac{N}{2}$','$\frac{3N}{4}$','Interpreter','latex','fontsize',16,'location','northwest')
legend boxoff
xlim([0 14])

ind = find(mean(squeeze(posx1))>=11);
ind = ind(1);
bmin = 0.1; bmax = 1;

figure(2)
colormap(turbo)
rbn = (bmin - bmax)/2*tanh(3*(conc1(:,:,ind) - 0.5)) + (bmax + bmin)/2;
%h = cline(posx1(:,:,ind),posy1(:,:,ind),rbn);
plot(posx1(:,:,ind),posy1(:,:,ind),'color', [0.4796 0.0158 0.0106],'linewidth',3)
hold on
plot3(posx1(1,:,ind),posy1(1,:,ind),10,'r.','markersize',25)
plot3(posx1(N/4,:,ind),posy1(N/4,:,ind),10,'b.','markersize',25)
plot3(posx1(N/2,:,ind),posy1(N/2,:,ind),10,'g.','markersize',25)
plot3(posx1(3*N/4,:,ind),posy1(3*N/4,:,ind),10,'m.','markersize',25)

axis equal
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'xcolor','white')
set(gca,'ycolor','white')
set(gca,'Visible','off')
view(2)
