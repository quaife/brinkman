addpath ../../src
set(0,'DefaultAxesFontSize',22)

irate = 10; % controls the speed of the visualization

if 0
%  file = 'extensional2VesJData.bin';
  file = '~/projects/brinkman/vesicle_code/results/extensional2Ves/adR4em1adS7em1Chi1em1_ra075/extensional2VesData.bin';
  ax = [-4 4 -5 5];
  options.confined = false;
end
if 0
  file = 'extensionalManyVesData.bin';
  ax = [-12 12 -2 2];
  options.confined = false;
end
if 0
  file = 'relaxation1VesData.bin';
  ax = ([-3 3 -3 3]);
  options.confined = false;
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
  file = 'shear1VesAData.bin.run1';
  ax = [-8 8 -5 5];
  options.confined = false;
end
if 1
  file = 'shear2VesHData.bin.run1';
%  file = '~/projects/brinkman/vesicle_code/results/shear2Ves/adR4em1adS7em1Chi5em1_ra050/shear2VesData.bin';
  ax = [-10 10 -3 3];
  options.confined = false;
end
if 0
  file = 'choke1VesData.bin';
  ax = [-10.5 10.5 -3.5 3.5];
  options.confined = true;
end

[posx,posy,ten,wallx,wally,ea,el,time,n,nv] = loadFile(file);
%[posx,posy,ten,~,~,wallx,wally,ea,el,time,n,nv] = loadFileOld(file);
% load positions, tension, stresses, errors, time, number of points, and
% number of vesicles
istart = 1;
ntime = numel(time);

%min_ten = floor(min(min(min(ten))));
%max_ten = floor(max(max(max(ten))));
%if (max_ss - min_ss < 1e-5)
%  min_ss = min_ss - 
%  max_ss = max_ss + 1;
%end

%istart = 19000;
iend = numel(time);

figure(1); clf
for k = istart:irate:iend
  xx = interpft(posx(:,:,k),96); yy = interpft(posy(:,:,k),96);  
  vec1 = [xx(:,:);xx(1,:)];
  vec2 = [yy(:,:);yy(1,:)];
  if 1
    clf
    plot(vec1,vec2,'r','linewidth',3)
    hold on;
    plot(vec1(1,:),vec2(1,:),'b.','markersize',20)
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
%  caxis([min_ten max_ten])
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
%  if options.savefig
%    filename = ['./frames/image', sprintf('%04d',count),'.png'];
%    count = count+1;
%    figure(1);
%    print(gcf,'-dpng','-r300',filename);
%  end
  pause(0.01)
%  pause
end



