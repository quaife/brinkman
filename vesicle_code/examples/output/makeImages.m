set(0,'DefaultAxesFontSize',22)

irate = 5; % controls the speed of the visualization

if 1
  file = 'shear2VesData.bin';
  ax = 5*[-1 1 -1 1];
  options.confined = false;
end
if 0
  file = 'shear4VesData.bin';
  ax = [-10 10 -5 5];
  options.confined = false;
end

[posx,posy,ten,shearStress,normalStress,...
    wallx,wally,ea,el,time,n,nv] = loadFile(file);
% load positions, tension, stresses, errors, time, number of points, and
% number of vesicles
ntime = numel(time);

min_ten = floor(min(min(min(ten))));
max_ten = floor(max(max(max(ten))));
min_ss = floor(min(min(min(shearStress))));
max_ss = floor(max(max(max(shearStress))));
min_ns = floor(min(min(min(normalStress))));
max_ns = floor(max(max(max(normalStress))));


figure(1); clf
for k = 1:irate:ntime
  xx = interpft(posx(:,:,k),96); yy = interpft(posy(:,:,k),96);  
  tt = interpft(ten(:,:,k),96);
  ss = interpft(shearStress(:,:,k),96);
  ns = interpft(normalStress(:,:,k),96);

  vec1 = [xx(:,:);xx(1,:)];
  vec2 = [yy(:,:);yy(1,:)];

  clf;
  subplot(1,3,1);hold on
  vec3 = [tt(:,:);tt(1,:)];
  for j = 1:nv
    h = cline(vec1(:,j),vec2(:,j),vec3(:,j));
    set(h,'LineWidth',4);
  end
  axis equal
  axis(ax);
  set(gca,'xtick',[]);
  set(gca,'ytick',[]);
  set(gca,'xcolor','w');
  set(gca,'ycolor','w');
  title('Tension')
  colorbar
  caxis([min_ten max_ten])

  subplot(1,3,2);hold on
  vec3 = [ss(:,:);ss(1,:)];
  for j = 1:nv
    h = cline(vec1(:,j),vec2(:,j),vec3(:,j));
    set(h,'LineWidth',4);
  end
  axis equal
  axis(ax);
  set(gca,'xtick',[]);
  set(gca,'ytick',[]);
  set(gca,'xcolor','w');
  set(gca,'ycolor','w');
  title('Shear Stress')
  colorbar
  caxis([min_ss max_ss])

  subplot(1,3,3);hold on
  vec3 = [ns(:,:);ns(1,:)];
  for j = 1:nv
    h = cline(vec1(:,j),vec2(:,j),vec3(:,j));
    set(h,'LineWidth',4);
  end
  axis equal
  axis(ax);
  set(gca,'xtick',[]);
  set(gca,'ytick',[]);
  set(gca,'xcolor','w');
  set(gca,'ycolor','w');
  title('Normal Stress')
  colorbar
  caxis([min_ns max_ns])
  
  titleStr = ['t = ' num2str(time(k),'%4.2e') ...
      ' eA = ' num2str(ea(k),'%4.2e') ...
      ' eL = ' num2str(el(k),'%4.2e')];
  suptitle(titleStr)
%  if options.savefig
%    filename = ['./frames/image', sprintf('%04d',count),'.png'];
%    count = count+1;
%    figure(1);
%    print(gcf,'-dpng','-r300',filename);
%  end
  pause(0.01)
end




