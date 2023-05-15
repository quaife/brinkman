usePlot = false;
addpath ..
%addpath D:\reserach\longchoke\bmaxp1\chip25\RAp6

% file = 'Chi2p5_shax4p37_scL0p397_Conc0p75_Beta0_n1024_nbd1024_dt5en6_bmax1_bmin0p1_eps0p04_a100_contracting_left.bin';
% [posx1,posy1,~,~,~,time1,~,~,~] = loadFile(file);

Nbd = 2048;
geomCenter = [0;0];
wallGeometry = 'contracting';
oc = curve(Nbd);
[~,Xwalls] = oc.initConfig(Nbd,false,...
             'scale', 1, ...
             'center', geomCenter, 'geometry', wallGeometry);
[xwalls,ywalls] = oc.getXY(Xwalls);

irate = 10; 
istart = 1;
iend = numel(time1);
ntime = iend;
ntime = numel(time1);
N = length(posx1(:,:,1));
% time = time(istart:iend);

% points that are more than this distance are excluded from the boundary
% layer points
tol = 0.3;

TopBLWidth = [];
BotBLWidth = [];
for k=istart:irate:iend
  xx = posx1(:,:,k);
  yy = posy1(:,:,k);
  if usePlot
    figure(1); clf; hold on;
    plot(xx,yy,'r');
    plot(xwalls,ywalls,'k');
    axis equal;
    xmin = min(xx) - 1;
    xmax = max(xx) + 1;
    %axis([xmin xmax -0.8 0.8])
  end
  dist2Top = zeros(size(xx));
  for j=1:length(xx)
      % START OF TOP BOUNDARY LAYER
      disty2Top = yy(j) - ywalls(1:Nbd/2);
      distx2Top = xx(j) - xwalls(1:Nbd/2);
      dist2Top(j) = min(sqrt(distx2Top.^2+disty2Top.^2));
  end
  % include a bit of overlap in case the first discretization point is a
  % local minimum
  dist2Top = [dist2Top; dist2Top(1:5)];
    
  % first derivative of distance from top
  Ddist2Top =  -diff(dist2Top);
  % second derivative of distance from top
  D2dist2Top = -diff(Ddist2Top);

  % find all indicies where the derivate changes sign, indicating that
  % it correspond to a local max or min.
  s = find(Ddist2Top(1:end-1).*Ddist2Top(2:end) < 0);
  % only care about local maxs of the distance function, meaning that it's
  % as close as possible.

  % remove points that are indexed number than N since this means its
  % one of the first discretization points that we have already counted.
  s = s(s<N);

  s2 = find(D2dist2Top(s) < 0);
  s = s(s2);

  % only care about points that are sufficiently close to the boundary
  s2 = find(dist2Top(s) < tol);
  s = s(s2);

  if numel(s) == 1
    % find the two points that are closest to the horizontal line y = tol
%    [~,ind] = sort(abs(posy(:,:,k) - tol));
%    ind = ind(1:2);
    ind = find(dist2Top < tol);
    ind = ind(ind<=N);
    [~,indRight] = max(posx1(ind,:,k));
    [~,indLeft] = min(posx1(ind,:,k));
    indRight_Top = ind(indRight);
    indLeft_Top  = ind(indLeft);
  else
    [~,indLeft]  = min(posx1(s,:,k));
    indLeft_Top = s(indLeft);
    [~,indRight] = max(posx1(s,:,k));
    indRight_Top = s(indRight);
  end
  % END OF TOP BOUNDARY LAYER

  % START OF BOTTOM BOUNDARY LAYER
  dist2Bot = zeros(size(xx));
  for j=1:length(xx)
      % START OF TOP BOUNDARY LAYER
      disty2Bot = yy(j) - ywalls(end/2+1:end);
      distx2Bot = xx(j) - xwalls(end/2+1:end);
      dist2Bot(j) = min(sqrt(distx2Bot.^2+disty2Bot.^2));
  end
  % include a bit of overlap in case the first discretization point is a
  % local minimum
  dist2Bot = [dist2Bot; dist2Bot(1:5)];

  % first derivative of distance from top
  Ddist2Bot =  -diff(dist2Bot);
  % second derivative of distance from top
  D2dist2Bot = -diff(Ddist2Bot);

  % find all indicies where the derivate changes sign, indicating that
  % it correspond to a local max or min.
  s = find(Ddist2Bot(1:end-1).*Ddist2Bot(2:end) < 0);
  % only care about local maxs of the distance function, meaning that it's
  % as close as possible.

  % remove points that are indexed number than N since this means its
  % one of the first discretization points that we have already counted.
  s = s(s<N);

  s2 = find(D2dist2Bot(s) > 0);
  s = s(s2);

  % only care about points that are sufficiently close to the boundary
  s2 = find(dist2Bot(s) < tol);
  s = s(s2);
  

  if numel(s) == 1
    % find the two points that are closest to the horizontal line y = -tol
%    [~,ind] = sort(abs(posy(:,:,k) + tol));
%    ind = ind(1:2);
    ind = find(dist2Bot < tol);
    ind = ind(ind<=N);
    [~,indRight] = max(posx1(ind,:,k));
    [~,indLeft] = min(posx1(ind,:,k));
    indRight_Bot = ind(indRight);
    indLeft_Bot  = ind(indLeft);
  else
    [~,indLeft]  = min(posx1(s,:,k));
    indLeft_Bot = s(indLeft);
    [~,indRight] = max(posx1(s,:,k));
    indRight_Bot = s(indRight);
  end
  % END OF BOTTOM BOUNDARY LAYER

  if usePlot
    figure(1); hold on;
    % plot the left and right end points for the boundary layers
    plot(xx(indLeft_Top),yy(indLeft_Top),'kx','markersize',20)
    plot(xx(indRight_Top),yy(indRight_Top),'kx','markersize',20)
    plot(xx(indLeft_Bot),yy(indLeft_Bot),'kx','markersize',20)
    plot(xx(indRight_Bot),yy(indRight_Bot),'kx','markersize',20)
  end


  % find index of points that are along the top boundary layer. Care is
  % taken in case the first discretization point is between the two
  % boundaries
  if (indRight_Top > indLeft_Top)
    ind_Top = [(indRight_Top:N)';(1:indLeft_Top)'];
    %ind_Top = (indLeft_Top:indRight_Top)';
  else
    ind_Top = (indRight_Top:indLeft_Top)';
  end
  ind_Top = ind_Top(dist2Top(ind_Top) < tol);

  % find index of points that are along the bottom boundary layer. Care
  % is taken in case the first discretization point is between the two
  % boundaries
  if (indRight_Bot > indLeft_Bot)
    ind_Bot = (indLeft_Bot:indRight_Bot)';
  else
    ind_Bot = [(indLeft_Bot:N)';(1:indRight_Bot)'] ;
  end
  ind_Bot = ind_Bot(dist2Bot(ind_Bot) < tol);

  if usePlot
    figure(1); hold on
    plot(xx(ind_Top),yy(ind_Top),'b.','markersize',20)
    plot(xx(ind_Bot),yy(ind_Bot),'g.','markersize',20)
    pause
  end


  TopBLWidth = [TopBLWidth;mean(dist2Top(ind_Top))];
  BotBLWidth = [BotBLWidth;mean(dist2Bot(ind_Bot))];
  if usePlot
    figure(3);
    hold on
    plot(abs(TopBLWidth),'b')
    plot(abs(BotBLWidth),'k')
  end

end



figure(4);
hold on
plot(mean(squeeze(posx1(:,:,istart:irate:iend))),smoothdata(abs(TopBLWidth),'gaussian',200),'g','linewidth',3)
plot(mean(squeeze(posx1(:,:,istart:irate:iend))),smoothdata(abs(BotBLWidth),'gaussian',200),'g--','linewidth',3)
% plot(mean(squeeze(posx1(:,:,istart:irate:iend))),abs(TopBLWidth),'r')
% plot(mean(squeeze(posx1(:,:,istart:irate:iend))),abs(BotBLWidth),'r--')
%xlim([0 18])
% legend('Top','bottom')