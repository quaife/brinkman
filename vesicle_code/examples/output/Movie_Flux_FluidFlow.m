addpath ../../src

set(0,'DefaultAxesFontSize',22)
set(gcf,'Position',get(0,'Screensize'));
ax = [-2.3 2.3 -1.75 1.75];

file = 'Fig3pt6_remake_Q2019_SPbpt1.bin';
name = 'Fig3pt6_remake_Q2019_SPbpt1';
[posx,posy,ten,~,~,~,~,time,n,nv] = loadFile(file);

oc = curve;
op = poten(n);

istart = 1;
irate = 10; % controls the speed of the visualization
iend = numel(time);
ntime = numel(time);

figure(1); clf
for k = istart:irate:iend 
  xx = posx(:,:,k);
  yy = posy(:,:,k);
  tens = ten(:,:,k); 
  vec1 = [xx(:,:);xx(1,:)];
  vec2 = [yy(:,:);yy(1,:)];

  ves = capsules([xx;yy],tens,[],1,1);
  f = ves.tracJump([xx;yy],tens);
  Pf = ves.normalProjection(f);
  G = op.stokesSLmatrix(ves);
  for j = 1:nv
    Gf(:,j) = G(:,:,j)*f(:,j); 
  end
  
  v = Gf + [yy;zeros(length(yy),nv)];

  [vx,vy] = oc.getXY(v);
  [~,xt,~] = oc.diffProp([xx;yy]);
  [tanx,tany] = oc.getXY(xt);
  nx = -tany;
  ny = tanx;

  vdotn = vx.*nx + vy.*ny;

  xarr = linspace(-2.3,2.3,100);
  yarr = linspace(-1.75,1.75,100);
  [Xarr,Yarr] = meshgrid(xarr,yarr);
  Xarr = Xarr(:);
  Yarr = Yarr(:);
  Tpoints = capsules([Xarr;Yarr],[],[],[],[]);

  kernel = @op.exactStokesSL;
  kernelDirect = @op.exactStokesSL;
  Galpert = op.stokesSLmatrix(ves);

  SLP = @(X) op.exactStokesSLdiag(ves,Galpert,X);
  [~,NearV2T] = ves.getZone(Tpoints,2);
  Fslp = op.nearSingInt(ves,f,SLP,...
      NearV2T,kernel,kernelDirect,Tpoints,false,false,false);
  
  
  if 1
    clf
    %plot(vec1,vec2,'r','linewidth',3)
    scale = 1;
    %shear
    %quiver(Xarr,Yarr,Fslp(1:end/2)+Yarr,Fslp(end/2+1:end),scale)
    %relaxed
    quiver(Xarr,Yarr,Fslp(1:end/2),Fslp(end/2+1:end),scale)

    hold on;
    %quiver(xx,yy,vx,vy)
    for j = 1:nv
        h = cline(xx(:,j),yy(:,j),vdotn(:,j)); 
        set(h,'linewidth',5)

    end
    plot(vec1(1,:),vec2(1,:),'b.','markersize',20)
    cb = colorbar;
    set(cb, 'ylim', [-2 20])
    %colorbar
   
    view(2); hold off;
    axis equal
    axis(ax)
    titleStr = ['t = ' num2str(time(k),'%4.2e')];
    title(titleStr)
  end
  figname=sprintf('%s_%s_%04.f%s',name,'frame',k);
  %saveas(gcf, figname)
  set(gcf,'Units','inches');
  screenposition = get(gcf,'Position');
  set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
  print(figname,'-dpdf','-painters')
  %pause(0.01)
end




