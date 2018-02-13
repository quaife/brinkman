set(0,'DefaultAxesFontSize',22)
options.tracers = false; % plot saved tracers
options.pressure = false; % plot saved pressures
options.stress = false; % plot the saved stresses
options.quiver = false;
options.jacobian = false;
options.dist = false;
options.angle = false;

irate = 1; % controls the speed of the visualization


if 0
  file = 'shear1VesData.bin';
  ax = 5*[-1 1 -1 1];
  options.confined = false;
end
if 1
  file = 'shear2VesData.bin';
  ax = [-10 10 -5 5];
  options.confined = false;
end
if 0
  file = 'shear4VesData.bin';
  ax = [-10 10 -5 5];
  options.confined = false;
end
if 0
  file = 'taylor1VesData.bin';
  ax = 2*[0 pi 0 pi];
  options.confined = false;
end

[posx,posy,ten,wallx,wally,ea,el,time,n,nv] = loadFile(file);
% load positions, tension, errors, time, number of points, and number
% of vesicles
ntime = numel(time);

figure(1); clf
for k = 1:irate:ntime
  xx = interpft(posx(:,:,k),96); yy = interpft(posy(:,:,k),96);  
  vec1 = [xx(:,:);xx(1,:)];
  vec2 = [yy(:,:);yy(1,:)];
  plot(vec1,vec2,'r-','linewidth',2)

  axis equal
  axis(ax);
  titleStr = ['t = ' num2str(time(k),'%4.2e') ...
      ' eA = ' num2str(ea(k),'%4.2e') ...
      ' eL = ' num2str(el(k),'%4.2e')];
  title(titleStr)
  
  if options.confined
    hold on
    vec1 = [wallx;wallx(1,:)];
    vec2 = [wally;wally(1,:)];
    plot(vec1,vec2,'k','linewidth',3);
    hold off
  end
  set(gca,'xtick',[]);
  set(gca,'ytick',[]);
  set(gca,'xcolor','w');
  set(gca,'ycolor','w');

%  if options.savefig
%    filename = ['./frames/image', sprintf('%04d',count),'.png'];
%    count = count+1;
%    figure(1);
%    print(gcf,'-dpng','-r300',filename);
%  end
  pause(0.01)
end

if (options.angle || options.jacobian)
  IA = zeros(ntime,1);
  jacobian = zeros(n,ntime);
  % compute inclination angle on an upsampled grid
  N = 1024;
  modes = [(0:N/2-1)';0;(-N/2+1:-1)'];

  for k = 1:ntime
    x = interpft(posx(:,:,k),N); y = interpft(posy(:,:,k),N);
    Dx = real(ifft(1i*modes.*fft(x)));
    Dy = real(ifft(1i*modes.*fft(y)));
    jac = sqrt(Dx.^2 + Dy.^2);
    tx = Dx./jac; ty = Dy./jac;
    nx = ty; ny = -tx;
    rdotn = x.*nx + y.*ny;
    rho2 = x.^2 + y.^2;

    J11 = 0.25*sum(rdotn.*(rho2 - x.*x).*jac)*2*pi/N;
    J12 = 0.25*sum(rdotn.*(-x.*y).*jac)*2*pi/N;
    J21 = 0.25*sum(rdotn.*(-y.*x).*jac)*2*pi/N;
    J22 = 0.25*sum(rdotn.*(rho2 - y.*y).*jac)*2*pi/N;

    J = [J11 J12; J21 J22];
    [V,D] = eig(J);

    [~,ind] = min(abs(diag(D)));
    IA(k) = atan(V(2,ind)/V(1,ind));
    jacobian(:,k) = interpft(jac,n);
  end
  IA(1) = pi/2;
end



