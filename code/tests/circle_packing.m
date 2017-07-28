function [xc,yc,rc] = circle_packing(n_bodies)
% input the number of desired bodies
% return the x and y coordinates of the circles and their radii

xc = []; yc = []; rc = [];

rmax = 0.7;
rcircle = 0.04;

rc = (rmax - rcircle)*rand;
theta = 2*pi*rand;
xc = rc*cos(theta);
yc = rc*sin(theta);
% first point

N = 128;
omega = (0:N-1)'*2*pi/N;
clf; hold on
fill(xc(1)+rcircle*cos(omega),yc(1)+rcircle*sin(omega),'k')
axis equal;
axis([-1 1 -1 1])
drawnow()

while (numel(xc) < n_bodies)
  rp = (rmax - rcircle)*rand;
  theta = 2*pi*rand;
  xp = rp*cos(theta);
  yp = rp*sin(theta);

  iinner = check_inner(xp,yp,rcircle,xc,yc,rcircle);
  % check if the new circle intersects any of the other circles
  if ~iinner
    xc = [xc xp];
    yc = [yc yp];
    rc = [rc rp];
    % if no conflicts, add the new centers and radii to the running list
    fill(xc(end)+rcircle*cos(omega),yc(end)+rcircle*sin(omega),'k')
    drawnow()
    disp(numel(xc))
  end
end



%fid = fopen('thlen.dat','w');
%fprintf(fid,'%d\n',N);
%fprintf(fid,'%d\n',n_bodies);
%for k = 1:n_bodies
%  fprintf(fid,'%20.16e\n',theta+pi/2+pi/N);
%  fprintf(fid,'%20.16e\n',[2*pi*radii(k),xc(k),yc(k)]);
%end
%fclose(fid);


end


%%%%%%%%%%%%%%%%%%%%
function iinner = check_inner(xp,yp,rp,xc,yc,rc)
%[xp yp rp]
%[xc yc rc]
%pause

buffer = 1.01;
% buffer == 1 => circles can be infintesimally close
dist = sqrt((xp - xc).^2 + (yp - yc).^2);
iinner = any(dist < buffer*rc + buffer*rp);

end

