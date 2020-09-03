function [theta2,x0,y0] = initiallreconstruction(m,theta,sl)
% This routine finds a new opening angle and single tracker point so
% that it has the same area as the initial vesicle, but theta is
% band-limited to only 9 frequencies

% m - number of points on vesicle
% theta - opening angle
% sl - length of vesicle

% theta2 - new band-limited opening angle
% x0 - x-coordinate of tracker point
% y0 - y-coordinate of tracker point

x0=1;
y0=0;
 
[x,y] = recon(m,x0,y0,sl,theta);
% construct the x and y coordinates of the shape with length sl and
% length theta. (x0,y0) is the first point on the curve and m is the
% number of points
xRef = x-mean(x); yRef = y-mean(y);

a = theta(1);

areasum = sum(sin(theta(1,1:m)).*x(1,1:m) - ...
              cos(theta(1,1:m)).*y(1,1:m))/2*sl/m;
% area of the shape

thetas(1:m) = theta(1:m) - a - (0:m-1)*2*pi/m;
% subtract off the initial value and the linear function so that it is
% periodic
ktheta = fft(thetas,m);
% compute fourier series of thetas
ktheta(4:m-4) = 0;
% remove everything with frequency greater than 4
dtheta = real(ifft(ktheta,m));
% compute filtered theta in physical space
theta2 = dtheta + a + (0:m-1)*2*pi/m;
% add back on the constant term and linear function
[x,y] = recon(m,x0,y0,sl,theta2);
% compute the new geometry with the filtered theta

area = sum(sin(theta2(1:m)).*x(1:m)-...
           cos(theta2(1:m)).*y(1:m))/2*sl/m;
% compute the area of the modified vesicle shape)

k=1;
while abs(area-areasum)/areasum>1e-10
  coefficient = area - areasum;
  dtheta = dtheta*(1 + coefficient/30);
  theta2 = dtheta + a + (0:m-1)*2*pi/m;
  % scale the angle theta by an amount proportional to the difference
  % between the desired and true area
 
  [x,y] = recon(m,x0,y0,sl,theta2);
  x = x - mean(x(1,1:m));
  y = y - mean(y(1,1:m));
  % compute new x and y coordinates
%  clf; hold on
%  plot(xRef,yRef);
%  plot(x,y,'r--')
%  axis equal
%  pause

  area = sum(sin(theta2(1,1:m)).*x(1,1:m)-...
             cos(theta2(1,1:m)).*y(1,1:m))/2*sl/m;
  % find new area

  k = k+1;
  if k > 100
    break
  end
  % stop if it takes too many iterations

end

x0=x(1);
y0=y(1);

end

