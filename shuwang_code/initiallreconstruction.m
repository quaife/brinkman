function [theta2,sl,k,x0,y0] = initiallreconstruction(m,theta,sl)
    
x0=1;
y0=0;
 
[x,y]=recon(m,x0,y0,sl,theta);

a=theta(1);

areasum = sum(sin(theta(1,1:m)).*x(1,1:m) - ...
              cos(theta(1,1:m)).*y(1,1:m))/2*sl/m;

thetas(1,1:m) = theta(1,1:m)-a-[0:m-1]*2*pi/m;
ktheta = fft(thetas,m);
ktheta(4:m-4) = 0;
dtheta = real(ifft(ktheta,m));
theta2 = dtheta+a+[0:m-1]*2*pi/m;
[x,y] = recon(m,x0,y0,sl,theta2);

area = sum(sin(theta(1,1:m)).*x(1,1:m)-...
            cos(theta(1,1:m)).*y(1,1:m))/2*sl/m;

k=1;
while abs(area-areasum)/areasum>1e-10
  coefficient = area-areasum;
  dtheta = dtheta*(1+coefficient/3);
  theta2 = dtheta + a + [0:m-1]*2*pi/m;
 
  [x,y] = recon(m,x0,y0,sl,theta2);
  x = x-mean(x(1,1:m));
  y = y-mean(y(1,1:m));

  area = sum(sin(theta2(1,1:m)).*x(1,1:m)-...
            cos(theta2(1,1:m)).*y(1,1:m))/2*sl/m;

  k=k+1;
  if k>100
    break
  end

end

x0=x(1);
y0=y(1);

end

