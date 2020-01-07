function [x,y,theta,rcon,sl] = initiall(...
    shortax,ngrid,concentration,symmetry,smallper)

tol = 1e-12;
N = (0:1:ngrid);
x0 = cos(2*pi*N/ngrid);
y0 = shortax*sin(2*pi*N/ngrid);

dx = fd1(x0,ngrid);
dy = fd1(y0,ngrid);

slo = sqrt(dx.^2+dy.^2);
slo(ngrid+1) = slo(1);
alpha = arc(slo,ngrid);

x = cos(pi*2*alpha);
y = shortax*sin(pi*2*alpha);
sxn = fd1(x,ngrid);
syn = fd1(y,ngrid);
sl = sum(sqrt(sxn(1:ngrid).^2+syn(1:ngrid).^2));
sl = sl/ngrid;
if abs(sxn(1)) > tol
  t0 = atan(abs(syn(1)/sxn(1)));
else
  t0 = pi/2;
end
theta = thetsolve(sxn,syn,sl,t0,ngrid);

if symmetry == -1
  rcon=concentration+smallper*random('unif',-5,5,1,ngrid+1);
elseif symmetry==0
  rcon = concentration + 3*smallper*cos(alpha*pi*2) + ...
    0.5*smallper*cos(3*alpha*pi*2) + ...
    0.5*smallper*cos(4*alpha*pi*2);
elseif symmetry==1
  rcon = concentration+5*smallper*cos(2*alpha);        
elseif symmetry==2
  rcon = concentration+5*smallper*sin(alpha*pi*2); 
end

end
