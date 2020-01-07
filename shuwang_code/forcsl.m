function fsl = forcsl(m,theta,un)

stheta = fd1x(theta,m);

temp(1,1:m) = stheta(1,1:m).*un(1,1:m);
temp(1,m+1) = temp(1,1);
tempi = fin(temp,m);

fsl = tempi(1,m+1);
end
      
