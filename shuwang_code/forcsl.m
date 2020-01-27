function fsl = forcsl(m,theta,un)
% compute fourier derivative of opening tangent angle using "trick" by
% subtracting a function that goes from 0 to 2*pi over the range of
% alpha values
stheta = fd1x(theta,m);

temp(1,1:m) = stheta(1,1:m).*un(1,1:m);
temp(1,m+1) = temp(1,1);
tempi = fin(temp,m);

fsl = tempi(1,m+1);
end
% This routine is equation (61) in the Sohn et al 2010 JCP paper
      
