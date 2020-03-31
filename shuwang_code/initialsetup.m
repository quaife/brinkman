function [x,y,theta,rcon,sl] = initialsetup( ...
  shortax,ngrid,concentra,oddeven)
% shortax - size of the short axis. Default value of long axis is 1 (I
% expect)
% ngrid - number of discretization points
% concentra - initial concentration of the lipid species?
% oddeven - sets the initial concetration profile?

% x - x coordinate
% y - y coordinate
% theta - opening angle
% rcon - lipid species concentration
% sl - length of vesicle

smallperturbation = 5e-2;
if concentra == 0
  smallperturbation = 0;
end
    
[x,y,theta,rcon,sl] = initiall(shortax,ngrid*4,concentra,...
    oddeven,smallperturbation);
% initialize the vesicle shape, opening angle, length, and concentration
% of lipid species

theta = theta(1:4:length(theta));
% upsampled by a factor of 4 when calling initiall. Now need to
% downsample back to ngrid
rcon = rcon(1:4:4*ngrid+1);
% downsample the concentration to the correct number of discretization
% points

figure(1); clf; hold on
plot(acurv(sl,theta(1:end-1),64),'b')

[theta,x0,y0] = initiallreconstruction(ngrid,theta,sl);
% Find a vesicle shape with a similar area (within 1e-10), but with a
% theta that is band-limited to only frequencies between -4 and 4.
size(theta)
plot(acurv(sl,theta(1:end),64),'r--')
pause

[x,y] = recon(ngrid,x0,y0,sl,theta);
% find x and y coordinates of the band-limited theta. Note that theta
% has already been downsampled, so x and y have ngrid discretization
% points

end
