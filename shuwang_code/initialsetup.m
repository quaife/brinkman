function [x,y,theta,rcon,sl] = initialsetup( ...
  shortax,ngrid,concentra,oddeven)

smallperturbation = 5e-2;
if concentra == 0
  smallperturbation = 0;
end
    
[x,y,theta,rcon,sl] = initiall(shortax,ngrid*4,concentra,...
    oddeven,smallperturbation);
 
% Here we set the position of the first point, which is used in the
% integral
ulam = 1;
theta = theta(1:4:length(theta));

[theta,sl,kstop,x0,y0] = initiallreconstruction(ngrid,theta,sl);

[x,y] = recon(ngrid,x0,y0,sl,theta);
rcon = rcon(1:4:4*ngrid+1);







end
