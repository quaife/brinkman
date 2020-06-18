classdef capsules < handle

properties
N; % number of discretization points
X; % positions on the vesicle
cur; % curvature
theta; % opening tangent angle
L; % length
rcon; % concentration field
x0; % single tracker point
y0; % single tracker point
bendsti; % max bending stiffness
bendratio; % bending ratio
viscIn;
viscOut;

end % properties

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = capsules(X,rcon,params)
% constructor

oc = curve;

o.N = size(X,1)/2;
o.X = X;
[o.L,o.theta,o.cur] = oc.computeOpeningAngle(o.N,X);
o.rcon = rcon;
o.x0 = X(1);
o.y0 = X(o.N + 1);
o.bendsti = params.bendsti;
o.bendratio = params.bendratio;
o.viscIn = params.viscosityInside;
o.viscOut = params.viscosityOutside;

end % vesicle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function smoothGeom(ves)
% rewrite of the routine intiallreconstruction. Goal is to find a
% vesicle shape with same area and length, but with a bandlimited
% opening angle

oc = curve;

N = ves.N;
X = ves.X;
theta = ves.theta;
L = ves.L;
x0 = ves.x0;
y0 = ves.y0;
X = oc.recon(N,x0,y0,L,theta);

areaRef = sum(sin(theta).*X(1:end/2) - cos(theta).*X(end/2+1:end))*...
      0.5*L/N;

theta0 = theta(1);
theta_periodicPart = theta - theta0 - (0:N-1)'*2*pi/N;
%  Compute the periodic part of the opening angle
thetah = fft(theta_periodicPart);
thetah(6:N-4) = 0;
% QUESTION: IN SHUWANG'S ORIGINAL CODE, THIS WAS NOT SET UP
% SYMETRICALLY AROUND THE ZERO MODE. CHECKED AND THE IMAG PART OF thetah
% IS NON-ZERO UNLESS THE COEFFICIENTS ARE ADJUSTED AS IN ABOVE
% remove all modes with Frequency above 4. Note that Shuwang's original
% code had an error since he did not truncate the positive and negative
% coefficients to the same level. Therefore, the imaginary part after
% taking an ifft was not 0.
theta_periodicPart = real(ifft(thetah));
theta = theta_periodicPart + theta0 + (0:N-1)'*2*pi/N;

X = oc.recon(N,x0,y0,L,theta);
% new shape with filtered opening angle

area = sum(sin(theta).*X(1:end/2) - cos(theta).*X(end/2+1:end))*...
      0.5*L/N;

iter = 1;
while abs(area - areaRef)/areaRef > 1e-10
  theta_periodicPart = theta_periodicPart * ...
        (1 + (area - areaRef)/3);
  theta = theta_periodicPart + theta0 + (0:N-1)'*2*pi/N;

  X = oc.recon(N,x0,y0,L,theta);
  % new shape with scaled periodic part of theta

  area = sum(sin(theta).*X(1:end/2) - cos(theta).*X(end/2+1:end))*...
        0.5*L/N;
  % new area

  iter = iter + 1;
  if iter > 100
    break
  end
  % break if it is taking too long
end

ves.X = X;
[ves.L,ves.theta,ves.cur] = oc.computeOpeningAngle(N,X);
% replace the geometry with the new shape


end % smoothGeom

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Eu,Esigma] = variationsNonStiff(ves)
% [Eu,Esigma] = variationsNonStiff() computes the non-stiff parts of the
% variations defined in equations (13) and (14) of the Sohn et al JCP
% paper. This corresonds to the last two terms in equation (13) and the
% first two terms in equaiton (14). Note the spotaneous curvature is 0,
% so the last term in equation (13) vanishes and the first two terms in
% equation (14) simplify

oc = curve;
N = ves.N;
IK = oc.modes(N);
rcon = ves.rcon;
cur = ves.cur;

b0 = ves.bendsti;
b1 = ves.bendsti * ves.bendratio;

rbn = b0 * (ones(N,1) - rcon) + b1*rcon;
Drbn = oc.diffFT(rbn,IK)/(2/pi);

Drbn_cur = oc.diffFT(rbn.*cur,IK)/(ves.L/2/pi);
DDrbn_cur = oc.diffFT(Drbn_cur,IK)/(ves.L/2/pi);

Esigma = -DDrbn_cur - 0.5*rbn.*cur.^3;
%Esigma is equation (14) with spotaneous curvature set to zero.

Eu = -Drbn.*cur.^2;
%Eu is the second term in equation (13) (differs by a negative
%sign - possibly from the negative sign in eq(23) which has a negative on 
%the variation for u). The last term drops since spontaneous curvature is 
%0. The first term is not in this routine since we are only calculating 
%variations due to changes in the vesicle shape and not the lipid species.
%** SHUWANG QUESTION: THIS IS A PLUS SIGN IN THE PAPER (EQUATION (13)),
% BUT IS A MINUS SIGN IN SHUWANG'S CODE **OLD COMMENT???

%clf
%plot(Drbn)
%max(Drbn)
%plot(DDrbn_cur)
%pause



end % variations

end % methods

end % classdef


