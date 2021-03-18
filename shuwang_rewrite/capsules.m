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
SPc;
ten;

end % properties

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = capsules(X,rcon,params)
% constructor for capsules class

oc = curve; %shorthand for curve class
o.N = size(X,1)/2; % shorthand for Number of discretization points 
o.X = X; % shorthand for discretization points
%compute the length, theta, and curvature of the vesicle
[o.L,o.theta,o.cur] = oc.computeOpeningAngle(o.N,X); 
o.rcon = rcon; % shorthand for concentration
o.x0 = X(1); % shorthand for x tracking point
o.y0 = X(o.N + 1); % shorthand for y tracking point

o.bendsti = params.bendsti;
o.bendratio = params.bendratio;
o.viscIn = params.viscosityInside;
o.viscOut = params.viscosityOutside;
o.SPc = params.SPcoeff;
%%% Curvature check:
    %disp('here')
    %curv = oc.acurv(o);
    %clf; hold on
    %plot(o.cur)
    %plot(curv,'r--')
    %norm(curv - o.cur,inf)
    %pause
%%%
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
% Changed this from 20
maxFreq = 4;
thetah(maxFreq:N-maxFreq) = 0;
% remove all modes with Frequency above 4. Note that Shuwang's original
% code had an error since he did not truncate the positive and negative
% coefficients to the same level. Therefore, the imaginary part after
% taking an ifft was not 0.
theta_periodicPart = real(ifft(thetah));
theta = theta_periodicPart + theta0 + 2*pi*(0:N-1)'/N;
% Reconstruct the vesicle shape
X = oc.recon(N,x0,y0,L,theta);
% new area with filtered opening angle
area = sum(sin(theta).*X(1:end/2) - cos(theta).*X(end/2+1:end))*...
      0.5*L/N;
%Adjust theshape until error in area is less than the tolerance
iter = 1;
while abs(area - areaRef)/areaRef > 1e-10
  theta_periodicPart = theta_periodicPart * ...
        (1 + (area - areaRef)/30);
  theta = theta_periodicPart + theta0 + 2*pi*(0:N-1)'/N;
  % new shape with scaled periodic part of theta
  X = oc.recon(N,x0,y0,L,theta);
  % new area
  area = sum(sin(theta).*X(1:end/2) - cos(theta).*X(end/2+1:end))*...
        0.5*L/N;
  % break if max iterations reached
  iter = iter + 1;
  if iter > 100
    break
  end
end
% update X
%X(1:end/2) = X(1:end/2) - mean(X(1:end/2));
%X(end/2+1:end) = X(end/2+1:end) - mean(X(end/2+1:end));
%ves.X = X;
% update x0 and y0
%ves.x0 = X(1);
%ves.y0 = X(1 + ves.N);
% replace the geometry with the new shape

[ves.L,ves.theta,ves.cur] = oc.computeOpeningAngle(N,X);

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
cur = oc.acurv(ves.N,ves.theta,ves.L);
b0 = ves.bendsti;
b1 = ves.bendsti * ves.bendratio;
% bending coefficient which depends on the lipid concentration that is
% stored in rcon. This is the variable b(u) in equation (10)
% rcon is the concentration u
rbn = b0 * (ones(N,1) - rcon) + b1*rcon;
% take the derivative of b(u)
Drbn = oc.diffFT(rbn,IK)/ves.L;
% Fourier modes
IK = 2*pi*1i*[0:1:ves.N/2 -ves.N/2+1:1:-1]';
% derivative of the curvature
Drbn_cur = oc.diffFT(rbn.*cur,IK)/ves.L; 
% second derivative of the curvature
DDrbn_cur = oc.diffFT(Drbn_cur,IK)/ves.L; 
%Esigma is equation (14) with spotaneous curvature set to zero.
Esigma = -DDrbn_cur - 0.5*rbn.*cur.^3;
%Eu is the second term in equation (13) (differs by a negative
%sign - possibly from the negative sign in eq(23) which has a negative on 
%the variation for u). The last term drops since spontaneous curvature is 
%0. The first term is not in this routine since we are only calculating 
%variations due to changes in the vesicle shape and not the lipid species.
Eu = -0.5*Drbn.*cur.^2;
%** SHUWANG QUESTION: THIS IS A PLUS SIGN IN THE PAPER (EQUATION (13)),
% BUT IS A MINUS SIGN IN SHUWANG'S CODE **OLD COMMENT???
% ADDED -

end % variations

end % methods

end % classdef


