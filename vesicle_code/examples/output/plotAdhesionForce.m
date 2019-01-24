addpath ../../src
addpath ..


if 1
  file = '~/projects/brinkman/vesicle_code/results/extensional2Ves/adR4em1adS7em1Chi7em2_ra070/extensional2VesData.bin';
end

oc = curve;
%
[posx,posy] = loadFile(file);
% load positions, tension, errors, time, number of points, and number of
% vesicles

adRange = 0.4;
adStrength = 0.7;
% adhesion range and strength
N = size(posx,1); % number of points per vesicle
nv = size(posx,2); % number of vesicles
ntime = size(posx,3); % number of time steps

m = size(posx,3);
% final time step
xx = posx(:,:,m);
yy = posy(:,:,m);
X = [xx;yy];
% vesicle position

vesicle = capsules(X,[],[],[],[]);
% create object for the vesicle
f = vesicle.adhesionTerm(adStrength,adRange);
% adhesion force of the vesicles
[fx,fy] = oc.getXY(f);
% seperate adhesion force into x and y coordinates

figure(1); clf
plot(xx,yy,'r')
hold on;

scale = 1/6;
% adjust if the quiver plot if the quiver plot has arrows that are too
% big or small
quiver(xx,yy,fx*scale,fy*scale,0,'b')
axis equal;





