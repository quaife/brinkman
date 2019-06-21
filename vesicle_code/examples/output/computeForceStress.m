addpath ../../src

if 0
  file = 'relaxation2VesData.bin';
end
if 0
  file = 'shear1VesData.bin';
end
if 0
  file = 'shear2VesData.bin';
  ax = [-8 8 -3 3];
end
if 0
  file = 'shear4VesData.bin';
end
if 1
  file = '~/projects/brinkman/vesicle_code/results/extensional2Ves/adR4em1adS7em1Chi1em1_ra070/extensional2VesData.bin';
end

oc = curve;

[posx,posy,ten,~,~,ea,el,time,n,nv] = loadFile(file);
% load positions, tension, errors, time, number of points, and number of
% vesicles

adRange = 0.4;
adStrength = 0.7;
adhesion = true;
% adhesion strength, range, and a flag to denote that adhesion is active

N = size(posx,1); % number of points per vesicle
nv = size(posx,2); % number of vesicles
ntime = size(posx,3); % number of time steps

Forcex = zeros(nv,ntime);
Forcey = zeros(nv,ntime);
%Torque = zeros(nv,ntime);

for k = 1:ntime
  X = [posx(:,:,k);posy(:,:,k)];
  sigma = ten(:,:,k);
  vesicle = capsules(X,sigma,[],1e-1,[]);
%  vesicle.center = 0*vesicle.center + 0;
%  vesicle.center = 2*ones(2,nv);

  tanx = vesicle.xt(1:end/2,:);
  tany = vesicle.xt(end/2+1:end,:);
  norx = +tany;
  nory = -tanx;

%  f = vesicle.tensionTerm(sigma);
%  f = vesicle.bendingTerm(X);
  f = vesicle.tracJump(X,sigma);
%  f = zeros(2*N,nv);
  if adhesion
    f = f + vesicle.adhesionTerm(adStrength,adRange);
  end

  Forcex(:,k) = sum(f(1:end/2,:).*vesicle.sa)*2*pi/N;
  Forcey(:,k) = sum(f(end/2+1:end,:).*vesicle.sa)*2*pi/N;

%  [cx,cy] = oc.getXY(vesicle.center);
%
%  for j=1:nv
%    Torque(j,k) = ...
%    sum((f(1:end/2,j).*(X(end/2+1:end,j) - cy(j)) - ...
%         f(end/2+1:end,j).*(X(1:end/2,j) - cx(j))).*...
%         vesicle.sa(:,j))*2*pi/N;
%  end
%
%  figure(1); clf;
%  subplot(2,2,1) 
%  plot(posx(:,1,k),posy(:,1,k),'b')
%  hold on
%  plot(posx(:,2,k),posy(:,2,k),'r')
%  hold on;
%  plot(cx,cy,'k.','markersize',15);
%  axis equal
%  axis(ax)
%
%%  subplot(2,2,3)
%%  plot(time(1:k),Torque(1,1:k),'b')
%%  hold on;
%%  plot(time(1:k),Torque(2,1:k),'r')
%%  xlim([0 time(end)])
%%%  ylim([min(Torque(1,:)) max(Torque(1,:))])
%%  ylabel('Torque')
%
%  subplot(2,2,2)
%  plot(time(1:k),Forcex(1,1:k),'b')
%  hold on;
%  plot(time(1:k),Forcex(2,1:k),'r')
%  xlim([0 time(end)])
%%  ylim([-max(abs(Forcex(1,:))) max(abs(Forcex(1,:)))])
%  ylabel('x-Force')
%
%  subplot(2,2,4)
%  plot(time(1:k),Forcey(1,1:k),'b')
%  hold on;
%  plot(time(1:k),Forcey(2,1:k),'r')
%  xlim([0 time(end)])
%%  ylim([min(Forcey(1,:)) max(Forcey(1,:))])
%  ylabel('y-Force')
%
%  pause(0.01);
end




