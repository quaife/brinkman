addpath ../../src
addpath ..

if 1
%  file = 'shear2VesHData.bin.run1';
  file = '~/projects/brinkman/vesicle_code/results/shear2Ves/adR4em1adS2p1e0Chi5em1_ra090/shear2VesData.bin';
  shearRate = 0.5;
end
if 0
  file = 'shear4VesData.bin';
end

oc = curve;

RA = linspace(0.5,0.99,4);
ra = zeros(size(RA));
eff_visc = zeros(size(RA));

k = 1;

%for k = 1:numel(RA)
%  shear1Ves_RA(RA(k));

  [posx,posy,ten,~,~,ea,el,time,n,nv] = loadFile(file);
%  [posx,posy,ten,~,~,~,~,ea,el,time,n,nv] = loadFileOld(file);
  % load positions, tension, errors, time, number of points, and number of
  % vesicles

  kappa = 1; % bending stiffness
  adStrength = 2.1;
%  adStrength = 7e-1;
  adRange = 4e-1;
  adhesion = true;
  % adhesion strength, range, and a flag to denote that adhesion is active

  N = size(posx,1); % number of points per vesicle
  nv = size(posx,2); % number of vesicles
  ntime = size(posx,3); % number of time steps

  stress11 = zeros(nv,ntime);
  stress12 = zeros(nv,ntime);
  stress21 = zeros(nv,ntime);
  stress22 = zeros(nv,ntime);
  %stress11ver2 = zeros(nv,ntime);
  %stress12ver2 = zeros(nv,ntime);
  %stress21ver2 = zeros(nv,ntime);
  %stress22ver2 = zeros(nv,ntime);
  %
  [ra,area,length] = oc.geomProp([posx(:,:,1);posy(:,:,1)]);

  for k = 1:ntime
    X = [posx(:,:,k);posy(:,:,k)];
    [~,~,cur] = oc.diffProp(X);
    sigma = ten(:,:,k);
%    sigma = sigma + 3/2*cur.^2;
% don't need to do this (I think) because the intrinsic viscosity is
% defined in terms of whatever we call the bending force plus the
% corresponding 'tension', even if it isn't the phyiscal tension
    vesicle = capsules(X,sigma,[],kappa,[]);

    [tanx,tany] = oc.getXY(vesicle.xt);
    norx = +tany;
    nory = -tanx;

    f = -vesicle.tracJump(X,sigma)*0;
    if adhesion
      f = f - vesicle.adhesionTerm(adStrength,adRange);
    end
    [fx,fy] = oc.getXY(f);

    for j = 1:nv
      stress11(j,k) = sum(...
          posx(:,j,k).*fx(:,j).*vesicle.sa(:,j))*2*pi/N;
      stress12(j,k) = sum(...
          posx(:,j,k).*fy(:,j).*vesicle.sa(:,j))*2*pi/N;
      stress21(j,k) = sum(...
          posy(:,j,k).*fx(:,j).*vesicle.sa(:,j))*2*pi/N;
      stress22(j,k) = sum(...
          posy(:,j,k).*fy(:,j).*vesicle.sa(:,j))*2*pi/N;

%      stress11(j,k) = 1/area(j)*sum(...
%          posx(:,j,k).*fx(:,j).*vesicle.sa(:,j))*2*pi/N;
%      stress12(j,k) = 1/area(j)*sum(...
%          posx(:,j,k).*fy(:,j).*vesicle.sa(:,j))*2*pi/N;
%      stress21(j,k) = 1/area(j)*sum(...
%          posy(:,j,k).*fx(:,j).*vesicle.sa(:,j))*2*pi/N;
%      stress22(j,k) = 1/area(j)*sum(...
%          posy(:,j,k).*fy(:,j).*vesicle.sa(:,j))*2*pi/N;
%
  %    stress11ver2(j,k) = 1/area*sum((...
  %        kappa*vesicle.cur(:,j).^2.*norx(:,j).*norx(:,j) + ...
  %        sigma(:,j).*tanx(:,j).*tanx(:,j)).* ...
  %        vesicle.sa(:,j))*2*pi/N;
  %    stress12ver2(j,k) = 1/area*sum((...
  %        kappa*vesicle.cur(:,j).^2.*norx(:,j).*nory(:,j) + ...
  %        sigma(:,j).*tanx(:,j).*tany(:,j)).* ...
  %        vesicle.sa(:,j))*2*pi/N;
  %    stress21ver2(j,k) = 1/area*sum((...
  %        kappa*vesicle.cur(:,j).^2.*nory(:,j).*norx(:,j) + ...
  %        sigma(:,j).*tany(:,j).*tanx(:,j)).* ...
  %        vesicle.sa(:,j))*2*pi/N;
  %    stress22ver2(j,k) = 1/area*sum((...
  %        kappa*vesicle.cur(:,j).^2.*nory(:,j).*nory(:,j) + ...
  %        sigma(:,j).*tany(:,j).*tany(:,j)).* ...
  %        vesicle.sa(:,j))*2*pi/N;
    end
  end

  time1 = 0; time2 = 100;

 time1 = 37; time2 = 100 
% for RA \approx 0.90 and Hamaker constant = 2.1

% Everything below here uses Hamaker constant 0.7
%  time1 = 34; time2 = 97; % for RA \approx 0.99
%  time1 = 34.5; time2 = 97.5; % for RA \approx 0.95
%  time1 = 35.5; time2 = 99; % for RA \approx 0.90
%  time1 = 23; time2 = 89; % for RA \approx 0.85
%  time1 = 24; time2 = 92; % for RA \approx 0.80
%  time1 = 26; time2 = 95; % for RA \approx 0.75
%  time1 = 25; time2 = 97; % for RA \approx 0.70
%  time1 = 23; time2 = 99; % for RA \approx 0.65
%  time1 = 32; time2 = 87; % for RA \approx 0.60
%  time1 = 35; time2 = 95.5; % for RA \approx 0.55
%  time1 = 16; time2 = 92; % for RA \approx 0.50
%  time1 = 40; time2 = 100; % for RA \approx 0.45
%  time1 = 38.5; time2 = 96; % for RA \approx 0.40
  [~,s1] = min(abs(time - time1));
  [~,s2] = min(abs(time - time2));

  eff_visc = 0;
  for k = 1:nv
    eff_visc = eff_visc + sum(diff(time(s1:s2))'.*stress12(k,s1+1:s2));
  end
  totalArea = area(1) + area(2);
  eff_visc = eff_visc/(time(s2) - time(s1))/shearRate/totalArea;
  clf
  plot(time(s1:s2),stress12(:,s1:s2))
%  eff_visc = 1/(time(end) - time(s))/shearRate*...
%      sum(diff(time(s:end))'.*stress12(s+1:end));

%end






