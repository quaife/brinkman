function [eulerX,eulerY,u,v] = computeEulerVelocity(...
    options,prams,radii,centers)

om = monitor(options,prams);
% monitor object

[Ninner,Nouter,nv,Xinner,Xouter,sigmaInner,sigmaOuter] = ...
    om.loadGeometry(options.dataFile);
% load all the information about the geometry and density function

if ((prams.Ninner == Ninner) + ...
      (prams.Nouter == Nouter) + (prams.nv == nv))~=3
  message = 'Saved data does not match input parameters';
  om.writeMessage(message);
  message = 'I am stopping';
  om.writeMessage(message);
  return
end

innerGeom = bodies(Xinner);
outerGeom = bodies(Xouter);
% build objects for the inner and outer boundaries

op = poten(innerGeom,options.fmm,false,true);
% object for evaluating layer potentials

xmin = prams.xmin; xmax = prams.xmax; nx = prams.nx;
ymin = prams.ymin; ymax = prams.ymax; ny = prams.ny;

[eulerX,eulerY] = meshgrid(...
    linspace(xmin,xmax,nx),linspace(ymin,ymax,ny));
% build the Eulerian grid where the solution will be computed

eX = eulerX(:); eY = eulerY(:);
% write as column vectors

vel = zeros(2*numel(eX),1);

tic
vel = op.layerEvalVelocity([eX;eY],...
      innerGeom,outerGeom,sigmaInner,sigmaOuter);
% compute velocity at points [eX;eY]

for k = 1:numel(eX)
  if any((eX(k) - centers(:,1)).^2 + ...
        (eY(k) - centers(:,2)).^2 <= radii.^2)
    vel(k) = 0;
    vel(k+numel(eX)) = 0;
  end
end
% assign a velocity of 0 to points that are inside a pore

u = reshape(vel(1:end/2),size(eulerX));
v = reshape(vel(end/2+1:end),size(eulerY));
% decompose velocity field in x and y coordinates
%
%om.writeStars
%message = '****   Velocity found on Eulerian Grid   ****';
%om.writeMessage(message);
%message = ['**** Required time was ' num2str(toc,'%4.2e') ...
%    ' seconds  ****'];
%om.writeMessage(message);
%om.writeStars
%om.writeMessage(' ');
%
%om.writeEulerVelocities(eulerX,eulerY,u,v);
%% save the velocity field on the Eulerian grid
%
%u = reshape(vel(1:end/2),size(eulerX));
%v = reshape(vel(end/2+1:end),size(eulerY));
%% put velocity field in meshgrid format


%eulerX = [];
%eulerY = [];
%u = [];
%v = [];
