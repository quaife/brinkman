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
function o = capsules(X,rcon,params,varargin)
% constructor for capsules class

o.N = size(X,1)/2; % shorthand for Number of discretization points 
oc = curve(o.N); %shorthand for curve class
o.X = X; % shorthand for discretization points
%compute the length, theta, and curvature of the vesicle
[o.L,o.theta,o.cur] = oc.computeOpeningAngle(o.N,X); 
o.rcon = rcon; % shorthand for concentration
o.x0 = X(1); % shorthand for x tracking point
o.y0 = X(o.N + 1); % shorthand for y tracking point

if nargin >= 3
    o.bendsti = params.bendsti;
    o.bendratio = params.bendratio;
    o.viscIn = params.viscosityInside;
    o.viscOut = params.viscosityOutside;
    o.SPc = params.SPcoeff;
end
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

N = ves.N;
oc = curve(N);
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
% Adjust the shape until error in area is less than the tolerance
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

N = ves.N;
oc = curve(N);
IK = oc.modes(N);
rcon = ves.rcon;
save('rcon.mat', 'rcon')
cur = oc.acurv(ves.N,ves.theta,ves.L);
b0 = ves.bendsti;
b1 = ves.bendsti * ves.bendratio;
% bending coefficient which depends on the lipid concentration that is
% stored in rcon. This is the variable b(u) in equation (10)
% rcon is the concentration u
rbn = b0 * (ones(N,1) - rcon) + b1*rcon;
% take the derivative of b(u)
Drbn = oc.diffFT(rbn)/ves.L;
% Fourier modes
% !!! IK = 2*pi*1i*[0:1:ves.N/2 -ves.N/2+1:1:-1]';
% derivative of the curvature
Drbn_cur = oc.diffFT(rbn.*cur)/ves.L; 
% second derivative of the curvature
DDrbn_cur = oc.diffFT(Drbn_cur)/ves.L; 
%Esigma is equation (14) with spotaneous curvature set to zero.
Esigma = -DDrbn_cur - 0.5*rbn.*cur.^3;
%Eu is the second term in equation (13) (differs by a negative
%sign - possibly from the negative sign in eq(23) which has a negative on 
%the variation for u). The last term drops since spontaneous curvature is 
%0. The first term is not in this routine since we are only calculating 
%variations due to changes in the vesicle shape and not the lipid species.
Eu = -0.5*Drbn.*cur.^2;
save('drbnCur.mat', 'Drbn','cur')
%** SHUWANG QUESTION: THIS IS A PLUS SIGN IN THE PAPER (EQUATION (13)),
% BUT IS A MINUS SIGN IN SHUWANG'S CODE **OLD COMMENT???
% ADDED -

end % variations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NearSelf,NearOther] = getZone(vesicle1,vesicle2,relate)
% [NearSelf,NearOther] = getZone(vesicle1,vesicle2,relate) constructs
% each vesicle, index of the closest point, nearest point on a local
% interpolant, and argument of that nearest point.  vesicle1 contains
% the source points (which are also target points) and vesicle2 
% contains additional target points.  The 
% values of relate corresond to
% relate == 1  => only require NearSelf  (ie. vesicle to vesicle)
% relate == 2  => only require NearOther (ie. vesicle to wall)
% relate == 3  => require both NearSelf and NearOther
% THIS ROUTINE HAS A LOOP OVER THE TARGET POINTS WHICH SHOULD BE REMOVED
NearSelf = [];
NearOther = [];

N1 = vesicle1.N; % number of source/target points
%nv1 = vesicle1.nv; % number of source/target vesicles
nv1 = 1; % number of source/target vesicles
X1 = vesicle1.X; % source and target points
oc = curve(N1);
[xsou,ysou] = oc.getXY(X1); 
% separate targets into x and y coordinates

h = vesicle1.L/N1; 
% smallest arclength over all vesicles
ptsperbox = 10; 
% Estimate for number of points per box.  This simply sets the 
% number of uniformly refined boxes we take.  Estimate is not very
% accurate.  What ptsperbox represents is the total number of points
% that could be put in each two-dimensional bin where no two are
% less than distance h from one another.  However, our points live
% on curves and thus will not fill up an entire bin

H = sqrt(ptsperbox)*h;
xmin = min(min(xsou));
xmax = max(max(xsou));
xmin = xmin - H;
xmax = xmax + H;
ymin = min(min(ysou));
ymax = max(max(ysou));
ymin = ymin - H;
ymax = ymax + H;
% Add a buffer around the points so that it is easier to
% work with vesicle2

Nx = ceil((xmax - xmin)/H);
Ny = ceil((ymax - ymin)/H);
% Find bounds for box that contains all points and add a buffer 
% so that all points are guaranteed to be in the box

Nbins = Nx * Ny; % Total number of bins

ii = ceil((xsou - xmin)/H);
jj = ceil((ysou - ymin)/H);
% Index in x and y direction of the box containing each point
bin = (jj-1)*Nx + ii;
% Find bin of each point using lexiographic ordering (x then y)

%figure(2);
%clf; hold on
%plot(xsou,ysou,'k.')
%axis equal
%axis([xmin xmin+Nx*H ymin ymin+Ny*H])
%set(gca,'xtick',linspace(xmin,xmin+Nx*H,Nx+1))
%set(gca,'ytick',linspace(ymin,ymin+Ny*H,Ny+1))
%grid on
%set(gca,'xticklabel',[])
%set(gca,'yticklabel',[])
%figure(1)
%pause
% DEBUG: This does a simple plot of the points with a grid that 
% aligns with the boundary of the boxes

fpt = zeros(Nbins,nv1);
lpt = zeros(Nbins,nv1);
% allocate space for storing first and last points
[binsort,permute] = sort(bin);
% build permute.  Need binsort to find first and last points
% in each box

for k = 1:nv1 % Loop over vesicles
  for j = 1:N1 % Loop over bins
    ibox = binsort(j,k);
    if (fpt(ibox,k) == 0)
      fpt(ibox,k) = j;
      lpt(ibox,k) = 1;
    else
      lpt(ibox,k) = lpt(ibox,k) + 1;
    end
  end
  lpt(:,k) = fpt(:,k) + lpt(:,k) - 1;
end
% Construct first and last point in each box corresponding
% to each vesicle.  The order is based on permute.  For example,
% permute(fpt(ibox,k)),...,permute(lpt(ibox,k)) is the set of 
% all points from vesicle k contained in box ibox

neigh = zeros(Nbins,9);

%Do corners first
neigh(1,1:4) = [1 2 Nx+1 Nx+2]; 
% bottom left corner
neigh(Nx,1:4) = [Nx Nx-1 2*Nx 2*Nx-1]; 
% bottom right corner
neigh(Nbins-Nx+1,1:4) = [Nbins-Nx+1 Nbins-Nx+2 ...
    Nbins-2*Nx+1 Nbins-2*Nx+2];
% top left corner
neigh(Nbins,1:4) = [Nbins Nbins-1 Nbins-Nx Nbins-Nx-1]; 
% top right corner

for j = 2:Nx-1
  neigh(j,1:6) = j + [-1 0 1 Nx-1 Nx Nx+1];
end
% neighbors of bottom row

for j = Nbins-Nx+2:Nbins-1
  neigh(j,1:6) = j + [-1 0 1 -Nx-1 -Nx -Nx+1];
end
% neighbors of top row

for j=Nx+1:Nx:Nbins-2*Nx+1
  neigh(j,1:6) = j + [-Nx -Nx+1 0 1 Nx Nx+1];
end
% neighbors of left column

for j=2*Nx:Nx:Nbins-Nx
  neigh(j,1:6) = j + [-Nx-1 -Nx -1 0 Nx-1 Nx];
end
% neighbors of right column

J = (Nx + 1:Nbins - Nx);
J = J(mod(J-1,Nx)~=0);
J = J(mod(J,Nx)~=0);
% J is the index of boxes that are not on the boundary
for j=J
  neigh(j,:) = j + [-Nx-1 -Nx -Nx+1 -1 0 1 Nx-1 Nx Nx+1];
end
% neighbors of interior points
% TREE STRUCTURE IS COMPLETE


if (relate == 1 || relate == 3)
  for k = 1:nv1
    distSS{k} = spalloc(N1,nv1,0);
    % dist(n,k,j) is the distance of point n on vesicle k to
    zoneSS{k} = spalloc(N1,nv1,0);
    % near or far zone
    nearestSS{k} = spalloc(2*N1,nv1,0);
    % nearest point using local interpolant
    icpSS{k} = spalloc(N1,nv1,0);
    % index of closest discretization point
    argnearSS{k} = spalloc(N1,nv1,0);
    % argument in [0,1] of local interpolant
  end
  % Represent near-singular integration structure so that we can use
  % sparse matricies.

  % begin classifying points where we are considering 
  % vesicle to vesicle relationships
  for k = 1:nv1
    boxes = unique(bin(:,k));
    % Find all boxes containing points of vesicle k
    boxes = neigh(boxes,:);
    % Look at all neighbors of boxes containing vesicle k
    boxes = unique(boxes(:));
    % Remove repetition
    boxes = boxes(boxes~=0);
    % Delete non-existent boxes that came up because of neigh

    K = [(1:k-1) (k+1:nv1)];
    for k2 = K
      istart = fpt(boxes,k2);
      iend = lpt(boxes,k2);
      istart = istart(istart ~= 0);
      iend = iend(iend ~= -1);
      % Find index of all points in neighboring boxes of vesicle k
      % that are in vesicle k2
      
      neighpts = zeros(sum(iend-istart+1),1);
      
      % Allocate space to assign possible near points
      is = 1;
      for j=1:numel(istart)
        ie = is + iend(j) - istart(j);
        neighpts(is:ie) = permute(istart(j):iend(j),k2);
        is = ie + 1;
      end
      % neighpts contains all points on vesicle k2 that are in 
      % neighboring boxes to vesicle k

      neighpts = sort(neighpts);
      % sorting should help speedup as we won't be jumping around
      % through different boxes

      n = 0;
      for i=1:numel(neighpts)
        ipt = neighpts(i);
        ibox = bin(ipt,k2);
        % box containing ipt on vesicle k2
        if (ibox ~= n)
          n = ibox;
          % Check if we need to move to a new box
          neighbors = neigh(ibox,:);
          % neighbors of this box
          neighbors = neighbors(neighbors~=0);
          % Remove non-existent neighbors
          istart = fpt(neighbors,k);
          iend = lpt(neighbors,k);
          istart = istart(istart ~= 0);
          iend = iend(iend ~= -1);
          % Find points on vesicle k in neighboring boxes
          neighpts2 = zeros(sum(iend-istart+1),1);
          is = 1;
          for j=1:numel(istart)
            ie = is + iend(j) - istart(j);
            neighpts2(is:ie) = permute(istart(j):iend(j),k);
            is = ie + 1;
          end
          % neighpts2 contains all points on vesicle k that 
          % are in neighboring box of ibox 
        end % decide if we need to switch boxes

        [d0,d0loc] = min((xsou(ipt,k2) - xsou(:,k)).^2 + ...
            (ysou(ipt,k2) - ysou(:,k)).^2);
        % Find minimum distance between ipt on vesicle k2 to
        % possible closest points on vesicle k
        d0 = sqrt(d0);
        % Save on not taking the square root on a vector but instead
        % on a single real number

        icpSS{k}(ipt,k2) = d0loc;
        if (d0 < 2*h);
          [distSS{k}(ipt,k2),nearestx,nearesty,argnearSS{k}(ipt,k2)] = ...
              vesicle1.closestPnt([xsou;ysou],xsou(ipt,k2),...
              ysou(ipt,k2),k,icpSS{k}(ipt,k2));
          nearestSS{k}(ipt,k2) = nearestx;
          nearestSS{k}(ipt+N1,k2) = nearesty;
          % Find closest point along a local interpolant using
          % Newton's method.

          if (distSS{k}(ipt,k2) < h)
            zoneSS{k}(ipt,k2) = 1;
          end
          % Point ipt of vesicle k2 is in the near zone of
          % vesicle k
        end
      end % ipt
    end % k2
  end % k

  NearSelf.dist = distSS;
  NearSelf.zone = zoneSS;
  NearSelf.nearest = nearestSS;
  NearSelf.icp = icpSS;
  NearSelf.argnear = argnearSS;
  % Store everything in the structure NearSelf.  This way it is 
  % much cleaner to pass everything around
end % relate == 1 || relate == 3

% Bin target points with respect to the source points
if (relate == 2 || relate == 3)
  N2 = vesicle2.N; % number of additional targets
%  nv2 = vesicle2.nv; % number of additional vesicles
  nv2 = 1; % number of additional vesicles
  X2 = vesicle2.X; % additional target points
  [xtar,ytar] = oc.getXY(X2);
  % separate additional target points into x and y coordinates
%  figure(2); clf
%  plot(xtar,ytar,'r.')
%  hold on
%  plot(xsou,ysou,'k.')
%  pause
% DEBUG: FOR SEEING TARGET AND SOURCE POINTS IN THE TREE STRUCTURE
% WHICH CAN BE PLOTTED ABOVE

  for k = 1:nv1
    distST{k} = spalloc(N1,nv2,0);
    % dist(n,k,j) is the distance of point n on vesicle k to
    zoneST{k} = spalloc(N1,nv2,0);
    % near or far zone
    nearestST{k} = spalloc(2*N1,nv2,0);
    % nearest point using local interpolant
    icpST{k} = spalloc(N1,nv2,0);
    % index of closest discretization point
    argnearST{k} = spalloc(N1,nv2,0);
    % argument in [0,1] of local interpolant
  end
  % Represent near-singular integration structure using sparse matricies

  itar = ceil((xtar - xmin)/H);
  jtar = ceil((ytar - ymin)/H);
  [indx,indy] = find((itar >= 1) & (itar <= Nx) & ...
      (jtar >= 1) & (jtar <= Ny));
  % Only have to consider xx(ind),yy(ind) since all other points
  % are not contained in the box [xmin xmax] x [ymin ymax]

  for k = 1:nv1 % loop over sources
    for nind = 1:numel(indx) 
      % loop over points that are not outside the box that surrounds
      % all target points with a sufficiently large buffer
      ii = indx(nind);
      jj = indy(nind);
      binTar = (jtar(ii,jj)-1)*Nx + itar(ii,jj);
      boxesTar = neigh(binTar,:);
      boxesTar = boxesTar(boxesTar~=0);
      istart = fpt(boxesTar,k);
      iend  = lpt(boxesTar,k);
      istart = istart(istart ~= 0);
      iend   = iend(iend ~= -1);
      
      neighpts = zeros(sum(iend-istart+1),1);
  
      % Allocate space to assign possible near points
      if numel(neighpts) > 0
        % it is possible of the neighboring boxes to contain
        % no points.
        is = 1;
        
        for j = 1:numel(istart)
          ie = is + iend(j) - istart(j);
          neighpts(is:ie) = permute(istart(j):iend(j),k);
          is = ie + 1;
        end
        % Set of potentially nearest points to 
        % (xtar(jj),ytar(jj))
        
        [d0,d0loc] = min((xtar(ii,jj) - xsou(neighpts,k)).^2 + ...
          (ytar(ii,jj) - ysou(neighpts,k)).^2);
        % find closest point and distance between (xtar(jj),ytar(jj))
        % and vesicle k.  Only need to look at points in neighboring
        % boxes
        
        d0 = d0.^0.5;
        icpST{k}(ii,jj) = neighpts(d0loc);

        if d0 < 2*h
          [distST{k}(ii,jj),nearestx,nearesty,argnearST{k}(ii,jj)] = ...
            vesicle1.closestPnt([xsou;ysou],xtar(ii,jj),...
            ytar(ii,jj),k,icpST{k}(ii,jj));
          nearestST{k}(ii,jj) = nearestx;
          nearestST{k}(ii+N2,jj) = nearesty;
          
          % DEBUG: CHECK THAT NEWTON'S METHOD HAS DONE A GOOD JOB
          % CONVERGING TO THE NEAREST POINT
          % compute distance and nearest point between 
          % (xtar(ii,jj),ytar(ii,jj)) and vesicle k
          if distST{k}(ii,jj) < h
            zoneST{k}(ii,jj) = 1;
            % (xtar(ii,jj),ytar(ii,jj)) is in the near zone of vesicle k
          end
        end % d0 < 2*h
      end % numel(neighpts) > 0
    end % ii and jj
  end % k

  NearOther.dist = distST;
  NearOther.zone = zoneST;
  NearOther.nearest = nearestST;
  NearOther.icp = icpST;
  NearOther.argnear = argnearST;
  % store near-singluar integration requirements in structure NearOther
end % relate == 2 || relate == 3

end % getZone

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dist,nearestx,nearesty,theta] = ...
    closestPnt(o,X,xtar,ytar,k,icp)
% [dist,nearestx,nearesty,theta] = closestPnt(X,xtar,ytar,k,icp)
% computes the closest point on vesicle k to (xtar,ytar)
% using a Lagrange interpolant.  icp is the index of the closest
% point on the discrete mesh which is used as an initial guess

N = size(X,1)/2; % Number of points per vesicle
A = poten.lagrangeInterp;
interpOrder = size(A,1);
% need interpolation matrix and its size

p = ceil((interpOrder+1)/2);
% Accommodate for either an even or odd number of interpolation points
pn = mod((icp-p+1:icp-p+interpOrder)' - 1,N) + 1;
% band of points around icp.  The -1,+1 combination sets index
% 0 to N as required by the code

px = A*X(pn,k); % polynomial interpolant of x-coordinate
py = A*X(pn+N,k); % polynomial interpolant of y-coordinate
Dpx = px(1:end-1).*(interpOrder-1:-1:1)';
Dpy = py(1:end-1).*(interpOrder-1:-1:1)';
D2px = Dpx(1:end-1).*(interpOrder-2:-1:1)';
D2py = Dpy(1:end-1).*(interpOrder-2:-1:1)';
% To do Newton's method, need two derivatives

theta = 1/2;
% midpoint is a good initial guess
for newton = 1:1
  zx = filter(1,[1 -theta],px);
  zx = zx(end);
  zy = filter(1,[1 -theta],py);
  zy = zy(end);
  Dzx = filter(1,[1 -theta],Dpx);
  Dzx = Dzx(end);
  Dzy = filter(1,[1 -theta],Dpy);
  Dzy = Dzy(end);
  D2zx = filter(1,[1 -theta],D2px);
  D2zx = D2zx(end);
  D2zy = filter(1,[1 -theta],D2py);
  D2zy = D2zy(end);
  % Using filter is the same as polyval, but it is much
  % faster when only requiring a single polyval such as here.

  newtonNum = (zx-xtar)*Dzx + (zy-ytar)*Dzy;
  % numerator of Newton's method
  newtonDen = (zx-xtar)*D2zx + (zy-ytar)*D2zy + ...
      Dzx^2 + Dzy^2;
  % denominator of Newton's method
  theta = theta - newtonNum/newtonDen;
  % one step of Newton's method
end
% Do a few (no more than 3) Newton iterations

nearestx = filter(1,[1,-theta],px);
nearestx = nearestx(end);
nearesty = filter(1,[1,-theta],py);
nearesty = nearesty(end);
dist = sqrt((nearestx - xtar)^2 + (nearesty - ytar)^2);
% Compute nearest point and its distance from the target point

end % closestPnt




end % methods

end % classdef


