classdef curve
% This class implements that basic calculus on the curve.
% The basic data structure is a matrix X in which the columns 
% represent periodic C^{\infty} closed curves with N points, 
% X(1:n,j) is the x-coordinate of the j_th curve and X(n+1:N,j) 
% is the y-coordinate of the j_th curve; here n=N/2
% X coordinates do not have to be periodic, but the curvature,
% normals, etc that they compute will be garbage.  This allows
% us to store target points for tracers or the pressure using
% this class and then using near-singular integration is easy
% to implement

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y] = getXY(o,X)
% [x,y] = getXY(X) get the [x,y] component of curves X
N = size(X,1)/2;
x = X(1:N,:);
y = X(N+1:end,:);

end % getXY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function V = setXY(o,x,y)
% V = setXY(x,y) set the [x,y] component of vector V on the curve
N = size(x,1);
V=zeros(2*N,size(x,2));
V(1:N,:) = x;
V(N+1:end,:) = y;

end % setXY


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Dx,Dy] = getDXY(o,X)
% [Dx,Dy] = getDXY(X), compute the derivatives of each component of X
% these are the derivatives with respect the parameterization not
% arclength
x = X(1:end/2,:);
y = X(end/2+1:end,:);
N = size(x,1);
nv = size(x,2);
IK = fft1.modes(N,nv);
Dx = fft1.diffFT(x,IK);
Dy = fft1.diffFT(y,IK);

end % getDXY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [jacobian,tangent,curvature] = diffProp(o,X)
% [jacobian,tangent,curvature] = diffProp(X) returns differential
% properties of the curve each column of the matrix X. Each column of 
% X should be a closed curve defined in plane. The tangent is the 
% normalized tangent.
%
% EXAMPLE:
%    n = 128; nv = 3;
%    X = boundary(n,'nv',nv,'curly');
%    c = curve;
%    [k t s] = c.diffProp(X);

N = size(X,1)/2;
nv = size(X,2);

% get the x y components
[Dx,Dy] = o.getDXY(X);

jacobian = sqrt(Dx.^2 + Dy.^2); 

if nargout>1  % if user requires tangent
  tangent = o.setXY(Dx./jacobian,Dy./jacobian);
end

if nargout>2  % if user requires curvature
  IK = fft1.modes(N,nv);
  DDx = curve.arcDeriv(Dx,1,ones(N,nv),IK);
  DDy = curve.arcDeriv(Dy,1,ones(N,nv),IK);
  curvature = (Dx.*DDy - Dy.*DDx)./(jacobian.^3);
end
% curvature of the curve

end % diffProp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [reducedArea,area,length] = geomProp(o,X)
% [reducedArea area length] = geomProp(X) calculate the length, area 
% and the reduced volume of domains inclose by columns of X. 
% Reduced volume is defined as 4*pi*A/L. 
% EXAMPLE:
%   X = boundary(64,'nv',3,'curly');
%   c = curve(X);
%   [rv A L] = c.geomProp(X);

if size(X,1) < 96
  X = [interpft(X(1:end/2,:),96);interpft(X(end/2+1:end,:),96)];
end
% if there are less than 96 points on the boundary, upsample to 96

[x,y] = o.getXY(X);
N = size(x,1);
[Dx,Dy] = o.getDXY(X);
length = sum(sqrt(Dx.^2 + Dy.^2))*2*pi/N;
area = sum(x.*Dy - y.*Dx)*pi/N;

reducedArea = 4*pi*area./length.^2;

end % geomProp



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,nv] = initConfig(o,N,varargin)       
% [X,nv] = initConfig(n,varargin) returns N coordinates of boundary
% points.
% X = BOUNDARY(N,OPTIONS) can be used to call for different
% configuration.  Available OPTIONS are
%
%   'nv'          - (followed by) number of vesicles (in 
%                   lineajr configuration),
%   'angle'       - inclination angle of the vesicle(s) form the vertical
%                   position, 
%   'curly'       - returns a somehow wiggly vesicle.
%   'couette'     - the vesicles inside the default geometry of 
%                   couette flow (confined). 
%   'scale'       - multiplies the size of the output boundary
%   'choke'       - returns a choked domain.
%   'slit'        - returns a slitted geometry
%   'couette'     - returns a domain for a couette apparatus.
%   'tube'        - returns a domain that is an ellongated ellipse

% EXAMPLE: X = boundary(64,'nv',3,'theta',pi/6);
%

options = varargin;
X = [];

if(any(strcmp(options,'nv')))
  nv = options{find(strcmp(options,'nv'))+1};
else
  nv = 1;
end
% number of vesicles

if(any(strcmp(options,'angle')))
  theta = options{find(strcmp(options,'angle'))+1};
else
  theta = zeros(nv,1);
end
% rotate the vesicle by a certain angle

if(any(strcmp(options,'center')))
  cen = options{find(strcmp(options,'center'))+1};
else
  cen = [(0:nv-1);zeros(1,nv)];
end
% pick the center of the vesicles

if(any(strcmp(options,'reducedArea')))
  ra = options{find(strcmp(options,'reducedArea'))+1};
else 
  ra = [];
end
% desired reduced area

if(any(strcmp(options,'scale')))
  scale = options{find(strcmp(options,'scale'))+1};
else
  scale = 1;
end
% scale the size of the vesicle (reduced area is invariant under
% scale)

if(any(strcmp(options,'folds')))
  folds = options{find(strcmp(options,'folds'))+1};
else
  folds = 1;
end
% number of folds for a star shape

t = (0:N-1)'*2*pi/N;
% Discretization in parameter space

if any(strcmp(options,'curly'))
  a = 1; b = 3*a; c = 0.85; 
  r = 0.5*sqrt( (a*cos(t)).^2 + (b*sin(t)).^2) + ...
      .07*cos(12*(t));
  x = scale*c*r.*cos(t);
  y = scale*r.*sin(t);
  X0 = [x;y];
  % radius of curly vesicle

elseif any(strcmp(options,'star'))
  radius = 1 + 0.4*cos(folds*t);
  X = scale*[radius.*cos(t);radius.*sin(t)];
  % a star that comes very close to intersecting itself at the origin

%  sa = o.diffProp(X);
%  t = o.arc(sa);
%  % Find discretization points to achieve equispaced in arclength
%  % discretization points
%  radius = 1 + 0.4*cos(folds*t);
%  X = scale*[radius.*cos(t);radius.*sin(t)];
%  % Reparameterize with points that are nearly equispaced in arclength

elseif any(strcmp(options,'choke'))
  a = 10; b = 3; c = 0.6; order = 8;
  % parameters for the boundary
  Nsides = ceil(0.5*b/(2*a+2*b)*N);
  Ntop = (N-4*Nsides)/2;
  t1 = linspace(0,0.2*pi,Nsides+1); t1 = t1(1:end-1)';
  t2 = linspace(0.2*pi,pi-0.2*pi,Ntop+1); t2 = t2(1:end-1)';
  t3 = linspace(pi-0.2*pi,pi+0.2*pi,2*Nsides+1); t3 = t3(1:end-1)';
  t4 = linspace(pi+0.2*pi,2*pi-0.2*pi,Ntop+1); t4 = t4(1:end-1)';
  t5 = linspace(2*pi-0.2*pi,2*pi,Nsides+1); t5 = t5(1:end-1)';
  t = [t1;t2;t3;t4;t5];
  % Parameterize t so that geometry is closer to equi-spaced in
  % arclength
%  t = (0:N-1)'*2*pi/N;
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*cos(t); y = b*r.*sin(t);
  ind = abs(x) < pi;
  y(ind) = y(ind).*(1-c*cos(x(ind)))/(1+c);
  X0 = [x;y];
  % choked domain.  a and b control the length and height.  c
  % controls the width of the gap, and order controls the
  % regularity

elseif any(strcmp(options,'chokeLong'))
  a = 20; b = 3; c = 0.6; order = 8;
  % parameters for the boundary
  t = (0:N-1)'*2*pi/N;
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*cos(t); y = b*r.*sin(t);
  ind = abs(x) < 10;
%  ind = abs(x) < 0;
  scalRat = 2*c/(1+c)*(0.5-0.5*cos(pi*x(ind(1:end/2))/10)).^10 + ...   
      (1-c)/(1+c);
  y(ind) = y(ind).*[scalRat;scalRat];
  X0 = [x;y];
  % choked domain.  a and b control the length and height.  c
  % controls the width of the gap, and order controls the
  % regularity

  sa = o.diffProp(X0);
  t = o.arc(sa);
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*cos(t); y = b*r.*sin(t);
  ind = abs(x) < 10;
%  ind = abs(x) < 0;
  scalRat = 2*c/(1+c)*(0.5-0.5*cos(pi*x(ind(1:end/2))/10)).^10 + ...   
      (1-c)/(1+c);
  y(ind) = y(ind).*[scalRat;scalRat];

%  t = (0:N-1)'*2*pi/N;
%  x = 10*cos(t); y = 10*sin(t);
  X0 = [x;y];
  % redistribute points so that it is equispaced in arclength

elseif any(strcmp(options,'chokeLonger'))
  a = 100; b = 25/2; c = 0.67; order = 8;
  % parameters for the boundary
  t = (0:N-1)'*2*pi/N;
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*cos(t); y = b*r.*sin(t);
  ind = abs(x) < 35;
  scalRat = 2*c/(1+c)*(0.5-0.5*cos(pi*x(ind(1:end/2))/35)).^10 + ...   
      (1-c)/(1+c);
  y(ind) = y(ind).*[scalRat;scalRat];
  X0 = [x;y];
  % choked domain.  a and b control the length and height.  c
  % controls the width of the gap, and order controls the
  % regularity

  sa = o.diffProp(X0);
  t = o.arc(sa);
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*cos(t); y = b*r.*sin(t);
  ind = abs(x) < 35;
  scalRat = 2*c/(1+c)*(0.5-0.5*cos(pi*x(ind(1:end/2))/35)).^10 + ...   
      (1-c)/(1+c);
  y(ind) = y(ind).*[scalRat;scalRat];

  X0 = [x;y];
  % redistribute points so that it is equispaced in arclength

elseif any(strcmp(options,'chokeLongest'))
%  a = 100; b = 25/2; c = 0.67; order = 8;
  a = 200; b = 25/2; c = 0.67; order = 8;
  % parameters for the boundary
  len = 100;
  % length of the constriction
  t = (0:N-1)'*2*pi/N;
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*cos(t); y = b*r.*sin(t);

  ind = abs(x) < len;
  scalRat = 2*c/(1+c)*(0.5-0.5*cos(pi*x(ind(1:end/2))/len)).^10 + ...   
      (1-c)/(1+c);
  y(ind) = y(ind).*[scalRat;scalRat];
  X0 = [x;y];
  % choked domain.  a and b control the length and height.  c
  % controls the width of the gap, and order controls the
  % regularity

  sa = o.diffProp(X0);
  t = o.arc(sa);
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*cos(t); y = b*r.*sin(t);
  ind = abs(x) < len;
  scalRat = 2*c/(1+c)*(0.5-0.5*cos(pi*x(ind(1:end/2))/len)).^10 + ...   
      (1-c)/(1+c);
  y(ind) = y(ind).*[scalRat;scalRat];

  X0 = [x;y];
  % redistribute points so that it is equispaced in arclength
%  plot(x,y,'b-o')
%  axis equal;
%  pause

elseif any(strcmp(options,'chokeMulti'))
  width = 1;
  Length = 10;
  gap = 0.5;

  NN = 32;
  L1 = Length + width;
  L2 = pi*(gap/2+width);
  L3 = Length;
  L4 = pi*gap/2;
  L5 = Length;
  L6 = pi*(gap/2+width);
  L7 = Length + width;
  L8 = width;
  L9 = Length + width;
  L10 = pi*gap/2;
  L11 = Length;
  L12 = pi*(gap/2+width);
  L13 = Length;
  L14 = pi*gap/2;
  L15 = Length + width;
  L16 = width;

  N1 = round(NN*L1);
  N2 = round(NN*L2);
  N3 = round(NN*L3);
  N4 = round(NN*L4);
  N5 = round(NN*L5);
  N6 = round(NN*L6);
  N7 = round(NN*L7);
  N8 = 5*round(NN*L8);
  N9 = round(NN*L9);
  N10 = round(NN*L10);
  N11 = round(NN*L11);
  N12 = round(NN*L12);
  N13 = round(NN*L13);
  N14 = round(NN*L14);
  N15 = round(NN*L15);
  N16 = 5*round(NN*L16);

  x1 = linspace(-(gap/2+width),Length,N1)';
  y1 = zeros(N1,1);

  t = linspace(-pi/2,pi/2,N2)';
  x2 = x1(end) + (gap/2+width)*(cos(t) - cos(t(1)));
  x2 = x2(2:end);
  y2 = y1(end) + (gap/2+width)*(sin(t) - sin(t(1)));
  y2 = y2(2:end);

  x3 = linspace(Length,0,N3)';
  x3 = x3(2:end);
  y3 = y2(end)*ones(N3,1);
  y3 = y3(2:end);

  t = linspace(-pi/2,-3*pi/2,N4)';
  x4 = x3(end) + (gap/2)*(cos(t) - cos(t(1)));
  x4 = x4(2:end);
  y4 = y3(end) + (gap/2)*(sin(t) - sin(t(1)));
  y4 = y4(2:end);

  x5 = linspace(0,Length,N5)';
  x5 = x5(2:end);
  y5 = (2*width+2*gap)*ones(N5,1);
  y5 = y5(2:end);

  t = linspace(-pi/2,pi/2,N6)';
  x6 = x5(end) + (gap/2+width)*(cos(t) - cos(t(1)));
  x6 = x6(2:end);
  y6 = y5(end) + (gap/2+width)*(sin(t) - sin(t(1)));
  y6 = y6(2:end);

  x7 = linspace(Length,-(gap/2+width),N7)';
  x7 = x7(2:end);
  y7 = y6(end)*ones(N7,1);
  y7 = y7(2:end);

  x8 = -(gap/2+width)*ones(N8,1);
  x8 = x8(2:end-1);
  y8 = linspace(4*width+3*gap,3*width+3*gap,N8)';
  y8 = y8(2:end-1);

  x9 = linspace(-(gap/2+width),Length,N9)';
  y9 = (3*width+3*gap)*ones(N9,1);

  t = linspace(+pi/2,-pi/2,N10)';
  x10 = x9(end) + (gap/2)*(cos(t) - cos(t(1)));
  x10 = x10(2:end);
  y10 = y9(end) + (gap/2)*(sin(t) - sin(t(1)));
  y10 = y10(2:end);

  x11 = linspace(Length,0,N11)';
  x11 = x11(2:end);
  y11 = y10(end)*ones(N11,1);
  y11 = y11(2:end);

  t = linspace(pi/2,3*pi/2,N12)';
  x12 = x11(end) + (gap/2+width)*(cos(t) - cos(t(1)));
  x12 = x12(2:end);
  y12 = y11(end) + (gap/2+width)*(sin(t) - sin(t(1)));
  y12 = y12(2:end);

  x13 = linspace(0,Length,N13)';
  x13 = x13(2:end);
  y13 = y12(end)*ones(N13,1);
  y13 = y13(2:end);

  t = linspace(+pi/2,-pi/2,N14)';
  x14 = x13(end) + (gap/2)*(cos(t) - cos(t(1)));
  x14 = x14(2:end);
  y14 = y13(end) + (gap/2)*(sin(t) - sin(t(1)));
  y14 = y14(2:end);

  x15 = linspace(Length,-(gap/2+width),N15)';
  x15 = x15(2:end);
  y15 = width*ones(N15,1);
  y15 = y15(2:end);

  x16 = -(gap/2+width)*ones(N16,1);
  x16 = x16(2:end-1);
  y16 = linspace(width,0,N16)';
  y16 = y16(2:end-1);

  xx = [x1;x2;x3;x4;x5;x6;x7;x8;x9;x10;x11;x12;x13;x14;x15;x16];
  yy = [y1;y2;y3;y4;y5;y6;y7;y8;y9;y10;y11;y12;y13;y14;y15;y16];
  z = xx+1i*yy;
  z = interpft(z,N);
  X = [real(z);imag(z)];
%  clf
%  plot(X0(1:end/2),X0(end/2+1:end),'k--')
%  axis equal;
%  pause

elseif any(strcmp(options,'slit'))
  a = 10; b = 3; c = 0.7; order = 8;
  % parameters for the boundary
  t = (0:N-1)'*2*pi/N;
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*cos(t); y = b*r.*sin(t);
  thresh = 0.5;
  ind = abs(x) < thresh;
  scalRat = 2*c/(1+c)*...
      (0.5-0.5*cos(pi*x(ind(1:end/2))/thresh)).^10 + ...   
      (1-c)/(1+c);
  y(ind) = y(ind).*[scalRat;scalRat];
  X0 = [x;y];

  sa = o.diffProp(X0);
  t = o.arc(sa);
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*cos(t); y = b*r.*sin(t);
  ind = abs(x) < thresh;
  scalRat = 2*c/(1+c)*...
      (0.5-0.5*cos(pi*x(ind(1:end/2))/thresh)).^10 + ...   
      (1-c)/(1+c);
  y(ind) = y(ind).*[scalRat;scalRat];
  X0 = [x;y];
  % choked domain.  a and b control the length and height.  c
  % controls the width of the gap, and order controls the
  % regularity

elseif any(strcmp(options,'contracting'))
  w = 0.5; % width of the opening
  ell1 = 3.0; % length before contracting region
  ell2 = 12; % length (in x) of contracting region
  ell3 = 3.0; % length after contracting region
  angd = 14; % angle in degrees of contracting region
  ang = angd*pi/180;

  W = 2*ell2*tan(ang) + w;  % total width of the geometry

  % Set up 10 marker points
  z1 = 0 - 1i*(W/2);
  z2 = z1 + ell1;
  z3 = z2 + (ell2 + 1i*ell2*tan(ang));
  z4 = z3 - 1i*ell2*tan(ang);
  z5 = z4 + ell3;
  z6 = z5 + 1i*W;
  z7 = z6 - ell3;
  z8 = z7 - 1i*ell2*tan(ang);
  z9 = z8 - (ell2 - 1i*ell2*tan(ang));
  z10 = z9 - ell1;

  dz = 0.1; % length spacing between successive nodes

  % setup roughly equispaced points along the geometry
  len1 = abs(z1 - z2);
  zz1 = linspace(z1,z2,round(len1/dz)); zz1 = zz1(1:end-1);
  len2 = abs(z2 - z3);
  zz2 = linspace(z2,z3,round(len2/dz)); zz2 = zz2(1:end-1);
  len3 = abs(z3 - z4);
  zz3 = linspace(z3,z4,round(len3/dz)); zz3 = zz3(1:end-1);
  len4 = abs(z4 - z5);
  zz4 = linspace(z4,z5,round(len4/dz)); zz4 = zz4(1:end-1);
  len5 = abs(z5 - z6);
  zz5 = linspace(z5,z6,round(len5/dz)); zz5 = zz5(1:end-1);
  len6 = abs(z6 - z7);
  zz6 = linspace(z6,z7,round(len6/dz)); zz6 = zz6(1:end-1);
  len7 = abs(z7 - z8);
  zz7 = linspace(z7,z8,round(len7/dz)); zz7 = zz7(1:end-1);
  len8 = abs(z8 - z9);
  zz8 = linspace(z8,z9,round(len8/dz)); zz8 = zz8(1:end-1);
  len9 = abs(z9 - z10);
  zz9 = linspace(z9,z10,round(len9/dz)); zz9 = zz9(1:end-1);
  len10 = abs(z10 - z1);
  zz10 = linspace(z10,z1,round(len10/dz)); zz10 = zz10(1:end-1);

  z = [zz1,zz2,zz3,zz4,zz5,zz6,zz7,zz8,zz9,zz10];
  winsize = 10; % Gaussian window size
  zPad = [z(end-winsize+1:end),z,z(1:winsize)];
  % pad so that periodicity is incorporated into smoothing
  z = smoothdata(zPad,'gaussian',winsize);
  % smooth padded data
  z = z(winsize+1:end - winsize);
  % remove the padding

  z = interpft(z,N);
  % interpolate with FFT to the requested number of discretization
  % points
  
  X0 = [real(z)';imag(z)'];
  % store as (x,y) rather than x + 1i*y

elseif any(strcmp(options,'doublechoke'))
  a = 10; b = 3; c = 0.6; order = 8;
  shift = pi/2 + 0.1;
  % parameters for the boundary
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*cos(t); y = b*r.*sin(t);
  ind = abs(x-shift) < pi/2;
  y(ind) = y(ind).*(1-c*cos(2*(x(ind)-shift)))/(1+c);
  ind = abs(x+shift) < pi/2;
  y(ind) = y(ind).*(1-c*cos(2*(x(ind)+shift)))/(1+c);
  X0 = [x;y];
  % choked domain.  a and b control the length and height.  c
  % controls the width of the gap, and order controls the
  % regularity.  shift controls the distance between the two
  % regions where the domain is restricted

elseif any(strcmp(options,'couette'))
  x = [20*cos(t)+cen(1,1) 10*cos(-t)+cen(1,2)];
  y = [20*sin(t)+cen(2,1) 10*sin(-t)+cen(2,2)];
  X = o.setXY(x,y);
  % annular domain

elseif any(strcmp(options,'doubleCouette'))
  x = [20*cos(t)+cen(1,1) 5*cos(-t)+cen(1,2) 5*cos(-t)+cen(1,3)];
  y = [20*sin(t)+cen(2,1) 5*sin(-t)+cen(2,2) 5*sin(-t)+cen(2,3)];
  X = o.setXY(x,y);
  % annular domain of genus 2

elseif any(strcmp(options,'quadCouette'))
  x = [20*cos(t)+cen(1,1) 5*cos(-t)+cen(1,2) 5*cos(-t)+cen(1,3) ...
                          5*cos(-t)+cen(1,4) 5*cos(-t)+cen(1,5)];
  y = [20*sin(t)+cen(2,1) 5*sin(-t)+cen(2,2) 5*sin(-t)+cen(2,3) ...
                          5*sin(-t)+cen(2,4) 5*sin(-t)+cen(2,5)];
  X = o.setXY(x,y);
  % annular domain of genus 4

elseif any(strcmp(options,'doubleFlower'))
  r = 17 + 2*cos(7*t);
  x = [r.*cos(t)+cen(1,1) 5*cos(-t)+cen(1,2) 5*cos(-t)+cen(1,3)];
  y = [r.*sin(t)+cen(2,1) 5*sin(-t)+cen(2,2) 5*sin(-t)+cen(2,3)];
  X = o.setXY(x,y);
  % annular domain of genus 2

elseif any(strcmp(options,'cylinder'))
  x = [20*scale*cos(t)+cen(1,1)];
  y = [20*scale*sin(t)+cen(2,1)];
  X = o.setXY(x,y);
  % single cylinder

elseif any(strcmp(options,'tube'))
  a = 10; b = 3; order = 8;
  % parameters for the boundary
  Nsides = ceil(0.5*b/(2*a+2*b)*N);
  Ntop = (N-4*Nsides)/2;
  t1 = linspace(0,0.2*pi,Nsides+1); t1 = t1(1:end-1)';
  t2 = linspace(0.2*pi,pi-0.2*pi,Ntop+1); t2 = t2(1:end-1)';
  t3 = linspace(pi-0.2*pi,pi+0.2*pi,2*Nsides+1); t3 = t3(1:end-1)';
  t4 = linspace(pi+0.2*pi,2*pi-0.2*pi,Ntop+1); t4 = t4(1:end-1)';
  t5 = linspace(2*pi-0.2*pi,2*pi,Nsides+1); t5 = t5(1:end-1)';
  t = [t1;t2;t3;t4;t5];
  % Parameterize t so that geometry is closer to 
  % equispaced in arclength
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*cos(t); y = b*r.*sin(t);
  X0 = [x;y];
  % rounded off cylinder.  a and b control the length and height 
  % and order controls the regularity
elseif any(strcmp(options,'porous'))
  a = 10; b = 3; order = 8;
  % parameters for the boundary
  Nsides = ceil(0.5*b/(2*a+2*b)*N);
  Ntop = (N-4*Nsides)/2;
  t1 = linspace(0,0.2*pi,Nsides+1); t1 = t1(1:end-1)';
  t2 = linspace(0.2*pi,pi-0.2*pi,Ntop+1); t2 = t2(1:end-1)';
  t3 = linspace(pi-0.2*pi,pi+0.2*pi,2*Nsides+1); t3 = t3(1:end-1)';
  t4 = linspace(pi+0.2*pi,2*pi-0.2*pi,Ntop+1); t4 = t4(1:end-1)';
  t5 = linspace(2*pi-0.2*pi,2*pi,Nsides+1); t5 = t5(1:end-1)';
  t = [t1;t2;t3;t4;t5];
  % Parameterize t so that geometry is closer to 
  % equispaced in arclength
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*cos(t); y = b*r.*sin(t);
  X0 = [x;y];
  % rounded off cylinder.  a and b control the length and height 
  % and order controls the regularity

  t = (0:N-1)'*2*pi/N;
  xpores = 0.5*cos(-t);
  ypores = 0.5*sin(-t);
  Xpores = [[xpores;ypores] [xpores-4;ypores+0.5] ...
            [xpores-2;ypores-0.5] [xpores-3;ypores-2.0]];
  
  X = [X0 Xpores];
else
  if ~isempty(ra)
    X0 = o.ellipse(N,ra);
    % build a vesicle of reduced area ra with N points
  else
    % this shape has reduced area very close to 0.65
    X0 = o.ellipse(N,0.65);
  end
  
  [ra,area,length] = o.geomProp(X0);

  if any(strcmp(options,'volFrac'))
    X0 = 1.55*X0;  
    ind = find(strcmp(options,'volFrac'));
    volFrac = options{ind+1};
    Xwalls = options{ind+2};
    fmm = options{ind+3};
    optionstt = options{ind+4};
    prams = options{ind+5};
    X = oc.fillDomain(X0,volFrac,Xwalls,fmm,optionstt,prams);
    nv = size(X,2);
  else
    X0 = 2*X0/sqrt(area/pi);
    % make the generic shape have an area equal to pi
%    scale = 2*scale/sqrt(area/pi);
%    X0 = scale*X0;
  end

  if size(cen,2) ~= nv
    b = 3*max(X0(1:N));
    cen = [b*(0:nv-1);zeros(1,nv)];
  end
  % centers if there are multiple vesicles.
end 
% end of building reference vesicles.  Only need to rotate
% and shift them as desired


if isempty(X)
  % if X has not been defined, we only have X0 which is 
  % the reference vesicle which not needs to be
  % rotated and translated
  X = zeros(2*N,nv);

  for k=1:nv
    X(1:N,k) = scale*(cos(theta(k)) * X0(1:N) - ...
      sin(theta(k)) * X0(N+1:2*N)) + cen(1,k);
    X(N+1:2*N,k) = scale*(sin(theta(k)) * X0(1:N) +  ...
      cos(theta(k)) * X0(N+1:2*N)) + cen(2,k);
  end
  % Rotate vesicles as requested 
end


end % initConfig

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X0 = ellipse(o,N,ra)
% X0 = o.ellipse(N,ra) finds the ellipse (a*cos(theta),sin(theta)) so
% that the reduced area is ra.  Uses N points.  Parameter a is found 
% by using bisection method

t = (0:N-1)'*2*pi/N;
a = (1 - sqrt(1-ra^2))/ra;
% initial guess using approximation length = sqrt(2)*pi*sqrt(a^2+1)

X0 = [a*cos(t);sin(t)];

[raNew,~,~] = o.geomProp(X0);

cond = abs(raNew - ra)/ra < 1e-4;
maxiter = 10;
iter = 0;
while (~cond && iter < maxiter);
  iter = iter + 1;
  if (raNew > ra)
    a = 0.9 * a;
  else
    a = 1.05 * a;
  end
  % update the major axis
  X0 = [cos(t);a*sin(t)];
  % Compute new possible configuration
  [raNew,~,~] = o.geomProp(X0);
  % compute new reduced area
  cond = abs(raNew - ra) < 1e-2;
  % check if the residual is small enough
end
% iteration quits if reduced area is achieved within 1% or 
% maxiter iterations have been performed

end % ellipse

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = fillDomain(o,Xref,volFrac,Xwalls,fmm,options,prams)
% fillDomain(X0,volFrac,Xwalls) uses a Monte-Carlo method to fill 
% the domain enclosed by Xwalls with dialations and translations
% of Xref.  Tries to obtain volume fraction volFrac
[xwalls,ywalls] = o.getXY(Xwalls);
% solid walls
[x0,y0] = o.getXY(Xref);
% reference solution
N = numel(x0);
%figure(1); clf; hold on
%plot(xwalls,ywalls,'k')
nvbd = size(xwalls,2);
% number of solid walls
[~,area,~] = o.geomProp(Xwalls);
areaGeom = sum(area);
% Total area occupied by the physical domain
[~,areaVes,~] = o.geomProp(Xref);
% area of reference vesicle

nv = ceil(volFrac*areaGeom/areaVes);
% total number of vesicles to achieve volume fraction
% volFrac

walls = capsules(Xwalls,[],[],0,0,true);
% build object for walls
radx = 1/2*(max(xwalls(:,1)) - min(xwalls(:,1)));
rady = 1/2*(max(ywalls(:,1)) - min(ywalls(:,1)));
% need radx and rady to decide where to randomly sample the geometry.
% This way we do not pick a center that is  outside of the solid walls
% on a regular basis

X = zeros(2*N,nv);
k = 1;
% counter for the number of successfully placed vesicles
tt = tstep(options,prams);
while k <= nv 
  cx = 2*radx*(rand-1/2);
  cy = 2*rady*(rand-1/2);
  
  phi = 2*pi*rand;
  % center
  xpot = cx + x0*cos(phi) + y0*sin(phi);
  ypot = cy - x0*sin(phi) + y0*cos(phi);
  % potential vesicle location

  accept = true; % tentatively accept the vesicle

  vesicle = capsules([xpot;ypot],[],[],[],[],true);
  [~,NearV2W] = vesicle.getZone(walls,3);
  % create capsule with the potential vesicle
  [~,icollisionWall] = vesicle.collision(walls,[],NearV2W,fmm,1,tt.op);
  if icollisionWall
    accept = false;
    % at least one of the vesicles's points is outside of one
    % of the solid wall components
  end
  
  if 1.1*sqrt(mean(xpot)^2 + mean(ypot)^2)>=max(xwalls(:,1)) || ...
          0.9*sqrt(mean(xpot)^2 + mean(ypot)^2)<=max(xwalls(:,2)) 
    accept = false;
  end
  % reject vesicle if it is outside the domain
  
  if accept 
    X(:,k) = [xpot;ypot];
    % if vesicle is not outside of physical walls, accept it as a
    % potential new vesicle.  It will be kept as long as it intersects
    % no other vesicles

    vesicle = capsules(X(:,1:k),[],[],0,0,true);
    % create an object with the current configuration of vesicles
    NearV2V = vesicle.getZone([],1);
    icollisionVes = vesicle.collision(...
        [],NearV2V,[],fmm,1,tt.op);
    % see if vesicles have crossed using collision detection code
    % Can be used with or without the fmm
    if icollisionVes
      X(:,k) = 0;
      % if they've crossed, reject it
    else
      k = k + 1;
      fprintf('%d Vesicles left to fill the domain to the desired volume fraction\n', nv-k+1)
      % if they haven't crossed, increase the total number of vesicles
    end
  end

end % while k <= nv

end % fillDomain


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xnew = correctAreaAndLength(o,X,timeTolerance,om)
% Xnew = correctAreaAndLength(X,a0,l0) changes the shape of the vesicle
% by finding the shape Xnew that is closest to X in the L2 sense and
% has the same area and length as the original shape

% tolConstraint (which controls area and length) comes from the
% area-length tolerance for time adaptivity.

a0 = om.area;
l0 = om.length;
[~,at,lt] = o.geomProp(X);
eAt = abs(at-a0)./a0;
eLt = abs(lt-l0)./l0;

N  = size(X,1)/2;

tolConstraint = 1e-4; % 1 percent error in constraints
% tolConstraint = timeTolerance;
tolFunctional = 1e-4; % Allowed to change shape by 1 percent

options = optimset('Algorithm','sqp','TolCon',tolConstraint,...
    'TolFun',tolFunctional,'display','off','MaxFunEvals',3000);

% Options for Algorithm are:
% 'active-set', 'interior-point', 'interior-point-convex' 'sqp'

Xnew = zeros(size(X));

for k = 1:size(Xnew,2);
    
  minFun = @(z) 1/N*min(sum((z - X(:,k)).^2));
  [Xnew(:,k),~,iflag] = fmincon(minFun,X(:,k),[],[],[],[],[],[],...
      @(z) o.nonlcon(z,a0(k),l0(k)),options);
  if iflag~=1 && iflag~=2
    message = ['Correction scheme failed, do not correct at this step'];
    om.writeMessage(message,'%s\n')
    Xnew(:,k) = X(:,k);
  end
  % if fmincon fails, keep the current iterate for this time step.
  % Hopefully it'll be corrected at a later step.
  
end
% Looping over vesicles, correct the area and length of each vesicle

end % correctAreaAndLength

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cIn,cEx] = nonlcon(o,X,a0,l0)
% [cIn,cEx] = nonlcon(X,a0,l0) is the non-linear constraints required
% by fmincon

[~,a,l] = o.geomProp(X);

cIn = [];
% new inequalities in the constraint
cEx = [(a-a0)/a0 (l-l0)/l0];
% want to keep the area and length the same

end % nonlcon

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alpha = arc(o,sa)
% alpha = arc(sa) computes the parameter values of alpha that can be
% used to distribute points evenly in arclength

N = numel(sa);
length = sum(sa)*2*pi/N; 
% total length of vesicle. This is spectrally accurate since it is the
% trapezoid rule applied to a periodic function
tol = 1e-13; % tolerance

IK = fft1.modes(N,1);
intsa = fft1.intFT(sa,IK);

f = (0:N-1)'*length/N - intsa;
% we are interested in fiding the N roots of f. First one is at 0
df = -sa;
% derivative of f so that we can apply Newton's method
fh = fft(f)/N;

dfh = fft(df)/N;
% Fourier coeffients of f and its derivative

alpha = zeros(N,1);

for j = 2:N
  alpha(j) = alpha(j-1);
  for k = 1:100 % apply Newton's method
    falpha = real(sum(fh.*exp(IK*alpha(j)))) + ...
        ((j-1)*length/N - alpha(j)/2/pi*length);
    % f at alpha
    dfalpha = real(sum(dfh.*exp(IK*alpha(j))));
    % derivative of f at alpha
    alpha(j) = alpha(j) - falpha/dfalpha;
    % apply Newton iteration
    if abs(abs(falpha) < tol)
      break
    end
  end
end


end % arc



    
end % methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = arcDeriv(f,m,sa,IK)
% f = arcDeriv(f,m,s,IK,col) is the arclength derivative of order m.
% f is a matrix of scalar functions (each function is a column)
% f is assumed to have an arbitrary parametrization
% sa = d s/ d a, where a is the aribtrary parameterization
% IK is the fourier modes which is saved and used to accelerate 
% this routine

for j=1:m
  f = sa.*ifft(IK.*fft(f));
end
f = real(f);

end % arcDeriv


end % methods (Static)

end % classdef
