classdef poten
%This class defines single layers for stokes kernels on 2D periodic curves.
%It also defines the matricies that map a density function defined on the
%boundary of a curve to the layer potential evaluated on the curve. This
%class also has the main routine that evaluates layer-potentials using 
%near-singular integration.

properties
N; % points per curve
end % properties

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function po = poten(N)
% o = poten(N) is a constructor the initializes the class
po.N = N;
end % constructor: poten

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function kmatrix = oddEvenMatrix(o)
% kmatrix = oddEvenMatrix builds a matrix that is 0 if the row and
% column have the same parity (both even or odd) and is 1 otherwise.
% This matrix is used to do odd-even integration with kernels that have
% a weak singularity at the diagonal terms such as the Stokes
% single-layer potential

kmatrix = ones(2,2*o.N);
kmatrix(1,1:2:end) = 0;
kmatrix(2,2:2:end) = 0;
kmatrix = repmat(kmatrix,o.N,1);

end % oddEvenMatrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function selfmatrix = StokesMatrixLogless(o,ves)
% selfmatrix = StokesMatrixLogless builds the Stokes matrix without the
% log singularity. Contains the rightmost kernel in equation (43).

X = ves.X;
% define points
pts = X(1:o.N) + 1i*X(o.N+1:end);
% r is the target minus source for the single-layer potential with ones
% on the main diagonal. The x-coordinate is in the real part and the y-
% coordinate is in the imaginary part
r = repmat(pts,[1 o.N]) - repmat(transpose(pts),[o.N 1]) + 1*eye(o.N);
% rt is (alpha - alpha')/2 in equation (43) where alpha and alpha' go
% from 0 to 2*pi. We are only going from 0 to pi to account for the
% divided by 2 in the argument of sin in equation (43). We put ones on
% the diagonal, the limit which is something smooth (see equation (44))
rt = repmat((1:o.N)'*pi/o.N , [1 o.N]) - ...
     repmat((1:o.N)*pi/o.N, [o.N 1]) + 1*eye(o.N);
% Compute the denominator of equation (43). Putting ones in the diagonal
% of r and denom results in taking log of something non-zero.
denom43 = 2*abs(sin(rt));
irr = 1./(conj(r).*r);    % 1/r^2
d1 = real(r); % x-coordinate of r
d2 = imag(r); % y-coordinate of r
Ilogr = log(abs(r)./denom43);  % log(1/r) diag block
salpha = ves.L/(2*pi); % limiting value for log part of kernel
Ilogr(1:o.N+1:end) = log(salpha); % correct the diagonal terms
oc = curve(o.N); %shorthand for curve class
% Compute the derivatives of x and y
[dx,dy] = oc.getDXY(ves.X);
% Divide by length to get the arclength derivative
dx = dx/ves.L; dy = dy/ves.L;
A11 = d1.^2.*irr; % (1,1) block ie. rx^2/(rx^2 + ry^2)
A12 = d1.*d2.*irr; % (1,2) block ie. rx*ry/(rx^2 + ry^2)
A22 = d2.^2.*irr; % (2,2) block ie. ry^2/(rx^2 + ry^2)
% Limiting values of the above 3 blocks
A11(1:o.N+1:end) = dx.^2;
A12(1:o.N+1:end) = dx.*dy;
A22(1:o.N+1:end) = dy.^2;
% Compute the selfmatrix
selfmatrix = [-Ilogr + A11 A12; A12 -Ilogr + A22];

end % StokesMatrixLogless

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Symm_sigma = IntegrateLogKernel(o,sigma)
%IntegrateLogKernel returns the symms operator applied to sigma, where
%sigma is a density function.

%take the fft of sigma
sigmah = fft(sigma); 
%define the Fourier coefficients
if length(sigma) == 2*o.N
    coeff = [[(0:o.N/2-1)';(-o.N/2:-1)'];[(0:o.N/2-1)';(-o.N/2:-1)']];
else
    coeff = [(0:o.N/2-1)';(-o.N/2:-1)'];
end
Symm_sigmah = -sigmah./abs(coeff);
Symm_sigmah(1) = 0;
Symm_sigma = real(ifft(Symm_sigmah));

end %IntegrateLogKernel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = StokesDLP(o,geom)

oc = curve(geom.N);
% Geometry positions
[x,y] = oc.getXY(geom.X);

% number of points on geometry
N = geom.N;

% tangent vector
tx = cos(geom.theta);
ty = sin(geom.theta);

% curvature
cur = geom.cur';

% target points
xtar = x(:,ones(N,1))';
ytar = y(:,ones(N,1))';

% source points
xsou = x(:,ones(N,1));
ysou = y(:,ones(N,1));

% tangent at sources
txsou = tx';
tysou = ty';

% rx and ry terms
diffx = xtar - xsou;
diffy = ytar - ysou;
% 1 over the distance to the power of 4
rho4 = (diffx.^2 + diffy.^2).^(-2);
% set diagonal terms to 0. Will fix with limiting value involving the
% curvature a bit later
rho4(1:N+1:N.^2) = 0;

kernel = diffx.*(tysou(ones(N,1),:)) - ...
        diffy.*(txsou(ones(N,1),:));
kernel = kernel.*rho4;

% (1,1) component
D11 = kernel.*diffx.^2;
% diagonal limiting term
D11(1:N+1:N.^2) = 0.5*cur.*txsou.^2;

% (1,2) component
D12 = kernel.*diffx.*diffy;
% diagonal limiting term
D12(1:N+1:N.^2) = 0.5*cur.*txsou.*tysou;

% (2,2) component
D22 = kernel.*diffy.^2;
% diagonal limiting term
D22(1:N+1:N.^2) = 0.5*cur.*tysou.^2;

% build matrix with four blocks
D = [D11 D12; D12 D22];
% scale with the arclength spacing and divide by pi
D = 1/pi*D*geom.L/N;


end % StokesDLP


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START OF ROUTINES THAT EVALUATE STOKES LAYER POTENTIALS AT TARGET
% POINTS THAT DIFFER FROM THE SOURCCE POINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vel = StokesSLPtar(o,ves,trac,Xtar)

oc = curve(ves.N);

[xsou,ysou] = oc.getXY(ves.X);
[fx,fy] = oc.getXY(trac);

[xtar,ytar] = oc.getXY(Xtar);
xtar = xtar(:); ytar = ytar(:);

velx = zeros(size(xtar));
vely = zeros(size(xtar));


for k = 1:numel(xtar)
  rx = xtar(k) - xsou;
  ry = ytar(k) - ysou;
  rho2 = rx.^2 + ry.^2;
  rdotf = rx.*fx + ry.*fy;

  velx(k) = sum(-0.5*log(rho2).*fx + rdotf./rho2.*rx)/(4*pi);
  vely(k) = sum(-0.5*log(rho2).*fy + rdotf./rho2.*ry)/(4*pi);
end

velx = velx * ves.L/ves.N;
vely = vely * ves.L/ves.N;
vel = [velx;vely];

end % StokesSLPtar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vel = StokesDLPtar(o,geom,eta,Xtar)
% evaluate the double-layer potential due to the geometry stored in geom
% using the density function eta, and evaluating at the target points in
% Xtar

vel = zeros(size(Xtar));

end % StokesDLPtar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF ROUTINES THAT EVALUATE STOKES LAYER POTENTIALS AT TARGET POINTS
% THAT DIFFER FROM THE SOURCCE POINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end % methods

end % classdef

