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
function logKernel = IntegrateLogKernel_old(o,sigma)
% logKernel = IntegrateLogKernel evaluates the integral of log|sin|: L

% See first term in equation (43)
bsigma = sum(sigma)/o.N;
%compute the fft of sigma
ssigma = fft(sigma,o.N); 
%define the fourier coefficients
coeff = -[1 1:o.N/2 o.N/2-1:-1:1]'*2*pi;
csigma = ssigma./coeff;
%zero out the zero mode
csigma(1,1) = 0;
ssigma = real(ifft(csigma,o.N));
%Construct the log kernel 
logKernel = ssigma/2+bsigma*log(1/2)/2/pi;

end %IntegrateLogKernel_old

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % methods

end % classdef

