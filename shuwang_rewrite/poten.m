classdef poten


properties
N

end % properties


methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function po = poten(N);
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
function selfmatrix = StokesMatrixLogless(o,X)
%selfmatrix = StokesMatrix builds the Stokes matrix without the log 
%singularity. Contains the rightmost kernel in equation (43).
pts = [X(1:o.N) + 1i*X(o.N+1:end)];
r = repmat(pts,[1 o.N]) - repmat(transpose(pts),[o.N 1]) + diag(ones(1,o.N));
% target minus source for the single-layer potential, but putting ones
% on the main diagonal. x-coordinate is in the real part and
% y-coordinate is in the imaginary part
rt = repmat([1:o.N]'*pi/o.N , [1 o.N]) - ...
     repmat([1:o.N]*pi/o.N, [o.N 1]) + diag(ones(1,o.N));
% only going from 0 to pi to account for the divided by 2 in the
% argument of sin in equation (43). ie. rt is (alpha - alpha')/2 in
% equation (43) where alpha and alpha' go from 0 to 2*pi. Again, putting
% ones in the diagonal instead the natural limit which is 0

denom43 = abs(sin(rt));
% denominator of equation (43)

irr = 1./(conj(r).*r);    % 1/r^2
d1 = real(r); % x-coordinate of rxxx
d2 = imag(r); % y-coordinate of r
Ilogr = -log(abs(r)./denom43);  % log(1/r) diag block
% Putting ones in the diagonal of r and denom results in taking log of
% something non-zero. Not clear how (or if) they use the true limiting
% term of this smooth kernel which is log(abs(r))

A12 = d1.*d2.*irr;   % off diag vel block
selfmatrix = [(Ilogr + d1.^2.*irr), A12;...
     A12, (Ilogr + d2.^2.*irr)]; 
% d1^2.*irr is the r1^2/r^2 in equation (6)
% d2^2.*irr is the r2^2/r^2 in equation (7)
% d1.*d2.*irr is the r1*r2/r^2 in equations (6) and (7)
kmatrix = o.oddEvenMatrix;
selfmatrix = 2*selfmatrix.*kmatrix;
% Expect the factor of 2 is because we are only using every other grid
% point, so the arclength spacing has doubled.

end %StokesMatrixLogless

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function logKernel = IntegrateLogKernel(o,sigma)
%logKernel = IntegrateLogKernel evaluates the integral of log|sin|: L
%See first term in equation (43)
bsigma = sum(sigma)/o.N;
%compute HI(ssigma)
ssigma = fft(sigma,o.N); 
%We may need to change these coeffs, I'm not sure.
coeff = -[1 1:o.N/2 o.N/2-1:-1:1]'*2*pi;
csigma = ssigma./coeff;
csigma(1,1) = 0;
ssigma = real(ifft(csigma,o.N));
%now ssigma=HD(ssigma) and we can construct the log kernel 
logKernel = ssigma/2+bsigma*log(1/2)/2/pi;

end %IntegrateLogKernel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function lambdaTilde = Stokes(
%end %Stokes




end % methods

end % classdef

