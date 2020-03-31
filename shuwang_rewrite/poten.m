classdef poten


properties
N

end % properties


methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = poten(N);
o.N = N;

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

end % methods

end % classdef

