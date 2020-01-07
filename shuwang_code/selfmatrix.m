function  A = selfmatrix(s,kmatrix)

N = numel(s.x); 
M = numel(s.x);
r = repmat(s.x,[1 N]) - repmat(s.x.',[M 1]) + diag(ones(1,N));
% C-# displacements mat

rt = repmat([1:N]'*pi/N , [1 N]) - ...
     repmat([1:N]*pi/N, [M 1]) + diag(ones(1,N));

denom = abs(sin(rt));

irr = 1./(conj(r).*r);    % 1/r^2
d1 = real(r);
d2 = imag(r);
Ilogr = -log(abs(r)./denom);  % log(1/r) diag block

A12 = d1.*d2.*irr;   % off diag vel block
A = [(Ilogr + d1.^2.*irr), A12;...
     A12, (Ilogr + d2.^2.*irr)]; 

A=A.*kmatrix;

end


