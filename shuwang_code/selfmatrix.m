function  A = selfmatrix(s,kmatrix)

N = numel(s.x); 
r = repmat(s.x,[1 N]) - repmat(s.x.',[N 1]) + diag(ones(1,N));
% target minus source for the single-layer potential, but putting ones
% on the main diagonal. x-coordinate is in the real part and
% y-coordinate is in the imaginary part

rt = repmat([1:N]'*pi/N , [1 N]) - ...
     repmat([1:N]*pi/N, [N 1]) + diag(ones(1,N));
% only going from 0 to pi to account for the divided by 2 in the
% argument of sin in equation (43). ie. rt is (alpha - alpha')/2 in
% equation (43) where alpha and alpha' go from 0 to 2*pi. Again, putting
% ones in the diagonal instead the natural limit which is 0

denom = abs(sin(rt));
% denominator of equation (43)

irr = 1./(conj(r).*r);    % 1/r^2
d1 = real(r); % x-coordinate of r
d2 = imag(r); % y-coordinate of r
Ilogr = -log(abs(r)./denom);  % log(1/r) diag block
% Putting ones in the diagonal of r and denom results in taking log of
% something non-zero. Not clear how (or if) they use the true limiting
% term of this smooth kernel which is log(abs(r))

A12 = d1.*d2.*irr;   % off diag vel block
A = [(Ilogr + d1.^2.*irr), A12;...
     A12, (Ilogr + d2.^2.*irr)]; 
% d1^2.*irr is the r1^2/r^2 in equation (6)
% d2^2.*irr is the r2^2/r^2 in equation (7)
% d1.*d2.*irr is the r1*r2/r^2 in equations (6) and (7)

A=A.*kmatrix;

end


