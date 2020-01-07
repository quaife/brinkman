function val = inclinationAngle(X)

n = length(X)/2;
c0 = fftshift(1i*(-n/2:1:n/2-1)'); D1 = @(f) real(ifft(c0.*fft(f)));
X1 = X(1:n); X2 = X(n+1:2*n);
sa = sqrt(D1(X1).^2 + D1(X2).^2);
nor = [D1(X2)./sa, -D1(X1)./sa];

c = [mean(X1), mean(X2)];
r = [X1-c(1), X2-c(2)];
rho2 = r(:,1).^2 + r(:,2).^2;
rn = sum(r.*nor, 2);

intg = @(f) 1/4*2*pi/n*sum(f.*sa.*rn);
J(1,1) = intg(rho2 - r(:,1).^2);
J(2,2) = intg(rho2 - r(:,2).^2);
J(1,2) = intg(-r(:,1).*r(:,2)); 
J(2,1) = J(1,2);

[V,D] = eig(J);
[~, ind] = min(D);

v = V(:,ind(1));
val = -atan2(v(1), v(2));

