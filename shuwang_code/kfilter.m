% fourier filter with krasny filtering with 1 input vector USING FAST.f
% FOURIER TRANSFORM

function u = kfilter(u,m)

%  figure(1); clf;
%  plot(u); hold on;
b(1,1:m) = u(1,1:m);
b(1,m+1) = b(1,1);
% filter fourier coeffs
b = fft(b,m);
dcut = 1e-10;  
for j=1:m
  testt = abs(b(1,j)/m);
  if testt < dcut
    b(1,j)=0;
  end
end 

for k = m/2:m/2+2
  b(1,k) = 0;
end 

c(1,1:m) = real(b(1,1:m)).*exp(-10*([0:2:m m-2:-2:2]/m).^25);
d(1,1:m) = imag(b(1,1:m)).*exp(-10*([1:2:m-1 0 m-1:-2:3]/m).^25);
d(1,m/2-1) = 0;
d(1,m/2+3) = 0;
b = real(ifft(c+1i*d,m));
u(1,1:m) = b(1,1:m);
%u(1,m+1) = b(1,1);
%plot(u,'r--')
%pause

end
