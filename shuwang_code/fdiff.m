function x = fdiff(a,n)

y = a;
y(1,n+1) = y(1,1);
z = fft(y(1,1:n+1),n);
s = 2*pi*1i*[0:1:n/2 -n/2+1:1:-1];
z = s.*z;
x(1,1:n) = real(ifft(z,n));
x(1,n+1) = x(1,1);

end
