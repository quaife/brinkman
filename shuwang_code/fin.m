% this code forms the indefinite integral of a periodic function
function yi=fin(y,ngrid)

n = ngrid;
N = (0:1:n);
a = [y y(1)];
b = fft(a,n);
aa = real(b(1));
b(1) = 0;
c = 1i*[1 1:1:n/2  -n/2+1:1:-1]*2*pi;
b = b./c;
%plot(b)
%pause
yi = real(ifft(b,n));
yi(n+1) = yi(1);
yi = aa/n*N/n+yi - yi(1);
% disp('Here')
% figure(1); clf
% plot(yi)
% figure(2); clf;
% semilogy(fftshift(abs(b)))
% pause
end
