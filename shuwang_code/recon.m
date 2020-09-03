% given an arclength sl, and a tan angle theta, this code reconstructs
% the curve and gives x,y need equal arclength to do
function [x,y] = recon(m,x0,y0,sl,theta)

temp(1,1:m) = sl*cos(theta(1,1:m));
temp1(1,1:m) = sl*sin(theta(1,1:m));       
y = fin(temp1(1,1:m),m);
x = fin(temp(1,1:m),m);
avx = x(1,m+1);
avy = y(1,m+1);

x(1,2:m+1) = x(1,2:m+1) - (1:m)*avx/m + x0;
y(1,2:m+1) = y(1,2:m+1) - (1:m)*avy/m + y0;

x(1,1) = x0;
y(1,1) = y0;
%if numel(theta) == 96
%  clf
%  semilogy(abs(fftshift(fft(theta - (0:95)*2*pi/96))))
%  hold on
%  semilogy(abs(fftshift(fft(x(1:end-1)))),'r')
%  pause
%end
