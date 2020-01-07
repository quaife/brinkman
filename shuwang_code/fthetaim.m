function fntheta = fthetaim(m,sl,theta,bendsti,un,utt)

ftheta = forctheta(m,sl,theta,un,utt);

%  now we must subtract off the stiffest part in Fourier space.
temp(1,1:m) = ftheta(1,1:m);
tempt(1,1:m) = theta(1,1:m)-pi*2*(0:m-1)/m;
temp(1,m+1) = temp(1,1);
tempt(1,m+1) = tempt(1,1);
temp = fft(temp,m);
tempt = fft(tempt,m);
% for nonlocal model

N = pi*2*[0 1:m/2 m/2-1:-1:1];
rlen = bendsti*(N/sl).^3/4;

temp2(1,1:m) = real(temp(1,1:m))+rlen(1,1:m).*real(tempt(1,1:m));
temp3(1,1:m) = imag(temp(1,1:m))+rlen(1,1:m).*imag(tempt(1,1:m));

temp4 = real(ifft(temp2+sqrt(-1)*temp3,m));

fntheta(1,1:m) = temp4(1,1:m);
%fntheta(1,m+1) = fntheta(1,1);
fntheta = kfilter(fntheta,m);

function ftheta = forctheta(m,sl,theta,un,utt)
stheta = fd1x(theta,m);

sun = fd1(un,m);
%c to Check T_s = -VH (local Laglange condition)
ftheta(1,1:m) = (-sun(1,1:m)+stheta(1,1:m).*utt(1,1:m))/sl;
ftheta(1,m+1) = ftheta(1,1);

end

end
