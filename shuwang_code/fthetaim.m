function fntheta = fthetaim(m,sl,theta,bendsti,un,utt)

ftheta = forctheta(m,sl,theta,un,utt);
% compute right hand side of equation (8) in Sohn et al JCP 2010

% now we must subtract off the stiffest part in Fourier space.
temp(1,1:m) = ftheta(1,1:m);
tempt(1,1:m) = theta(1,1:m)-pi*2*(0:m-1)/m;
temp(1,m+1) = temp(1,1);
tempt(1,m+1) = tempt(1,1);
temp = fft(temp,m);
tempt = fft(tempt,m);
% for nonlocal model, but everything is nonlocal in our formulation!

N = pi*2*[0 1:m/2 m/2-1:-1:1];
rlen = bendsti*(N/sl).^3/4;
% To the power of 3 term comes from the third-derivative of the function
% that L is being applied to

%temp2(1,1:m) = real(temp(1,1:m)) + rlen(1,1:m).*real(tempt(1,1:m));
%temp3(1,1:m) = imag(temp(1,1:m)) + rlen(1,1:m).*imag(tempt(1,1:m));
%temp4 = real(ifft(temp2 + 1i*temp3,m));

temp4 = real(ifft(temp(1,1:m) + rlen(1,1:m).*tempt(1,1:m)));

fntheta(1,1:m) = temp4(1,1:m);
%fntheta(1,m+1) = fntheta(1,1);
fntheta = kfilter(fntheta,m);
% Krasney filter applied to smooth out spurious high frequency terms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ftheta = forctheta(m,sl,theta,un,utt)

% compute first derivative of the opening angle using Fourier
% differentiation. This is actually the curvature
stheta = fd1x(theta,m);
% the `s' prefix corresponds to arclenght

% compute the derivative of the normal velocity
sun = fd1(un,m);

% Check T_s = -V * curvature (local inextensibility condition) (see
% equation (20) in Sohn et al JCP 2010)
ftheta(1,1:m) = (-sun(1,1:m) + stheta(1,1:m).*utt(1,1:m))/sl;
% see first half of equation (8) in Sohn et al JCP 2010.
ftheta(1,m+1) = ftheta(1,1);

end

end
