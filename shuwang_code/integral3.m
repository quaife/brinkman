% integral of log|sin|: L
function vel = integral3(sigma,m)

% make sure velocity is zero initially
vel = zeros(1,m+1);
bsigma = sum(sigma(1,1:m))/m;
% compute HI(ssigma)
ssigma = fft(sigma,m);  
coeff = -[1 1:m/2 m/2-1:-1:1]*2*pi;
csigma = ssigma./coeff;
csigma(1,1) = 0;
ssigma = real(ifft(csigma,m));


%now ssigma=HD(ssigma) and now construct velocity
vel(1,1:m) = ssigma(1,1:m)/2+bsigma*log(1/2)/2/pi;
vel(1,m+1) = vel(1,1);

end

