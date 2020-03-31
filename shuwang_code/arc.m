function slo = arc(ds,n)

slo = zeros(1,n+1);

% integrate ds and set-up functions for interpolation and iteration
N = n;
tol = 1e-13;
dsi = fin(ds,n);
% get the integration of ds;
 
ave = dsi(N+1);
% ave is the length of the vesicle
f = (0:1:N)*ave/N - dsi;
fp = -ds;
a = f/N;
c = fp/N;
x = fft(a(1:N+1),N);
y = fft(c(1:N+1),N);
x(n+1) = x(1);
y(n+1) = y(1);

%%%%% compute the interpolation
slo(1) = 0;
for j = 2:N+1
  ao=slo(j-1)+(ave/ds(j))/N;
  for k=1:100
    aoo=2*pi*ao;
    v = synth1(aoo,N,x);
    vp = synth1(aoo,N,y);
    v = v + ((j-1)/N-ao)*ave;
    dalp = v/vp;
    an = ao - dalp;
    err = abs(dalp);
    if err <= tol
      slo(j) = an;
      break
    end
    ao=an;
  end
end


end
