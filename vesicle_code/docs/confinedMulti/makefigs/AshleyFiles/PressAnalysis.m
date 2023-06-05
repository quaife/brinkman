meanpress = [];

for k = 1:509
%  figure(1); clf; hold on;
%  plot(posx,posy,'r');
%  plot(xx(10:end-9,k),yy(10:end-9,k),'k');
%  plot(xwalls,ywalls,'k');
%  axis([0.22 5.34 -1 1])
%  axis equal;
%%  disp(k)
%
%  figure(2); clf;
%  plot(yy(10:end-9,k),pp(10:end-9,k))
  meanpress = [meanpress; sum(pp(10:end-9,k))*1e-2];
%  pause(1)

end
