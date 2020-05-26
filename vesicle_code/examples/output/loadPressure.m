function [pressx,pressy,press] = loadFile(file)
fid = fopen(file,'r');
val = fread(fid,'double');
fclose(fid);
n = val(1);
val = val(2:end);

ntime = (numel(val) - 2*n)/n;
% 2 positions
if ntime ~= ceil(ntime);
  disp('PROBLEM WITH VALUES FOR n AND nv');
end

press = zeros(n,ntime);

istart = 1;
iend = istart + n - 1;
pressx = val(istart:iend);
% load x positions
istart = iend + 1;
iend = istart + n - 1;
pressy = val(istart:iend);
% load y positions
istart = iend + 1;

for m = 1:ntime
  iend = istart + n - 1;
  press(:,m) = val(istart:iend);
  istart = iend + 1;
end
% load pressure


