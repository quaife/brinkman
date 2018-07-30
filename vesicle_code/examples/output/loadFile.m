function [posx,posy,ten,shearStress,normalStress,...
    wallx,wally,ea,el,time,n,nv] = loadFile(file)
fid = fopen(file,'r');
val = fread(fid,'double');
fclose(fid);
n = val(1);
nv = val(2);
nbd = val(3);
nvbd = val(4);
walls = val(5:5+2*nbd*nvbd-1);
val = val(5+2*nbd*nvbd:end);

ntime = numel(val)/(5*n*nv+3);
if ntime ~= ceil(ntime);
  disp('PROBLEM WITH VALUES FOR n AND nv');
end

wallx = zeros(nbd,nvbd);
wally = zeros(nbd,nvbd);
for k = 1:nvbd
  istart = (k-1)*nbd+1;
  iend = k*nbd;
  wallx(:,k) = walls(istart:iend);
  wally(:,k) = walls((istart:iend)+nbd*nvbd);
end
posx = zeros(n,nv,ntime);
posy = zeros(n,nv,ntime);
ten = zeros(n,nv,ntime);
shearStress = zeros(n,nv,ntime);
normalStress = zeros(n,nv,ntime);
time = zeros(ntime,1);
ea = zeros(ntime,1);
el = zeros(ntime,1);

istart = 1;
for m = 1:ntime
  for k=1:nv
    iend = istart + n - 1;
    posx(:,k,m) = val(istart:iend);
    istart = iend + 1;
  end
  % load x positions

  for k=1:nv
    iend = istart + n - 1;
    posy(:,k,m) = val(istart:iend);
    istart = iend + 1;
  end
  % load y positions

  for k=1:nv
    iend = istart + n - 1;
    ten(:,k,m) = val(istart:iend);
    istart = iend + 1;
  end
  % load tensions

  for k = 1:nv
    iend = istart + n - 1;
    shearStress(:,k,m) = val(istart:iend);
    istart = iend + 1;
  end
  % shear stress

  for k = 1:nv
    iend = istart + n - 1;
    normalStress(:,k,m) = val(istart:iend);
    istart = iend + 1;
  end
  % normal stress

  ea(m) = val(istart);
  el(m) = val(istart+1);
  time(m) = val(istart+2);
  istart = istart + 3;

end



