function [posx,posy,conc,ea,el,time,xvel,yvel,ten] = loadFile(file)
fid = fopen(file,'r');
val = fread(fid,'double');
fclose(fid);
n = val(1);
nv = 1;
%nv = val(2);
%nbd = val(3);
%nvbd = val(4);
%walls = val(5:5+2*nbd*nvbd-1);
%if numel(val) > 1e7
%  val = val(end-100*(3*n+3)+1:end);
%  % If file is huge, only look at last 20 time steps or will run out of
%  % RAM
%else
  %val = val(2+2*nbd*nvbd:end);
  val = val(2:end);
%end
ntime = numel(val)/(6*n*nv+3);

% 2 positions, concentration, two stresses,tension
% error in area, error in length
if ntime ~= ceil(ntime);
  disp('PROBLEM WITH VALUES FOR n AND nv');
end
ntime = floor(ntime);

% wallx = zeros(nbd,nvbd);
% wally = zeros(nbd,nvbd);
% for k = 1:nvbd
%   istart = (k-1)*nbd+1;
%   iend = k*nbd;
%   wallx(:,k) = walls(istart:iend);
%   wally(:,k) = walls((istart:iend)+nbd*nvbd);
% end
posx = zeros(n,nv,ntime);
posy = zeros(n,nv,ntime);
conc = zeros(n,nv,ntime);
ten  = zeros(n,nv,ntime);
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
    conc(:,k,m) = val(istart:iend);
    istart = iend + 1;
  end
  % load tensions
  ea(m) = val(istart);
  el(m) = val(istart+1);
  time(m) = val(istart+2);
  
  istart = istart + 3;
  for k=1:nv
    iend = istart + n - 1;
    xvel(:,k,m) = val(istart:iend);
    istart = iend + 1;
  end
  
  for k=1:nv
    iend = istart + n - 1;
    yvel(:,k,m) = val(istart:iend);
    istart = iend + 1;
  end
  
  for k=1:nv
    iend = istart + n - 1;
    ten(:,k,m) = val(istart:iend);
    istart = iend + 1;
  end
  
end



