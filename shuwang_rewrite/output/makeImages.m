addpath ..
set(0,'DefaultAxesFontSize',22)
options.savefig = false;

irate = 1; % controls the speed of the visualization

file = 'relaxation1VesData.bin';
[posx,posy,curv,ea,el,time,xvel,yvel] = loadFile(file);
% load positions, curvature, errors, time, and velocities
pause
istart = 1;
iend = numel(time);
ntime = iend;
ntime = numel(time);
%time = time(istart:iend);
