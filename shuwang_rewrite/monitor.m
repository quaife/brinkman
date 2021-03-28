classdef monitor
% Used for doing input/output of data.
% Can do plotting, write to files, compute error in area
% and length, and write to console

properties
area            % area of initial vesicles 
length          % length of initial vesicles
reducedArea     % reduced area
verbose         % write data to console
saveData        % save the vesicle positions, conc, etc
usePlot         % plot the vesicles
dataFile        % name of the file to write the data
logFile         % name of the file to write the log
T               % time horizion
m               % number of time steps
N               % points per vesicle
cls             % concentration

end % properties

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = monitor(X,params,options)
% monitor(X,Xwalls,options,prams) saves options and parameters needed
% by the class and also computes the initial error in area and length
% so that the errors in area and length can be monitored throughout the
% simulation.
% This is the constructor
o.N = params.N; % points per vesicle
o.T = params.T; % time horizon
o.m = params.T/params.dt; % number of time steps
o.cls = params.concentra; % concentration of lipid species
o.verbose = options.verbose; % write data to console
o.saveData = options.saveData; % save the data
o.usePlot = options.usePlot; % plot the data
o.dataFile = options.dataFile; % data file name
o.logFile = options.logFile; % log file name
oc = curve; %shorthand for curve class
% compute the area, length, and reduced area of initial shape
[o.reducedArea,o.area,o.length] = oc.geomProp(X);

end % constructor: monitor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initializeFiles(o,X,conc,time,vel,ten)
% initializeFiles(X,sig,eta,RS,Xwalls) does the initial writing of data
% to files and the console.  It first deletes any previous data and then
% writes the number of points, tracer initial conditions, pressure
% targets X and Xwalls are the vesicles and solid walls.

N = o.N; % points per vesicle
oc = curve; % shorthand for curve class
% Save data to files if desired
if o.saveData
  fid = fopen(o.dataFile,'w');
  %write number of points and vesicles to data file
  fwrite(fid,N,'double');
  fclose(fid);
  fid = fopen(o.logFile,'w');
  fclose(fid);
end
o.writeMessage(' ','%s\n')
message = ['************* PHYSICAL PARAMETERS *************'];
o.writeMessage(message,'%s\n')
% write number of points
message = [num2str(N) ' points per vesicle'];
o.writeMessage(message,'%s\n')
%write the timestep size
message = ['Time step size is ' ...
  num2str(o.T/o.m)];
o.writeMessage(message,'%s\n')
% write time step size or that we are doing adaptive time stepping
message = ['The concentration of lipid species is: ' ...
    num2str(o.cls,'%4.2e')];
o.writeMessage(message,'%s\n')
% compute the errors in length and area  
[ea,el] = o.errors(X);
% Save data
if o.saveData
  o.writeData(X,conc,ea,el,time,vel,ten, N);
  % save initial configuartion
  message = ['Initial Area is:                ' ...
      num2str(o.area(1),'%10.2e')];
  o.writeMessage(message,'%s\n')
  message = ['Initial Length is:              ' ...
      num2str(o.length(1),'%10.2e')];
  o.writeMessage(message,'%s\n')
  message = ['Initial Reduced Area is:        ' ...
      num2str(o.reducedArea(1),'%10.2e')];
  o.writeMessage(message,'%s\n\n')
end
% write initial reduced area, area, and length to log file
o.writeStars
o.writeMessage(' ','%s\n')

end % initializeFiles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function terminate = outputInfo(o,X,conc,time,vel,ten)
% computes the error in area and length and write messages to the data
% file, the log file, and the console.

%Compute the error in area and length 
[ea,el] = o.errors(X);
if(ea > 1)
    terminate = 1;%'True';
else
    terminate = 0;%'False';
end

% Begin plotting
if o.usePlot
  o.plotData(X,time,ea,el,vel,conc);
  pause(0.01)
end
% End plotting

% Begin saving data
if o.saveData
  % don't want to save initial small time steps, but still want to check
  % the error in area and length so that the simulation is killed early
  % on if need be
  o.writeData(X,conc,ea,el,time,vel,ten, o.N);
end
% End saving data

% Begin sending messages to log file and console

message = ['t = ' num2str(time,'%4.2e') ...
      ' of T = ' num2str(o.T,'%4.2e')]; 
o.writeMessage(message,'%s\n')
message1 = ['Max error in area is   ' num2str(ea,'%4.2e')];
message2 = ['Max error in length is ' num2str(el,'%4.2e')];
o.writeMessage(message1,'%s\n')
o.writeMessage(message2,'%s\n')
o.writeMessage(' ','%s\n')
% End sending data to files and console

end % outputInfo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ea,el] = errors(o,X)
% function [ea,el] = errors(X) computes the errors in area and length
% of the new vesicle position

oc = curve;
[~,a,l] = oc.geomProp(X);
% compute new areas and length

ea = max(abs(a./o.area - 1));
el = max(abs(l./o.length - 1));
% compute error in area and length

end % errors

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeStars(o)
% writeStars writes a message of stars to the console and the log file
% depending on verbose and saveData

messageStars = '*********************************************';
o.writeMessage(messageStars,'%s\n')

end % writeStars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeMessage(o,message,format)
% function writeMessage(message,format) appends message to o.fileName
% with format

if o.saveData
  fid = fopen(o.logFile,'a');
  fprintf(fid,format,message);
  fclose(fid);
end
% save to log file
if o.verbose
  disp(message)
end
% write to console

end % writeMessage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotData(o,X,time,ea,el,vel,conc)
% plotData(X,Xwalls,time,ea,el) plots the current configuration with
% title X is the vesicle position time is the current time, ea and el
% are the errors in area and length
%TODO: Edit this so that ea and el are computed instead of inputs

oc = curve; % shorthand for curve class
 [x,y] = oc.getXY(X); % seperate x and y coordinates
% figure(1); clf; % hold on

%----------- General Position Plot ------------------------------------
%plot([x;x(1,:)],[y;y(1,:)],'r','linewidth',2)

%----------- Velocity Quiver Plot -------------------------------------
%[xvel,yvel] = oc.getXY(vel);
%quiver(x,y,xvel,yvel)

%----------- Cline Bending Stiffness Plot ---------------------------------
clf;
rbn = 1 * (ones(o.N,1) - conc) + 0.1*conc;
h = cline([x;x(1,:)],[y;y(1,:)],[rbn;rbn(1,:)]);
set(h,'linewidth',3)
colorbar
%Title and axis settings for all plots
titleStr = ['t = ' num2str(time,'%4.2e') ...
  ' eA = ' num2str(ea,'%4.2e') ...
  ' eL = ' num2str(el,'%4.2e')];
title(titleStr)
axis equal

end % plotData

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeData(o,X,conc,ea,el,time,vel,ten,N)
% writeData(X,sigma,ea,el,time,res) writes the position, concentration,
% errors, and time to a binary file.  Matlab can later read this file to
% postprocess the data
 
oc = curve;
[x,y] = oc.getXY(X);
[xvel,yvel] = oc.getXY(vel);
output = [x(:);y(:);conc(:);ea;el;time;xvel(:);yvel(:);ten(:)];
% format that postProcess/loadfile.m reads the output
fid = fopen(o.dataFile,'a');
fwrite(fid,output,'double');
fclose(fid);

end % writeData


end % methods


end % classdef

