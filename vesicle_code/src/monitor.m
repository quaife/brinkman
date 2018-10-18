classdef monitor
% Used for doing input/output of data.
% Can do plotting, write to files, compute error in area
% and length, and write to console

properties
area            % area of initial vesicles 
length          % length of initial vesicles
reducedArea     % reduced area
verbose         % write data to console
saveData        % save the vesicle positions, tension, etc
usePlot         % plot the vesicles
track           % place tracker points on the vesicles
quiver          % plot the velocity field on the vesicles
axis            % axis of the plot
dataFile        % name of the file to write the data
logFile         % name of the file to write the log
T               % time horizion
m               % number of time steps
errorTol        % error tolerance for errors in area and length
tracers         % use tracers to monitor passive particles
timeAdap        % use time adaptivity
N               % points per vesicle
nv              % number of vesicles
Nbd             % points per solid wall (if any)
nvbd            % number of solid wall components (if any)
order           % time stepping order
orderGL         % order of Gauss-Lobatto quadrature
pressure        % calculate pressure
rtolArea        % allowable relative area change
rtolLength      % allowable relative length change
adhesion        % flag to tell if adhesion is in model
adStrength      % strength of adhesion (W_0)
adRange         % range of adhesion (d_0)

end % properties

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = monitor(X,Xwalls,options,prams)
% monitor(X,Xwalls,options,prams) saves options and parameters needed
% by the class and also computes the initial error in area and length
% so that the errors in area and length can be monitored throughout the
% simulation.
% This is the constructor

o.N = prams.N;                      % points per vesicle
o.nv = prams.nv;                    % number of vesicles
o.Nbd = prams.Nbd;                  % points per solid wall
o.nvbd = prams.nvbd;                % number of solid walls
o.verbose = options.verbose;        % write data to console
o.saveData = options.saveData;      % save the data
o.usePlot = options.usePlot;        % plot the data
o.track = options.track;            % include tracker points
o.quiver = options.quiver;          % include velocity vectors
o.axis = options.axis;              % axis of plot
o.dataFile = options.dataFile;      % data file name
o.logFile = options.logFile;        % log file name
o.T = prams.T;                      % time horizon
o.m = prams.m;                      % number of time steps
o.errorTol = prams.errorTol;        % error tolerance
o.timeAdap = options.timeAdap;      % time adpativity
o.order = options.order;            % time stepping order
o.orderGL = options.orderGL;        % order of Gauss-Lobatto 
                                    % quadtrature
o.pressure = options.pressure;      % calculate pressure
o.rtolArea = prams.rtolArea;        % allowable error in area
o.rtolLength = prams.rtolLength;    % allowable error in length

o.adhesion = options.adhesion;
o.adStrength = prams.adStrength;
o.adRange = prams.adRange;

oc = curve;
[o.reducedArea,o.area,o.length] = oc.geomProp(X);
% area, length, and reduced area of initial shape

end % constructor: monitor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initializeFiles(o,X,sig,...
    eta,RS,Xwalls,Xtra,pressTar)
% initializeFiles(X,sig,eta,RS,Xwalls,Xtra,pressTar) does the initial
% writing of data to files and the console.  It first deletes any
% previous data and then writes the number of points, tracer initial
% conditions, pressure targets X and Xwalls are the vesicles and solid
% walls, Xtra and pressTar are initial conditions for the tracers and
% pressure/stress target locations

N = o.N; % points per vesicle
nv = o.nv; % number of vesicles
Nbd = o.Nbd; % points per solid wall
nvbd = o.nvbd; % number of solid walls
oc = curve;

[xx,yy] = oc.getXY(Xwalls);
% seperate x and y coordinates of solid walls
if o.saveData
  fid = fopen(o.dataFile,'w');
  fwrite(fid,[N;nv;Nbd;nvbd],'double');
  % write number of points and vesicles to data file
  fwrite(fid,[xx(:),yy(:)],'double');
  % write the solid walls to the data file
  fclose(fid);
  fid = fopen(o.logFile,'w');
  fclose(fid);
end

o.writeMessage(' ','%s\n')
message = ['************* PHYSICAL PARAMETERS *************'];
o.writeMessage(message,'%s\n')

message = [num2str(N) ' points per vesicle'];
o.writeMessage(message,'%s\n')
% write number of points

message = [num2str(nv) ' total vesicles'];
o.writeMessage(message,'%s\n')
% write number of vesicles

if ~o.timeAdap
  message = ['Time step size is ' ...
      num2str(o.T/o.m)];
  if o.orderGL > 0
    message = ['Monitoring residual using ' ...
      num2str(o.orderGL) ' Gauss-Lobatto points'];
  else
    message = ['Monitoring residual using ' ...
      num2str(-o.orderGL) ' equi-spaced points'];
  end
else
  if o.orderGL > 0
    message = ['Adaptive time stepping using ' ...
      num2str(o.orderGL) ' Gauss-Lobatto points'];
  else
    message = ['Adaptive time stepping using ' ...
      num2str(-o.orderGL) ' equi-spaced points'];
  end
end
o.writeMessage(message,'%s\n')
% write time step size or that we are doing
% adaptive time stepping

message = ['Time stepping order is ' ...
    num2str(o.order)];
o.writeMessage(message,'%s\n')
% write time step order

if o.timeAdap
  message = ['Allowable error in area is:   ' ...
      num2str(o.rtolArea,'%4.2e')];
  o.writeMessage(message,'%s\n')
  message = ['Allowable error in length is: ' ...
      num2str(o.rtolLength,'%4.2e')];
  o.writeMessage(message,'%s\n\n')
end
% write the allowable tolerances for the residual,
% error, and area

if o.adhesion
  message = ['Adhesion strength is:         ' ...
      num2str(o.adStrength,'%4.2e')];
  o.writeMessage(message,'%s\n')
  message = ['Adhesion range is:            ' ...
      num2str(o.adRange,'%4.2e')];
  o.writeMessage(message,'%s\n\n')
end

if o.saveData
  o.writeData(X,sig,0,0,0,0);
  % save initial configuartion

  message = ['Initial Areas are:            ' ...
      num2str(o.area(1),'%10.2e')];
  o.writeMessage(message,'%s\n')
  message = ['Initial Lengths are:          ' ...
      num2str(o.length(1),'%10.2e')];
  o.writeMessage(message,'%s\n')
  message = ['Initial Reduced Areas are:    ' ...
      num2str(o.reducedArea(1),'%10.2e')];
  o.writeMessage(message,'%s\n\n')
end
% write initial reduced area, area, and length to log file

o.writeStars
o.writeMessage(' ','%s\n')

end % initializeFiles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function terminate = outputInfo(o,X,sigma,u,eta,RS,...
    Xwalls,Xtra,time,iter,dtScale,res,iflag)
% computes the error in area and length and write messages to the data
% file, the log file, and the console.  Tells the simulation to stop if
% error in area or length is too large.  X, sigma, u are the position,
% tension, and velocity of vesicles, Xwalls is the parameterization of
% the solid walls, time is the current time, iter is the number of
% GMRES iterations, dtScale is the amount the time step is scaled by,
% res is the residual incurred from the previous time step to the
% current, and iflag is the success of gmres

errorTol = o.errorTol;
terminate = false;
[ea,el] = o.errors(X);
if max(ea,el) > errorTol
  message = 'ERROR IS TOO LARGE!  I AM STOPPING!';
  o.writeMessage(message,'%s\n');
  message = ['Max error in area is   ',num2str(ea,'%4.2e')];
  o.writeMessage(message,'%s\n');
  message = ['Max error in length is ',num2str(el,'%4.2e')];
  o.writeMessage(message,'%s\n');
  terminate = true;
  return
end
% check if error in area or length is greater than threshold

% Begin plotting
if o.usePlot
  o.plotData(X,time,ea,el);
  pause(0.01)
end
% End plotting

% Begin saving data
if o.saveData
  % don't want to save initial small time steps, but still
  % want to check the error in area and length so that the
  % simulation is killed early on if need be
  o.writeData(X,sigma,ea,el,time,res);
end
% End saving data


% Begin sending messages to log file and console
message = ['GMRES required ' num2str(iter) ...
    ' iterations to couple vesicles and solid walls'];
o.writeMessage(message,'%s\n')

if iflag == 1
  message = 'GMRES DID NOT CONVERGE: didn''t achieve tolerance';
elseif iflag == 2
  message = 'GMRES DID NOT CONVERGE: preconditioner ill-conditioned';
elseif iflag == 3
  message = 'GMRES DID NOT CONVERGE: successive iterates were the same';
end
if iflag ~=0
  o.writeStars
  o.writeMessage(message,'%s\n')
  o.writeStars
end

message = ['t = ' num2str(time,'%4.2e') ...
      ' of T = ' num2str(o.T,'%4.2e')]; 
o.writeMessage(message,'%s\n')

if o.timeAdap
  message = ['Time step size scaling ' num2str(dtScale,'%4.2e')];
  o.writeMessage(message,'%s\n')
end
% write the new time step size if doing time adaptivity

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeStars(o)
% writeStars writes a message of stars to the console and the log file
% depending on verbose and saveData

messageStars = '*********************************************';
o.writeMessage(messageStars,'%s\n')

end % writeStars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
function plotData(o,X,time,ea,el)
% plotData(X,time,ea,el) plots the current configuration with title X is
% the vesicle position time is the current time, ea and el are the
% errors in area and length

N = size(X,1)/2; % Number of points per vesicle
oc = curve;
[x,y] = oc.getXY(X);
% seperate x and y coordinates

figure(1); clf; hold on
plot([x;x(1,:)],[y;y(1,:)],'r','linewidth',2)
% Plot all vesicles

titleStr = ['t = ' num2str(time,'%4.2e') ...
  ' eA = ' num2str(ea,'%4.2e') ...
  ' eL = ' num2str(el,'%4.2e')];
title(titleStr)
% Title
axis equal
axis(o.axis)
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ycolor','w')
set(gca,'ycolor','w')

end % plotData

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeData(o,X,sigma,ea,el,time,res)
% writeData(X,sigma,ea,el,time,res) writes the position, tension,
% errors, and time to a binary file.  Matlab can later read this file to
% postprocess the data
 
oc = curve;
[x,y] = oc.getXY(X);
output = [x(:);y(:);sigma(:);ea;el;time];
% format that postProcess/loadfile.m reads the output
fid = fopen(o.dataFile,'a');
fwrite(fid,output,'double');
fclose(fid);

end % writeData


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function summary(o,matvecs,derivs,fmms,...
        naccept,nreject,countGMRES,totTime)
% summary(matvecs,derivs,fmms,naccept,nreject,totTime) writes a final
% summary if the simulation was successful

messageStar = '******************************************';
message = ['Number of times in TimeMatVec is ' ...
    num2str(matvecs)];
o.writeStars
o.writeMessage(message,'%s\n')
o.writeStars
% print the total number of matvecs

message = ['Number of times in computeDeriv is ' ...
    num2str(derivs)];
o.writeStars
o.writeMessage(message,'%s\n')
o.writeStars
% print the total number of times differential operators are
% formed in matrix form

message = ['Number of times calling Stokes FMM is ' ...
    num2str(fmms)];
o.writeStars
o.writeMessage(message,'%s\n')
o.writeStars
% print the total number of fmm calls
  
message = ['Number of total GMRES iterations is ' ...
    num2str(countGMRES)];
o.writeStars
o.writeMessage(message,'%s\n')
o.writeStars
%print the total number of GMRES iterations

if o.timeAdap  
  message = ['Number of accepted time steps is ' num2str(naccept)];
  o.writeStars
  o.writeMessage(message,'%s\n')
  message = ['Number of rejected time steps is ' ...
      num2str(nreject)];
  o.writeStars
  o.writeMessage(message,'%s\n')
  o.writeStars
end
% print the total number of accepted and rejected time steps 

message = ['Total time is ' num2str(totTime,'%5.2e') ];
o.writeMessage(message,'%s\n')
o.writeStars
% Save total time spent

end % summary


end % methods


end % classdef

