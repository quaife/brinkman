function stokesSolver(Xinner,Xouter,options,prams)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% tracers = stokesSolver(Xouter,Xinner,options,prams) returns the time
% history of a set of tracers
% INPUTS
% Xouter - parameterization of the outer boundary
% Xinner - parameterization of the inner boundaries
% prams - set of parameters
% options - set of options
% OUTPUTS:
% none 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global iteration
iteration = 0;
om = monitor(options,prams);
% object for doing I/O
om.welcomeMessage(options,prams);
% output parameters
om.clearFiles;
% empty the dat and bin files

innerGeom = bodies(Xinner);
outerGeom = bodies(Xouter);
clear Xinner
clear Xouter
% create objects for the inner and outer boundaries.
% The outer boundary will need more points so need
% two classes.

rhs = [zeros(2*innerGeom.N*innerGeom.nv,1); outerGeom.u(:)];
% right-hand side which corresponds to no-slip on the solid walls
%rhs = [innerGeom.u(:);outerGeom.u(:)];
%% right-hand side with an analytic solution which is a pipe flow
%% everywhere

op = poten(innerGeom,options.fmm,true,false);
% object for evaluating layer potentials
[~,NearO2I] = outerGeom.getZone(innerGeom,2);
[NearI2I,NearI2O] = innerGeom.getZone(outerGeom,3);

% single-layer potential of the inner geometry
D = op.stokesDLmatrix(outerGeom);
N0 = op.stokesN0matrix(outerGeom);
OuterPreco = inv(-1/2*eye(2*outerGeom.N) + D + N0);
% exact inverse of double-layer potential term for the outer boundary

tic
%[eta,flag,relres,iter,relresvec] = gmres(...
%    @(X) op.matVecMultiply(X,innerGeom,outerGeom,...
%    NearI2I,NearI2O,NearO2I),...
%    rhs,[],prams.gmresTol,prams.maxIter);
[eta,flag,relres,iter,relresvec] = gmres(...
    @(X) op.matVecMultiply(X,innerGeom,outerGeom,...
    NearI2I,NearI2O,NearO2I),...
    rhs,[],prams.gmresTol,prams.maxIter,...
    @(X) op.matVecInvBD(X,innerGeom,OuterPreco));

om.writeStars
message = ['****     GMRES took ' num2str(toc,'%4.2e') ...
    ' seconds     ****'];
om.writeMessage(message,'%s\n');
if (flag ~= 0)
  message = 'GMRES HAD A PROBLEM';
  om.writeMessage(message,'%s\n');
  message = ['flag = ' num2str(flag)];
  om.writeMessage(message,'%s\n');
  message = ['Relative residual = ' num2str(relres,'%4.2e')];
  om.writeMessage(message,'%s\n');
  message = ['Inner iterations = ' num2str(iter(2))];
  om.writeMessage(message,'%s\n');
else
  message = ['****     GMRES required ' num2str(iter(2),'%3d'),...
      ' iterations    ****'];
  om.writeMessage(message,'%s\n');
end
%om.writeMessage(message,'%s\n');
om.writeStars
om.writeMessage(' ');

sigmaInner = zeros(2*prams.Ninner,prams.nv);
% initialize space for desnity function along inner boundaries
for k = 1:prams.nv
  istart = 2*(k-1)*prams.Ninner + 1;
  iend = istart + 2*prams.Ninner - 1;
  sigmaInner(:,k) = eta(istart:iend);
end
% unstack the density function at the inner boundaries

istart = 2*prams.nv*prams.Ninner + 1;
iend = istart + 2*prams.Nouter - 1;
sigmaOuter = eta(istart:iend);
% unstack the density function at the outer boundary

om.writeGeometry(innerGeom.X,outerGeom.X,sigmaInner,sigmaOuter);
% write the density function to the data file for post processing
om.writeStars
message = '****    Density function and geometry    ****';
om.writeMessage(message);
message = '****         written to bin file         ****';
om.writeMessage(message);
om.writeStars
om.writeMessage(' ');



