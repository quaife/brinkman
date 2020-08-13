classdef tstep < handle
% This class implements the velocity evaulation schemes and the time
% stepping methods


properties

ves; % vesicle structure
shearRate; % background shear rate
bendsti; % maximum bending stiffness
bendratio; % ratio between max and min bending stiffness
kmatrix; % matrix for doing odd-even integration
viscIn;
viscOut;
gmresTol;
gmresMaxIter;

end % properties

methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = tstep(params,ves)
%o = tstep(prams,ves) is a constructor that initializes the class. 
%Take all elements of prams needed by the time stepper
o.ves = ves;
o.shearRate = params.shearRate;
o.bendsti = params.bendsti;
o.bendratio = params.bendratio;
o.viscIn = params.viscosityInside;
o.viscOut = params.viscosityOutside;
op = poten(params.N);
o.kmatrix = op.oddEvenMatrix;
o.gmresTol = params.gmresTol; %GMRES tolerance
o.gmresMaxIter = params.gmresMaxIter; %maximum number of GMRES iterations

end % tstep: constructor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [un, ut] = usetself(o)
% Compute x and y component of the velocity. Naming convention needs
% fixed

ves = o.ves;
x0 = ves.x0;
y0 = ves.y0;

N = ves.N;
L = ves.L;
cur = ves.cur;
theta = ves.theta;
op = poten(N);
oc = curve;

%ves.X = oc.recon(N,x0,y0,L,theta);
[Eu,Esigma] = ves.variationsNonStiff;

%Compute the force in eq (33)
tau = [[+Esigma.*sin(theta) - Eu.*cos(theta)]; ...
       [-Esigma.*cos(theta) - Eu.*sin(theta)]];
% MISSING u_s TERM???
%  plot(tau)
%  pause
%Construct Stokes matrix without the log singularity. ie. A3 only
%contains the rightmost kernel in equation (43)
StokesMat = op.StokesMatrixLogless(ves.X);
%form the velocity, k, on the interface corresponding to v^u in eq (33)
%If [[P^u n]]_sigma = tau, then k = stokesMatrix*tau
k = StokesMat*tau;
%ulam is the viscosity contrast
ulam = ves.viscIn/ves.viscOut;
%calculate constants that multiply the weakly singular and regular parts of 
%the integral operators
%c1 = L/pi/N/(1+ulam)/ves.viscOut/2;
%c2 = L/(1+ulam)/ves.viscOut;
c1 = 1/(4*pi)*L/N;
c2 = -L/(8*pi);

%LogKernel1 and LogKernel2 are the log kernels in the left term in the
%right hand side of equation (43) integrated against the x and y
%components of the density functions tau that involves Eu and Esigma
LogKernel1 = op.IntegrateLogKernel(tau(1:N));
LogKernel2 = op.IntegrateLogKernel(tau(N+1:end));

%sigma1 and sigma2 are the solution of equation (33) (still unsure
%about the u_s term). Also, it includes the background shear flow which
%is not stated in equation (33), but instead in equations (6) and (7)
sigma1 = k(1:N)*c1 + LogKernel1*c2 + ves.X(N+1:end)*o.shearRate;
sigma2 = k(N+1:end)*c1 + LogKernel2*c2;

% Calculate v dot n in eq (40)
vdotn = sigma1.*sin(theta) - sigma2.*cos(theta);
% Calculate v dot s in eq (40)
vdots = sigma1.*cos(theta) + sigma2.*sin(theta);

%Calculate d/ds(v dot s) in eq (40)
IK = oc.modes(N);
dvdotsds = oc.diffFT(vdots,IK)/L;
%pause
% compute the initial curvature using the tangent angle
%ves.cur = oc.acurv(ves);
% NOTE: make sure this is the actual curvature. Compare with the
% analytic expression for the ellipse that can easily be computed.
%Calculate the right hand side of equation (40)
% plot(dvdotsds/L)
% hold on
% plot(vdotn.*ves.cur)
% pause
% [dvdotsds/L vdotn.*ves.cur]
rhs = -(dvdotsds + ves.cur.*vdotn);
% disp('plotting rhs')
% plot(rhs)
% pause
%The velocity components in eq (40) are nonlocal linear functions of 
%lambdaTilde. 
%semilogy(abs(fftshift(fft(rhs)))/numel(rhs))
% figure(1)
% disp('plotting crv')
%  plot(ves.cur)
% disp('plotting theta')
% figure(2)
%  plot(ves.theta)
%  pause(1)

%Solve the linear system for LambdaTilde in (39) using GMRES.
%Each iteration of GMRES requires a solution of Stokes equation
%LambdaTilde is the lambda with a tilde in eq (39)
% [lambTil,flag,relres,iter,resvec] = ...
%       gmres(@(x) o.matvec40(x,StokesMat),rhs,[],o.gmresTol,...
%           o.gmresMaxIter);
%flag
%relres
%figure(1); clf;
%%plot(rhs)
%figure(2); clf
%plot(StokesMat*[rhs;rhs])
%surf(StokesMat)
%shading interp
%plot(o.matvec40(lambTil,StokesMat) - rhs)
%plot(o.matvec40(rhs,StokesMat))
%%plot(lambTil)
%pause
%[lambTil,flag,relres,iter,resvec] = ...
%      gmres(@(x) o.matvec40(x,StokesMat),rhs,[],o.gmresTol,...
%              o.gmresMaxIter,@(x) o.preconditioner(x));
[lambTil,flag,relres,iter,resvec] = ...
      gmres(@(x) o.matvec40(x,StokesMat),rhs,[],o.gmresTol,...
              o.gmresMaxIter);
%       flag
%       relres
%       iter
%       pause
%calculate the Fourier derivative of lambdaTilde
dlamTilds = oc.diffFT(lambTil,IK)/L;
% plot(dlamTilds)
% pause
% We can now calculate the traction jump in first part of equation (39).
% This comes from applying the product rule and using Frenet-Seret.
tracJump = [[-lambTil.*ves.cur.*sin(theta) + dlamTilds.*cos(theta)];...
            [+lambTil.*ves.cur.*cos(theta) + dlamTilds.*sin(theta)]];

% Adding the jump conditions in eq (39) to (33) which is in the variable
% tau
tau = -tau - tracJump;
%   disp('HERE2')
%   plot(tau)
%   pause
%Calculating u tilde in equations (38) through (40) without the weakly
%singular log kernel
k = StokesMat*tau;

%compute the weakly singluar log kernel part
force1 = op.IntegrateLogKernel(tau(1:N));
force2 = op.IntegrateLogKernel(tau(N+1:end));
%  disp('HERE')
%  plot(force1)
%  hold on
%  plot(force2)
%  pause
% % plot([force1 force2])
% pause
%Calulating  u in eqatuion (31) by adding the results from the
%non-singular and weakly singular integral operators
% un = k(1:N)*c1 - tau(1:N)*c2 + ves.X(N+1:end)*o.shearRate;
% ut = k(N+1:end)*c1 - tau(N+1:end)*c2;
un = k(1:N)*c1 + force1*c2 + ves.X(N+1:end)*o.shearRate;
ut = k(N+1:end)*c1 + force2*c2;

%  disp('Un Ut')
%  plot(un)
%  hold on
%  plot(ut)
%  pause

end % usetself

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = matvec40(o,rhs,StokesMat)   
%This function cooresponds to the matvec in equation 40 and returns
%%(u \cdot s)_s + kappa * (u \cdot n)  as x

ves = o.ves;
theta = ves.theta;
cur = ves.cur;
L = ves.L;
N = ves.N;
ulam = ves.viscIn/ves.viscOut;
oc = curve;
op = poten(N);
IK = oc.modes(N);
%The the derivative of the rhs of eq(40)
drhsds = oc.diffFT(rhs,IK)/L;
%Compute the forces
tau1 = -rhs.*cur.*sin(theta) + drhsds.*cos(theta);
tau2 = +rhs.*cur.*cos(theta) + drhsds.*sin(theta);
%form the velocity on the interface
krhs = StokesMat*[tau1;tau2];
%LogKerneltau1 and LogKerneltau2 are the log kernels integrated against the 
%density functions tau1 and tau2.
LogKerneltau1 = op.IntegrateLogKernel(tau1);
LogKerneltau2 = op.IntegrateLogKernel(tau2);
%calculate constants that multiply the weakly singular and regular parts of 
%the integral operators
%c1 = L/pi/N/(1+ulam)/ves.viscOut/2;
%c2 = L/(1+ulam)/ves.viscOut;
c1 = 1/(4*pi)*L/N;
c2 = -L/(8*pi);

% (sigma1,sigma2) is \tilde{u} in equation (39) and (40)
sigma1 = krhs(1:N)*c1 + LogKerneltau1*c2;
sigma2 = krhs(N+1:end)*c1 + LogKerneltau2*c2;

% Build left hand side of equation (40)
vdotn = sigma1.*sin(theta) - sigma2.*cos(theta);
vdots = sigma1.*cos(theta) + sigma2.*sin(theta);
dvdotsds = oc.diffFT(vdots,IK)/L;

%x is (u \cdot s)_s + kappa * (u \cdot n) in eq ()
x = (dvdotsds + cur.*vdotn);

end %matvec40



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = preconditioner(o,rhs)
% preconditioner (equation (63) and (64)) for equation (40) which is the
% Schur complement for a Lagrange multiplier needed to enforce local
% inextensibility
ves = o.ves;
N = ves.N;
x = fft(rhs,N)./[1 1:1:N/2 N/2-1:-1:1]';
x = real(ifft(x,N));
end %preconditioner

end % methods

end % classdef
