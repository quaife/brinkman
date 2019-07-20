%function [] = extensional2Ves

clear all; clc

fprintf('Two elliptical vesicles in a extensional flow.\n');

% Physics parameters
prams.N = 64;               % points per vesicle
prams.nv = 2;               % number of vesicles
prams.T = 100;               % time horizon (two tumbling)
prams.m = 100;             % number of time steps
prams.kappa = [1 1];         % bending coefficient
prams.viscCant = [1 1];         % viscosity contrast
options.farField = 'extensional'; % background velocity
options.farFieldSpeed = 0.07;
aptions.order = 1;          % time stepping order
options.vesves = 'implicit';
% Discretization of vesicle-vesicle interactions.
% Eithe 'explicit' or 'implicit'
options.inextens = 'method1';
options.near = true;        % near-singular integration
options.fmm = false;
options.verbose = true;
options.antiAlias = false;
prams.gmresMaxIter = 3*prams.N;
prams.gmresTol = 1e-6;
prams.errorTol = 1;

% ADD-ONS
options.correctShape = false;
options.adhesion = false;
prams.adRange = 4e-1;
prams.adStrength = 7e-1;

% TIME ADAPTIVITY (parameters for new implementation)
options.timeAdap = false;

prams.rtolArea = 1e-2;
prams.rtolLength = 1e-2;
prams.betaUp = 1.2;
prams.betaDown = 0.5;
prams.alpha = 0.9;

options.orderGL = 2;
options.nsdc = 0;
options.expectedOrder = 1;

% Plot on-the-fly
options.usePlot = true;
options.axis = [-4 4 -2 2];
options.track = false;
% Save vesicle information and create a log file
options.logFile = 'output/extensional2Ves.log';
% Name of log file for saving messages
options.dataFile = 'output/extensional2VesData.bin';
% Name of binary data file for storing vesicle information

options.saveError = true;
options.errorFile = 'output/extensional2VesError.bin';
% Name of binary data file for storing truncation errors after each step

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

posx1 = load('posx1_RA070_rotated.dat');
posy1 = load('posy1_RA070_rotated.dat');
posx2 = load('posx2_RA070_rotated.dat');
posy2 = load('posy2_RA070_rotated.dat');


if 0
ysep = mean(posy2) - mean(posy1);
posx1 = posx1 - mean(posx1);
posx2 = posx2 - mean(posx2);
posy1 = posy1 - mean(posy1);
posy2 = posy2 - mean(posy2);
% center everything at the origin so that the problem is symmetric

%ysep = 0.60; % approximate mean value of the minimum seperation
posy1 = posy1 - ysep/2;
posy2 = posy2 + ysep/2;
end

sig = 2e-2;
modes = (-64:63)';
gauss = exp(-2*sig^2*pi^2*modes.^2);
% Gaussian filter in Fourier space

posx1 = real(ifft(ifftshift(fftshift(fft(posx1)).*gauss)));
posx2 = real(ifft(ifftshift(fftshift(fft(posx2)).*gauss)));
posy1 = real(ifft(ifftshift(fftshift(fft(posy1)).*gauss)));
posy2 = real(ifft(ifftshift(fftshift(fft(posy2)).*gauss)));
% apply a Gaussian filter to smooth out the initial shape

posx1 = interpft(posx1,prams.N);
posx2 = interpft(posx2,prams.N);
posy1 = interpft(posy1,prams.N);
posy2 = interpft(posy2,prams.N);

X = [posx1 posx2; posy1 posy2];

oc = curve;
centerx = [0 3];
centery = zeros(1,2);
ang = pi/2*ones(2,1);
ra = 0.95;
scale = 0.5*sqrt(ra);
X = oc.initConfig(prams.N,'nv',prams.nv,...
    'reducedArea',ra,...
    'center',[centerx;centery],...
    'angle',ang,...
    'scale',scale);
% Initial configuration of reduce area 0.65 and aligned

Xfinal = Ves2D(X,[],prams,options);
% Run vesicle code

