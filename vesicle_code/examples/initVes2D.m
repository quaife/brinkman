function [options,prams] = initVes2D(options,prams)
% Set a path pointing to src directory and set options and
% prams to default values if they are left blank

P = path; ii = find(pwd == filesep); ii = ii(end);
subPath = pwd; subPath = [subPath(1:ii) 'src'];
if isempty(strfind(P, subPath)),addpath(subPath);end


PramList = {'N','nv','T','m','Nbd','nvbd','kappa','viscCont',...
    'gmresTol','gmresMaxIter','errorTol','rtolArea','rtolLength',...
    'betaUp','betaDown','alpha','adStrength','adRange',...
    'fluxCoeff','fluxShape','dtMax','dtMin','saveRate'};
defaultPram.N = 64;
defaultPram.nv = 1;
defaultPram.T = 1;
defaultPram.m = 100;
defaultPram.Nbd = 0;
defaultPram.nvbd = 0;
defaultPram.kappa = 1e-1*ones(1,prams.nv);
defaultPram.viscCont = ones(1,prams.nv);
defaultPram.gmresTol = 1e-12;
defaultPram.gmresMaxIter = 200;
defaultPram.errorTol = 1e-1;
defaultPram.rtolArea = 1e-2;
defaultPram.rtolLength = 1e-2;
defaultPram.betaUp = 1.2;
defaultPram.betaDown = 0.5;
defaultPram.alpha = 0.9;
defaultPram.adStrength = 4e-1;
defaultPram.adRange = 8e-1;
defaultPram.fluxCoeff = 0; 
defaultPram.fluxShape = zeros(prams.N,prams.nv); 
defaultPram.dtMax = 1e0;
defaultPram.dtMin = 1e-5;
defaultPram.saveRate = 1;

for k = 1:length(PramList)
  if ~isfield(prams,PramList{k})
    eval(['prams.' PramList{k} '=defaultPram.' PramList{k} ';'])
    % Set any unassigned parameters to a default value
  end
end

OptList = {'expectedOrder','farField',...
    'farFieldSpeed','near',...
    'fmm','fmmDLP','confined','usePlot',...
    'axis','saveData','logFile','dataFile','verbose',...
    'profile','collision','timeAdap',...
    'pressure','SDCcorrect','orderGL','nsdc','adhesion','expForce',...
    'fmmPrecision','semipermeable','fluxShape',...
    'xshift','xshiftLoc','xshiftVec'};

defaultOpt.expectedOrder = 2;
defaultOpt.farField = 'shear';
defaultOpt.farFieldSpeed = 1;
defaultOpt.near = true;
defaultOpt.fmm = false;
defaultOpt.fmmDLP = false;
defaultOpt.confined = false;
defaultOpt.usePlot = true;
defaultOpt.axis = [-5 5 -5 5];
defaultOpt.saveData = true;
defaultOpt.logFile = 'output/example.log';
defaultOpt.dataFile = 'output/exampleData.bin';
defaultOpt.verbose = true;
defaultOpt.profile = false;
defaultOpt.collision = false;
defaultOpt.timeAdap = false;
defaultOpt.pressure = false;
defaultOpt.SDCcorrect = false;
defaultOpt.orderGL = 2;
defaultOpt.nsdc = 0;
defaultOpt.adhesion = false;
defaultOpt.expForce = false;
defaultOpt.fmmPrecision = 4;
defaultOpt.semipermeable = false;
defaultOpt.fluxShape = 1;
defaultOpt.xshift = false;
defaultOpt.xshiftLoc = 0;
defaultOpt.xshiftVec = 0;

for k = 1:length(OptList)
  if ~isfield(options,OptList{k})
    eval(['options.' OptList{k} '=defaultOpt.' OptList{k} ';'])
    % Set any unassigned options to a default value
  end
end

if ~options.confined
  prams.Nbd = 0;
  prams.Nbdcoarse = 0;
  prams.nvbd = 0;
end
% If the geometry is unbounded, make sure to set the number
% of points and number of components of the solid walls
% to 0.  Otherwise, later components will crash

if numel(prams.kappa) ~=prams.nv
  prams.kappa = prams.kappa*ones(1,prams.nv);
end
if numel(prams.viscCont) ~=prams.nv
  prams.viscCont = prams.viscCont*ones(1,prams.nv);
end


