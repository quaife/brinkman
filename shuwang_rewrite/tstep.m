classdef tstep < handle
% This class implements the velocity evaulation schemes and the time
% stepping methods


properties

ves; % vesicle structure
shearRate; % background shear rate
bendsti; % maximum bending stiffness
bendratio; % ratio between max and min bending stiffness
kmatrix; % matrix for doing odd-even integration



end % properties

methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = tstep(params,ves)
  
o.ves = ves;
o.shearRate = params.shearRate;
o.bendsti = params.bendsti;
o.bendratio = params.bendratio;
op = poten(params.N);
o.kmatrix = op.oddEvenMatrix;




end % tstep: constructor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = usetself(o)

ves = o.ves;
X = ves.X;

[Eu,Esigma] = ves.variationsNonStiff;

end % usetself


end % methods

end % classdef
