function ves = Ves2D(ves,alpha,params,options)
% Ves2D does time stepping on the intial configuration X
% with parameters and options defined in prams and options.

tt = tstep(params,ves); %set up tstep class
om = monitor(ves.X,params,options); %set up monitor class
%Take first step
[ves,ux_old,uy_old,L,Ln,dcur0,fntheta,N2Hat] = tt.FirstSteps(ves,params,...
                                                        options,om);
%Begin time step routine
ves = tt.TimeStepLoop(ves,params,om,ux_old,uy_old,L,Ln,dcur0,...
                      fntheta,N2Hat);
end %Ves2D