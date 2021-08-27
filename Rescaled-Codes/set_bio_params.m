function [a,R,gk,ek,gna,ena,gl,el,ni,mi,hi,C,vstart,b]=set_bio_params()
a   = 0.05*1e-6;            % radius                       [m]
R   = 300*1e-2;               % resistance                   [ohm.m]
C   = 0.75e-2;                  % capacitance                  [F/m2]
gk  = 5*1e1;                 % potassium ion conductance    [S/m2]
ek  = -90*1e-3;              % potassium reversal potential [V]
gna = 50*1e1;                % sodium ion conductance       [S/m2]
ena = 50*1e-3;               % sodium reversal potential    [V]
gl  = 0.0*1e1;               % leak conductance             [S/m2]
el  = -70*1e-3;              % leak potential               [m]
vstart = 10*1e-3;            % clamp voltage                [V]

%compute diff coefficient
b=a/(2*R*C);

ni = 0.0376969;   % these are ion gating variables dimensionless
mi = 0.0147567;   % these are ion gating variables dimensionless
hi = 0.995941;    % these are ion gating variables dimensionless
end

