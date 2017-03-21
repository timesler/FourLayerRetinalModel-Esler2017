function [p] = NTESparams(mode)
%%NTESPARAMS
% NTESPARAMS returns the neurite dimensions and physiological parameters
% defined below. Parameters can either be returned as double or single
% precision floats, using mode = 'double' and mode = 'single',
% respectively.
%
%%%%%%%%%%%%%%%%%%%%%%%%% Created by: Tim Esler, 2017 %%%%%%%%%%%%%%%%%%%%%%%%%

p.b = 0.5e-6;                       % NTES radius (m)
p.d = 0.03e-6;                      % Width of extracellular sheath (m)
p.a = p.b-p.d;                      % Neurite radius (m)

p.C_m = 0.01;                       % Membrane capacitance (F/m^2)
p.R_m = 1;                          % Membrane unit area resistance (ohm.m^2)

p.rho_i = 0.7;                      % Intracellular resistivity (ohm.m)
p.rho_e = 0.7;                      % Extracellular resistivity (ohm.m)
p.r_m = p.R_m/(2*pi*p.a);           % Membrane unit length resistance (ohm.m)
p.r_i = p.rho_i/(pi*p.a^2);         % Intracellular resistance (ohm/m)
p.r_e = p.rho_e/(pi*(p.b^2-p.a^2)); % Extracellular resistance (ohm/m)

if strcmpi(mode,'single')
    p = structfun(@single,p,'UniformOutput',0);
end

end