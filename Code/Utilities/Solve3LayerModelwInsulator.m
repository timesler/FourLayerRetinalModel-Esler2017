%%SOLVE3LAYERMODELWINSULATOR
% SOLVE3LAYERMODELWINSULATOR use MATLAB's symbolic mathematics toolbox to
% derive expressions for the integration coefficients used to solve the
% system of partial differential equations which define extracellular
% potential in the four layer model of the retina. See
% CellComp4Layer_Ve_f_Plane for an implementation of these coefficients.
%
% This script was used to conduct the analysis presented in:
%
%	T. Esler, R.R. Kerr, B. Tahayori, D.B. Grayden, H. Meffin, and A.N. 
%	Burkitt (2017), "Minimizing activation of overlying axons with epiretinal
%	stimulation: The role of fiber orientation and electrode configuration",
%	PLOS ONE.
%
%%%%%%%%%%%%%%%%%%%%%%%%% Created by: Tim Esler, 2017 %%%%%%%%%%%%%%%%%%%%%%%%%

%% Define symbols

syms eta_V eta_F eta_L
syms Yii h_F d_E
syms sigma_V xi_T_f sigma_L_y
syms m

%% Define linear system

A = [1, 1, -1, -1, 0; ...
    0, 0, exp(eta_F*h_F), exp(-eta_F*h_F), -exp(-eta_L*h_F); ...
    sigma_V*eta_V, -sigma_V*eta_V, -xi_T_f*eta_F, xi_T_f*eta_F, 0; ...
    0, 0, xi_T_f*eta_F*exp(eta_F*h_F), -xi_T_f*eta_F*exp(-eta_F*h_F), sigma_L_y*eta_L*exp(-eta_L*h_F)
    sigma_V*eta_V*exp(-eta_V*(Yii+d_E)), -sigma_V*eta_V*exp(eta_V*(Yii+d_E)), 0, 0, 0];

B = [-m/(2*eta_V)*exp(-eta_V*Yii); ...
    0; ...
    m*sigma_V/2*exp(-eta_V*Yii); ...
    0; ...
    -m*sigma_V/2*exp(-eta_V*d_E)];

%% Solve and simplify

C = A\B;

C = simplify(C,'Steps',1000);
C = simplify(C,'Steps',1000);
C = simplify(C,'Steps',1000);

C = limit(C, d_E, 0, 'right');

C = simplify(C,'Steps',1000);
C = simplify(C,'Steps',1000);
C = simplify(C,'Steps',1000);
