function [Ve_L_f,Ve_T_X_f,Ve_T_Y_f,Ve_T_Z_f,neur_eq] = CellComp3Layer_Ve_f_Plane(Xi, Yi, Zi, I_M, I_D, x_max, z_max, t_max, d_x, d_z, d_t, h_F, Ya, rot, Ri)
%%CELLCOMP3LAYER_VE_F_PLANE
% CELLCOMP3LAYER_VE_F_PLANE returns the electric field in the Fourier
% domain due to epiretinal stimulation with point or disk electrodes in a
% plane in the x-z dimensions (where y is in the direction normal to the
% surface of the retina, toward the electrode). The retina is modelled
% using a 3-layer description, including Vitreous, Nerve Fibre
% Layer, and the Ganglion Cell Layer (+ others)
% (see diagram below for geometry).
%
% All units are S.I.
%       
%                   [|||]        [|||]        [|||]
%                              Electrodes
%
%       1.Vitreous
%       _____________________________________________________
%       2.Nerve fibre layer
%       _____________________________________________________
%       3.'Other' cell layer (incl. GCL)
%       _____________________________________________________
%
%
% The extracellular potential for the NFL is calculated using a modified version
% of the self-consistent, linear, sub-threshold model presented in:
%
%   B. Tahayori, H. Meffin, E.N. Sergeev, I.M.Y. Mareels, A.N. Burkitt, and
%   D.N. Grayden (2014), "Modelling extracellular electrical stimulation:
%   IV. Effect of the cellular composition of neural tissue on its
%   spatio-temporal filtering properties", J. Neural Eng. 11.
%
% This script was used to conduct the analysis presented in:
%
%	T. Esler, R.R. Kerr, B. Tahayori, D.B. Grayden, H. Meffin, and A.N. 
%	Burkitt (2017), "Minimizing activation of overlying axons with epiretinal
%	stimulation: The role of fiber orientation and electrode configuration",
%	PLOS ONE.
%
% INPUTS:
%
% Xi, Yi and Zi         the x-, y-, and z- coordinates of the point source
%                       electrodes. (m)
% I_M, I_D              the stimulation amplitude and duration for each
%                       electrode. (A, s)
% z_max, x_max, t_max   z-, x- and time-extent of the simulation. (m, m, s)
% d_z, d_x, d_t         z , x and t step sizes. (m, m, s)
% h_F                   Nerve fibre layer thickness (m)
% Ya                    Depth below the retina surface at which we want to 
%						calculate extracellular potential. (m)
% rot                   coordinate rotation for a fibre which is rotated w.r.t.
%						the fibre bindle (only set to a nonzero value if the
%                       output will be used to calculate the membrane
%                       potential of a rotated fiber). (rad)
% Ri                    Radius of disk for each electrode (set to 0 for
%                       point source). (m)
%
% OUTPUTS:
%
% Ve_L_f                The calculated longitudinal extracellular potential for
%                       a full plane in the Fourier domain. This is equal to
%						the extracellular potential itself. (V)
% Ve_T_#_f              The calculated transverse extracellular potential
%                       in the #-direction for a full plane in the Fourier
%                       domain. (V)
% neur_eq               The transformation matrix for converting the
%                       longitudinal extracellular potential into
%                       longitudinal membrane potential. Multiplying 
%						Ve_L_f by this pointwise will give the membrane
%						potential in the x, z, and t Fourier domain.
%
%
% EXAMPLE USAGE:
%
% Xi = 0e-6;
% Yi = 200e-6;
% Zi = 0e-6;
% I_M = -50e-6;
% I_D = 200e-6;
% x_max = 2000e-6;
% z_max = 2000e-6;
% t_max = 600e-6;
% d_x = 10e-6;
% d_z = 10e-6;
% d_t = 4e-6;
% h_F = 100e-6;
% Ri = 50e-6;
% Ya = -10e-6;
%
% [Ve_L_f, ~, ~, ~, n_e] = CellComp3Layer_Ve_f_Plane(...
% 		Xi, Yi, Zi, I_M, I_D, ...                   % Electrode parameters
% 		x_max, z_max, t_max, d_x, d_z, d_t, ...     % Spatial sampling
% 		h_F, Ya, rot, Ri);							% Simulation geometry
%
%%%%%%%%%%%%%%%%%%%%%%%%% Created by: Tim Esler, 2017 %%%%%%%%%%%%%%%%%%%%%%%%%

%% Error checking

if (length(Xi)-length(Yi)+1)*(length(Xi)-length(Zi)+1)*...
        (length(Xi)-length(I_M)+1)*(length(Xi)-length(I_D)+1) ~= 1
    error('Inconsistent number of electrodes in electrode specifications (Xi, Yi, Zi, I_M, I_D)');
end

%% Define parameters

p = NTESparams('double');         % Call parameter function

% Unpack parameters
b = p.b;                % NTES radius (m)
d = p.d;                % Width of extracellular sheath (m)
a = p.a;                % Neurite radius (m)

C_m = p.C_m;            % Membrane capacitance (F/m^2)
R_m = p.R_m;            % Membrane unit area resistance (ohm.m^2)

rho_i = p.rho_i;        % Intracellular resistivity (ohm.m)
rho_e = p.rho_e;        % Extracellular resistivity (ohm.m)
r_m = p.r_m;            % Membrane unit length resistance (ohm.m)
r_i = p.r_i;            % Intracellular resistance (ohm/m)
r_e = p.r_e;            % Extracellular resistance (ohm/m)

%% Define sampling in time-space and Fourier domains

% Sampling space Fourier domains
kz_max = pi/d_z;
nz = (fix(z_max/d_z));
d_kz = kz_max/nz;
kzp = -kz_max:d_kz:kz_max;
kzp(fix(length(kzp)/2+1)) = kzp(fix(length(kzp)/2+1))/1000;

kx_max = pi/d_x;
nx = (fix(x_max/d_x));
d_kx = kx_max/nx;
kxp = -kx_max:d_kx:kx_max;
kxp(fix(length(kxp)/2+1)) = kxp(fix(length(kxp)/2+2))/1000;

% Sampling time Fourier domain
w_max = pi/d_t;
nt = (fix(t_max/d_t));
d_w = w_max/nt;
w = -w_max:d_w:w_max;
w(fix(length(w)/2+1)) = 1e-20;

% Sampling space domain
Zp = -z_max:d_z:z_max;

Xp = -x_max:d_x:x_max;

% Sampling time domain
T = -t_max:d_t:t_max;

% Create sample mesh
[kxp_m,kzp_m] = ndgrid(kxp,kzp);

% Apply rotation in Fourier space
kx_m = kxp_m*cos(rot) + kzp_m*sin(rot);
kz_m = -kxp_m*sin(rot) + kzp_m*cos(rot);

Ve_T_X_f = zeros(length(Xp),length(Zp),length(T));
Ve_T_Y_f = zeros(size(Ve_T_X_f));
Ve_T_Z_f = zeros(size(Ve_T_Y_f));
Ve_L_f = zeros(size(Ve_T_Y_f));
neur_eq = zeros(size(Ve_T_Y_f));

for j = 1:length(w)
    w_m = w(j);
    
    %% Define electrotonic length constants, time constants
    % in the Fourier domain
    
    tau_m = R_m*C_m;                % Membrane time constant (s)
    
    % Electrotonic length constants (static and frequency-dependent, for both
    % current density (J) and voltage (V) boundary conditions)
    L_0J = sqrt(r_m/(r_e+r_i));
    L_0V = sqrt(r_m/r_i);
    L_J_m = L_0J./sqrt(1+1i*w_m*tau_m);
    L_V_m = L_0V./sqrt(1+1i*w_m*tau_m);
    
    %% Define admittivities and conductivities
    
    % Vitreous layer
    sigma_V = 1.78;                        % Ohm.m
    
    % Lower layers
    sigma_L_xz = 0.1;                     % Ohm.m
    sigma_L_y  = 0.1;
    
    % Nerve fibre layer (admittivity)
    xi_L_f_m = 1/rho_i * (1+(kz_m.^2).*L_J_m.^2)./(1+(kz_m.^2).*L_V_m.^2);
    xi_T_f = d/b/rho_e;
    chi_f_m = sqrt(xi_L_f_m/xi_T_f);    % Anisotropy ratio
    
    %% Define eta coefficients
    
    eta_V = sqrt((kx_m.^2*sigma_V    + kz_m.^2*sigma_V   )./sigma_V  );
    eta_L = sqrt((kx_m.^2*sigma_L_xz + kz_m.^2*sigma_L_xz)./sigma_L_y);
    eta_F = sqrt((kx_m.^2*xi_T_f     + kz_m.^2.*xi_L_f_m )./xi_T_f   );
    
    %% Iterate through the point sources
    
    B_den = 2*xi_T_f*eta_F.*(sigma_V*eta_V+sigma_L_y*eta_L).*cosh(eta_F*h_F) ...
        + 2*(sigma_V*eta_V*sigma_L_y.*eta_L+xi_T_f^2*eta_F.^2).*sinh(eta_F*h_F);
    
    for i = 1:length(I_M)
        %% Define point source stimulation in the time Fourier domain
        
        % Monophasic
        % I_hat = I_M(i)/sqrt(2*pi)*I_D(i)*exp(-1i*I_D(i)*w_m/2) ...
        %     .*sinc(I_D(i)*w_m/2/pi);
        
        % Biphasic
        % Apply Lanczos sigma factor to reduce Gibbs phenomenon
        sig = sinc(T(j)/d_t*pi*0.5/length(T));
        sig_m = permute(repmat(sig,length(Xp),1,length(Zp)),[1 3 2]);
        I_hat = 1i*2*I_D(i)*I_M(i)/sqrt(2*pi)*exp(-1i*I_D(i)*w_m) ...
            .*sinc(I_D(i)*w_m/2/pi).*sin(I_D(i)*w_m/2).*sig_m;
        
        % Sinusoidal
        %     I = I_M(i)*sqrt(2)*sin(pi*T/I_D(i)).*heaviside(T).*heaviside(-T+2*I_D(i));
        %     I_hat_vec = d_t/sqrt(2*pi)*ifftshift(fft(fftshift(I)));
        %     I_hat = permute(repmat(I_hat_vec,length(Xp),1,length(Zp)),[1 3 2]);
        
        %% Calculate extracellular voltage contributed by this source
        % in the Fourier domain
        
        % Define B.C. coefficients
        
        if Ri(i) == 0
            m = I_hat./(2*pi*sigma_V);
        elseif Ri(i) > 0
            m = I_hat.*besselj(1,sqrt(kx_m.^2 + kz_m.^2).*Ri(i))./...
                (pi*sigma_V*Ri(i)*sqrt(kx_m.^2 + kz_m.^2));
        else
            error('Electrode radius should be >= 0');
        end
        
        B_1 = m.*sigma_V./(B_den)...
            .*(xi_T_f*eta_F - sigma_L_y*eta_L).*exp(-eta_F*h_F-eta_V*Yi(i));
        B_2 = m.*sigma_V./(B_den)...
            .*(xi_T_f*eta_F + sigma_L_y*eta_L).*exp(eta_F*h_F-eta_V*Yi(i));
        
        e_shift = exp(-1i*kx_m*Xi(i)-1i*kz_m*Zi(i));
        
        % Extracellular voltage in Fourier domain for each layer
        % Vitreous
        if Ya > 0 && Ya <= Yi(i)
            A_1 = B_1 + B_2 - m./(2*eta_V).*exp(-eta_V*Yi(i));
            Ve_L_f(:,:,j) = Ve_L_f(:,:,j) + (A_1.*exp(-eta_V*Ya) + ...
                m./(2*eta_V).*exp(eta_V*(Ya-Yi(i)))).*e_shift;
            Ve_T_Y_f(:,:,j) = Ve_T_Y_f(:,:,j) + (-A_1.*eta_V.*exp(-eta_V*Ya) + ...
                m./2.*exp(eta_V*(Ya-Yi(i)))).*e_shift;
            % Nerve Fiber Layer
        elseif Ya <= 0 && Ya >= -h_F
            Ve_L_f(:,:,j) = Ve_L_f(:,:,j) + (B_1.*exp(-eta_F*Ya) + B_2.*exp(eta_F*Ya)).*e_shift;
            Ve_T_Y_f(:,:,j) = Ve_T_Y_f(:,:,j) + (-B_1.*eta_F.*exp(-eta_F*Ya) + ...
                B_2.*eta_F.*exp(eta_F*Ya)).*e_shift;
            % Lower Layers
        elseif Ya < -h_F
            C_2 = (B_1.*exp(eta_F*h_F) + B_2.*exp(-eta_F*h_F))./exp(-eta_L*h_F);
            Ve_L_f(:,:,j) = Ve_L_f(:,:,j) + C_2.*exp(eta_L*Ya).*e_shift;
            Ve_T_Y_f(:,:,j) = Ve_T_Y_f(:,:,j)+C_2.*eta_L.*exp(eta_L*Ya).*e_shift;
        else
            error('Ya (depth of plane of analysis) should be between -h_F and Yi');
        end
    end
    
    neur_eq(:,:,j) = -kzp_m.^2.*L_V_m.^2./(1+kzp_m.^2.*L_V_m.^2);
    
    %% Calculate membrane voltage in the time-space domain
    
    Ve_T_X_f(:,:,j) = b/2*1i*kxp_m.*Ve_L_f(:,:,j);
    Ve_T_Y_f(:,:,j) = b/2*Ve_T_Y_f(:,:,j);
    Ve_T_Z_f(:,:,j) = b/2*1i*kzp_m.*Ve_L_f(:,:,j);
    
end

end