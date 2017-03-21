%%SIMULATIONPARAMETERS
% SIMULATIONPARAMETERS defines the simulation geometry and parameters for a
% number of different simulations
%
%%%%%%%%%%%%%%%%%%%%%%%%% Created by: Tim Esler, 2017 %%%%%%%%%%%%%%%%%%%%%%%%%

if simulation == 1
    
    %% One Electrode Simulation
    
    % Define neurite rotations in xz plane
    d_rot = pi/18;
    rot = [0 0:d_rot:pi/2+pi/10000];
    
    % Define source locations
    Xi = 0e-6;
    Yi = 200e-6;
    Zi = 0e-6;
    
    % Define electrode radius
    Ri = 50e-6;
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M = -50e-6;
    I_D = 200e-6;
    
    % Should we scale the output so the AOP just hits threshold, or just use the
    % I_m value supplied?
    scaleOut = true;
    
    % Folder name
    folderName = 'Figures';
    
    % Plot filename
    figName = 'RotNeurite_OneElectrode';
    
elseif simulation == 2
    
    %% Two Electrode Simulation
    
    % Define neurite rotations in xz plane
    d_rot = pi/18;
    rot = [0 0:d_rot:pi/2+pi/10000];
    
    % Define source locations
    Xi = [0e-6 0e-6];
    Yi = [200e-6 200e-6];
    Zi = [-100e-6 100e-6];
    
    % Define electrode radius
    Ri = [50e-6 50e-6];
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M = [-50e-6 -50e-6];
    I_D = [200e-6 200e-6];
    
    % Should we scale the output so the AOP just hits threshold, or just use the
    % I_m value supplied?
    scaleOut = true;
    
    % Folder name
    folderName = 'Figures';
    
    % Plot filename    
    figName = 'RotNeurite_TwoElectrode';
    
elseif simulation == 4
    
    %% Four Electrode Simulation
    
    % Define neurite rotations in xz plane
    d_rot = pi/18;
    rot = [0 0:d_rot:pi/2+pi/10000];
    
    % Define source locations
    Xi = [0e-6 0e-6 0e-6 0e-6];
    Yi = [200e-6 200e-6 200e-6 200e-6];
    Zi = [-300e-6 -100e-6 100e-6 300e-6];
    
    % Define electrode radius
    Ri = [50e-6 50e-6 50e-6 50e-6];
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M = [-50e-6 -50e-6 -50e-6 -50e-6];
    I_D = [200e-6 200e-6 200e-6 200e-6];
    
    % Should we scale the output so the AOP just hits threshold, or just use the
    % I_m value supplied?
    scaleOut = true;
    
    % Folder name
    folderName = 'Figures';
    
    % Plot filename
    figName = 'RotNeurite_FourElectrode';
    
end