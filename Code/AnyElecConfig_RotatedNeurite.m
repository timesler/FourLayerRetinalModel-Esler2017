%%ANYELECCONFIG_ROTATEDNEURITE_4L
% ANYELECCONFIG_ROTATEDNEURITE_4L Simulates the epiretinal stimulation of 
% retinal ganglion cell axons in the NFL (axons of passage) and in the
% ganglion cell layer (axon initial segments).
%
% This script defines the:
%   * Simulation numerics (i.e. spatial and temporal extent, sampling rates)
%   * Geometry of the stimulating electrode(s)
%   * Location and orientation of the simulated axons
%   * Stimulation waveform (i.e. amplitude and phase duration for a
%   biphasic pulse)
%
% These parameters are then passed to the simulation function
% CellComp4Layer_Ve_f_Plane which performs the specified simulations (see
% "help CellComp4Layer_Ve_f_Plane" for more inforation).
% CellComp4Layer_Ve_f_Plane returns the extracellular potential in the
% specified plane in the Fourier domain. In order to calculate the longitudinal
% and transverse components of the membrane potential, this script applies
% the supplied neurite equation (neur_eq) before calculating the inverse
% Fourier transform.
%
%%%%%%%%%%%%%%%%%%%%%%%%% Created by: Tim Esler, 2017 %%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters

% Define simulation size and step sizes
x_max = 2000e-6;
z_max = 2000e-6;
t_max = 600e-6;
d_x = 10e-6;
d_z = 10e-6;
d_t = 4e-6;

Z = -z_max:d_z:z_max;
T = -t_max:d_t:t_max;
X = -x_max:d_x:x_max;
X_Centre = find(X == 0);

simulation = 4; % 1 = one electrode stimulation
                % 2 = two electrode stimulation 
                % 4 = four electrode stimulation

SimulationParameters % this script defines appropriate simulation parameters

% Define membrane thresholds
Th = [12.09e-3 6.30e-3*ones(1,length(rot)-1)];

%% Calculate the longitudinal and transverse components separately

Vm_L_rot = zeros(length(rot),length(Z),length(T));
Vm_T_rot = zeros(length(rot),length(Z),length(T));

h = waitbar(0,'Calculating Rotated Membrane Potential');


for i = 1:length(rot)
    
    [Ve_L_f, Ve_T_X_f, Ve_T_Y_f, ~,n_e] = CellComp4Layer_Ve_f_Plane(...
        Xi, Yi, Zi, I_M, I_D, ...                   % Electrode parameters
        x_max, z_max, t_max, d_x, d_z, d_t, ...     % Spatial sampling
        h_F, Ya(i), rot(i), Ri ...
        );
    
    % Get Vm_L
    Vm_L_f = n_e.*Ve_L_f;
    Vm_L = real(fftshift(ifftn(ifftshift((2*pi)^(3/2)/d_z/d_x/d_t*Vm_L_f))));
    % Get Vm_T
    Ve_T_X = real(fftshift(ifftn(ifftshift((2*pi)^(3/2)/d_z/d_x/d_t*Ve_T_X_f))));
    Ve_T_Y = real(fftshift(ifftn(ifftshift((2*pi)^(3/2)/d_z/d_x/d_t*Ve_T_Y_f))));
    
    Vm_L_rot(i,:,:) = squeeze(Vm_L(X_Centre,:,:));
    Vm_T_rot(i,:,:) = 2*sqrt(squeeze(Ve_T_X(X_Centre,:,:)).^2 + ...
        squeeze(Ve_T_Y(X_Centre,:,:)).^2);
    
    waitbar(i/length(rot),h)
    
end

close(h)

%% Sum the longitudinal and tranverse membrane potentials

Vm_rot = Vm_L_rot + Vm_T_rot;

Vm_rot_AOP = Vm_rot(1,:,:);
Vm_L_rot_AOP = Vm_L_rot(1,:,:);
Vm_T_rot_AOP = Vm_T_rot(1,:,:);

Vm_rot = Vm_rot(2:end,:,:);
Vm_L_rot = Vm_L_rot(2:end,:,:);
Vm_T_rot = Vm_T_rot(2:end,:,:);

rot_AOP = rot(1);
Ya_AOP = Ya(1);
rot = rot(2:end);
Ya = Ya(2:end);

%% Plot geometry of stimulation

ColorSet = varycolor(length(rot)+1);
ColorSet = [ColorSet; [1 .5 0]];

P4 = figure;
ax = axes;
set(ax, 'ColorOrder', ColorSet);
hold on
for i = 1:length(rot)
    h1 = plot3((-z_max:d_z:z_max)*1e6.*cos(rot(i)), ...
        (-z_max:d_z:z_max)*1e6.*sin(rot(i)), ...
        ones(size(-z_max:d_z:z_max))*Ya(i)*1e6, ...
        '-d','LineWidth',2,'MarkerSize',1);
end
grid
h2 = plot3(Zi*1e6,Xi*1e6,Yi*1e6,'o','MarkerFaceColor','k','MarkerSize',5);
h1 = plot3((-z_max:d_z:z_max)*1e6.*cos(rot_AOP), ...
    (-z_max:d_z:z_max)*1e6.*sin(rot_AOP), ...
    ones(size(-z_max:d_z:z_max))*Ya_AOP*1e6, ...
    '-d','LineWidth',2,'MarkerSize',1);
title('Geometry of simulations')
xlabel('Z (\mum)')
ylabel('X (\mum)')
zlabel('Y (\mum)')
axis equal
xlim([-z_max z_max]*1e6)
ylim([-z_max z_max]*1e6)
zlim([-max(abs([Yi(:);Ya(:)])) max(abs([Yi(:);Ya(:)]))]*1e6)
view(-45,12)
grid minor
hold off

annotation(P4,'textbox',[0 0 1 0.05],'String',...
    sprintf('t_{max} = %s, d_t = %s, d_x = %s, d_z = %s, Xi = %s, Yi = %s, Zi = %s, I_M = %s, I_D = %s, Ya = %s, h_F = %s',...
    mat2str(t_max),mat2str(d_t),mat2str(d_x),mat2str(d_z),mat2str(Xi),mat2str(Yi),mat2str(Zi),mat2str(I_M),mat2str(I_D),mat2str(Ya),mat2str(h_F)),...
    'LineStyle','none')

set(P4, 'Colormap', ColorSet);
cb = colorbar('SouthOutside');
caxis([0 120])
set(cb,'Ticks',([rot rot(end)+d_rot rot(end)+2*d_rot]+d_rot/2)*180/pi,'TickLabels',[num2cell(rot*180/pi) 'Electrode' 'AOP'])
ylabel(cb,'Neurite rotation')

%% Plots of membrane potential over time at neurite centre

% Longitudinal
P2 = figure;
ax = subplot(3,1,1);
hold on
set(ax, 'ColorOrder', ColorSet);
plot(T*1e6,squeeze(Vm_L_rot_AOP(1,round((z_max)/d_z+1),:))*1e3)
for i = 1:length(rot)
    plot(T*1e6,squeeze(Vm_L_rot(i,round((z_max)/d_z+1),:))*1e3)
end
title({'Membrane activation versus time at neurite center for various neurite orientations';...
    ' ';'V_{m,L}'})
ylabel('V_{m,L}')
xlim([0 T(end)*1e6]);
hold off
box off

% Transverse
ax = subplot(3,1,2);
hold on
set(ax, 'ColorOrder', ColorSet);
plot(T*1e6,squeeze(Vm_T_rot_AOP(1,round((z_max)/d_z+1),:))*1e3)
for i = 1:length(rot)
    plot(T*1e6,squeeze(Vm_T_rot(i,round((z_max)/d_z+1),:))*1e3)
end
title('V_{m,T}')
ylabel('V_{m,T}')
xlim([0 T(end)*1e6]);
hold off
box off

% Combined
ax = subplot(3,1,3);
hold on
set(ax, 'ColorOrder', ColorSet);
plot(T*1e6,squeeze(Vm_rot_AOP(1,round((z_max)/d_z+1),:))*1e3)
for i = 1:length(rot)
    plot(T*1e6,squeeze(Vm_rot(i,round((z_max)/d_z+1),:))*1e3)
end
title('V_{m}')
xlabel('Time (\mus)')
ylabel('V_{m}')
xlim([0 T(end)*1e6]);
hold off
box off

set(P2, 'Colormap', ColorSet(1:end-1,:));
cb = colorbar('East');
caxis([0 110])
set(cb,'Ticks',([rot rot(end)+d_rot]+d_rot/2)*180/pi,'TickLabels',['AOP' num2cell(rot*180/pi)],'YAxisLocation','right')
ylabel(cb,'Neurite rotation')

%% Plots of membrane potential along each neurite at end of catodic phase

figSize = [5 2 8 6]*3;
fontSize = 6*3;

fig = figure('Units','centimeters','Position',figSize);
set(fig, 'PaperPosition',figSize)

if scaleOut
    scale = Th(1)/max(Vm_rot_AOP(:));
    Vm_rot_AOP = Vm_rot_AOP*scale;
    Vm_rot = Vm_rot*scale;
    disp(['Total stimulus current: ' num2str(sum(I_M)*scale*1e6) 'uA'])
end

ax = subplot(1,1,1);
set(ax, 'ColorOrder', ColorSet);
hold on
plot(Z*1e6,squeeze(Vm_rot_AOP(1,:,floor((t_max+I_D(1))/d_t)))*1e3,'LineWidth',1.5)
for i = 1:length(rot)
    plot(Z*1e6,squeeze(Vm_rot(i,:,floor((t_max+I_D(1))/d_t)))*1e3,'LineWidth',1.5)
end
plot(Z*1e6,Th(1)*ones(1,length(Z))*1e3,'--','Color',ColorSet(1,:));
plot(Z*1e6,Th(2)*ones(1,length(Z))*1e3,'--k');
xlabel('Neurite axis (\mum)')
ylabel('V_{m} (mV)')
xlim([-1000 1000])
ylim([-6 14])
set(ax,'XTick',[-1000 0 1000])
set(ax,'YTick',[-4 0 4 8 12])
set(gca,'fontsize',fontSize)
ax.Box = 'off';

outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

saveas(fig, [fileparts(pwd()) '/' folderName '/' figName],'fig')
saveas(fig, [fileparts(pwd()) '/' folderName '/' figName],'jpeg')