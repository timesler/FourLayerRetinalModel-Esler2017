# FourLayerRetinalModel-Esler2017
Code to simulate epiretinal stimulation of the retina using a four-layer model.

## CellComp4Layer_Ve_f_Plane.m

*CellComp4Layer_Ve_f_Plane* returns the electric field in the Fourier domain due to epiretinal stimulation with point or disk electrodes in a plane in the x-z dimensions (where y is in the direction normal to the surface of the retina, toward the electrode). The retina is modelled using a 4-layer description, including Insulator, Vitreous, Nerve Fibre Layer, and the Ganglion Cell Layer (+ others). The insulator layer is modelled using a no-current boundary condition within the vitreous (see diagram below for geometry).
 
All units are S.I.
  
        _____________________________________________________
        1.Insulator
        _____________________________________________________
                    [|||]        [|||]        [|||]
                     Electrodes (at layer boundary)
 
        2.Vitreous
        _____________________________________________________
        3.Nerve fibre layer
        _____________________________________________________
        4.'Other' cell layer (incl. GCL)
        _____________________________________________________


The extracellular potential for the NFL is calculated using a modified version of the self-consistent, linear, sub-threshold model presented in:

* B. Tahayori, H. Meffin, E.N. Sergeev, I.M.Y. Mareels, A.N. Burkitt, and D.N. Grayden (2014), "Modelling extracellular electrical stimulation: IV. Effect of the cellular composition of neural tissue on its spatio-temporal filtering properties", J. Neural Eng. 11.

This script was used to conduct the analysis presented in:

\<Enter correct reference here\>

### Inputs

| Input | Description |
|---|---|
| Xi, Yi and Zi | x-, y-, and z- coordinates of the point source electrodes. (m) |
| I_M, I_D | Stimulation amplitude and duration for each electrode. (A, s) |
| z_max, x_max, t_max | z-, x- and time-extent of the simulation. (m, m, s) |
| d_z, d_x, d_t | z , x and t step sizes. (m, m, s) |
| h_F | Nerve fibre layer thickness (m) |
| Ya | Depth below the retina surface at which we want to calculate extracellular potential. (m) |
| rot | Coordinate rotation for a fibre which is rotated w.r.t. the fibre bindle (only set to a nonzero value if the output will be used to calculate the membrane potential of a rotated fiber). (rad) |
| Ri | Radius of disk for each electrode (set to 0 for point source). (m) |

### Outputs

| Output | Description |
| --- | ---|
| Ve_L_f | The calculated longitudinal extracellular potential for a full plane in the Fourier domain. This is equal to the extracellular potential itself. (V) |
| Ve_T_#_f | The calculated transverse extracellular potential in the #-direction for a full plane in the Fourier domain. (V) |
| neur_eq | The transformation matrix for converting the longitudinal extracellular potential into longitudinal membrane potential. Multiplying Ve_L_f by this pointwise will give the membrane potential in the x, z, and t Fourier domain. |

### Example usage

    Xi = 0e-6;
    Yi = 200e-6;
    Zi = 0e-6;
    I_M = -50e-6;
    I_D = 200e-6;
    x_max = 2000e-6;
    z_max = 2000e-6;
    t_max = 600e-6;
    d_x = 10e-6;
    d_z = 10e-6;
    d_t = 4e-6;
    h_F = 100e-6;
    Ri = 50e-6;
    Ya = -10e-6;

    [Ve_L_f, Ve_T_X_f, Ve_T_Y_f, Ve_T_Z_f, n_e] = CellComp4Layer_Ve_f_Plane(...
        Xi, Yi, Zi, I_M, I_D, ...                % Electrode parameters
        x_max, z_max, t_max, d_x, d_z, d_t, ...  % Spatial sampling
        h_F, Ya, rot, Ri);                       % Simulation geometry

A further example can be found in *AnyElecConfig_RotatedNeurite_4L*.

## CellComp3Layer_Ve_f_Plane.m

NOTE: Identical in implementation to CellComp3Layer_Ve_f_Plane but with a 3-layer model instead of 4, to simulate stimulation where no insulator is used.

*CellComp3Layer_Ve_f_Plane* returns the electric field in the Fourier domain due to epiretinal stimulation with point or disk electrodes in a plane in the x-z dimensions (where y is in the direction normal to the surface of the retina, toward the electrode). The retina is modelled using a 3-layer description, including Vitreous, Nerve Fibre Layer, and the Ganglion Cell Layer (+ others) (see diagram below for geometry).
 
All units are S.I.

                    [|||]        [|||]        [|||]
                     Electrodes (at layer boundary)
 
        2.Vitreous
        _____________________________________________________
        3.Nerve fibre layer
        _____________________________________________________
        4.'Other' cell layer (incl. GCL)
        _____________________________________________________


The extracellular potential for the NFL is calculated using a modified version of the self-consistent, linear, sub-threshold model presented in:

* B. Tahayori, H. Meffin, E.N. Sergeev, I.M.Y. Mareels, A.N. Burkitt, and D.N. Grayden (2014), "Modelling extracellular electrical stimulation: IV. Effect of the cellular composition of neural tissue on its spatio-temporal filtering properties", J. Neural Eng. 11.

This script was used to conduct the analysis presented in:

\<Enter correct reference here\>

### Inputs

| Input | Description |
|---|---|
| Xi, Yi and Zi | x-, y-, and z- coordinates of the point source electrodes. (m) |
| I_M, I_D | Stimulation amplitude and duration for each electrode. (A, s) |
| z_max, x_max, t_max | z-, x- and time-extent of the simulation. (m, m, s) |
| d_z, d_x, d_t | z , x and t step sizes. (m, m, s) |
| h_F | Nerve fibre layer thickness (m) |
| Ya | Depth below the retina surface at which we want to calculate extracellular potential. (m) |
| rot | Coordinate rotation for a fibre which is rotated w.r.t. the fibre bindle (only set to a nonzero value if the output will be used to calculate the membrane potential of a rotated fiber). (rad) |
| Ri | Radius of disk for each electrode (set to 0 for point source). (m) |

### Outputs

| Output | Description |
| --- | ---|
| Ve_L_f | The calculated longitudinal extracellular potential for a full plane in the Fourier domain. This is equal to the extracellular potential itself. (V) |
| Ve_T_#_f | The calculated transverse extracellular potential in the #-direction for a full plane in the Fourier domain. (V) |
| neur_eq | The transformation matrix for converting the longitudinal extracellular potential into longitudinal membrane potential. Multiplying Ve_L_f by this pointwise will give the membrane potential in the x, z, and t Fourier domain. |

### Example usage

    Xi = 0e-6;
    Yi = 200e-6;
    Zi = 0e-6;
    I_M = -50e-6;
    I_D = 200e-6;
    x_max = 2000e-6;
    z_max = 2000e-6;
    t_max = 600e-6;
    d_x = 10e-6;
    d_z = 10e-6;
    d_t = 4e-6;
    h_F = 100e-6;
    Ri = 50e-6;
    Ya = -10e-6;

    [Ve_L_f, Ve_T_X_f, Ve_T_Y_f, Ve_T_Z_f, n_e] = CellComp3Layer_Ve_f_Plane(...
        Xi, Yi, Zi, I_M, I_D, ...                % Electrode parameters
        x_max, z_max, t_max, d_x, d_z, d_t, ...  % Spatial sampling
        h_F, Ya, rot, Ri);                       % Simulation geometry

## AnyElecConfig_RotatedNeurite.m

*AnyElecConfig_RotatedNeurite* simulates the epiretinal stimulation of retinal ganglion cell axons in the NFL (axons of passage) and in the ganglion cell layer (axon initial segments).

This script defines the:
  * Simulation numerics (i.e. spatial and temporal extent, sampling rates)
  * Geometry of the stimulating electrode(s)
  * Location and orientation of the simulated axons
  * Stimulation waveform (i.e. amplitude and phase duration for a biphasic pulse)

These parameters are then passed to the simulation function *CellComp4Layer_Ve_f_Plane* which performs the specified simulations (type "help CellComp4Layer_Ve_f_Plane" for more inforation). This script can also be used to call the 3-layer version, *CellComp3Layer_Ve_f_Plane*. *CellComp4Layer_Ve_f_Plane* and *CellComp3Layer_Ve_f_Plane* return the extracellular potential in the specified plane in the Fourier domain. In order to calculate the longitudinal and transverse components of the membrane potential, this script applies the supplied neurite equation (neur_eq) before calculating the inverse Fourier transform.
