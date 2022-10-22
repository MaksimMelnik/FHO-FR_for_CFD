# FHO-FR_for_CFD
The repository contains various functions for the Forced Harmonic Oscillator Free-Rotation model [1]. The functions are developed for the implementation into MATLAB CFD codes.

## Functions descriptions
- `P_VT_FHO_FR_MA_ij.m` returns P_VT for molecule-atom collision, depending on *M1*, *Coll*, *i<sub>1</sub>*, *f<sub>1</sub>*, *E*, $\varepsilon_1$, y, $\vartheta_1$, $\phi_1$.  
- `P_VT_FHO_FR_avg_angles.m` returns P_VT molecule-atom collision, depending on *M1*, *Coll*, *i<sub>1</sub>*, *f<sub>1</sub>* and *E*. This fuction averages `P_VT_FHO_FR_MA_ij.m` over $\varepsilon_1$, y, $\vartheta_1$ and $\phi_1$ using numeric multiple intergation `integralN`. NOTE: before using be sure you downloaded [integralN](https://www.mathworks.com/matlabcentral/fileexchange/47919-integraln-m).
- `particles_data_ini.m` initializes particles and collisions variables in the format prepared for present functions.  
- `par_data.mat` containes some particles and collisions variables.  
- `levels_e_ex.m` calculates vibrational energies of states.  
##### Examples
- `gbuild_p_VT_FHO_FR_avg_angles.m` script provides an example of `P_VT_FHO_FR_avg_angles.m` usage and plotting data vs [1].

Here *M1* is the molecule under consideration variable, *Coll* is the collision variable, *i<sub>1</sub>* is the initial vibrational level of *M1*, *f<sub>1</sub>* is the final vibrational level of *M1*, *E* is the collision energy, $\varepsilon_1$ is the fruction of vibrational energy for *M1*, *y* is the collision parameter, $\vartheta_1$ is the rotation angle of *M1*, $\phi_1$ is the angle between the plane of rotation and the radius vector for *M1*.

## to do
- add P_VT for MM
- add trapz integral averaging of P_VT MM over $\varepsilon_1$, y, $\vartheta_1$, $\phi_1$, $\varepsilon_2$, $\vartheta_2$, $\phi_2$
- add usage examples for P_VT MM
- add fast trapz integral for P_VT MA
- restruct repository


## References
[1]  [I. V. Adamovich and J. W. Rich, J. Chem. Phys., Vol. 109, No. 18, 8 November 1998](https://doi.org/10.1063/1.477417).
