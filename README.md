# FHO-FR_for_CFD
The repository contains various functions for the Forced Harmonic Oscillator Free-Rotation model [1, 2]. The functions are developed for the implementation into MATLAB CFD codes.

## Functions descriptions
- `P_VT_FHO_FR_MA_ij.m` returns P_VT for molecule-atom collision, depending on *M<sub>1</sub>*, *Coll*, *i<sub>1</sub>*, *f<sub>1</sub>*, *E*, $\varepsilon_1$, *y*, $\vartheta_1$, $\phi_1$ [1].
- `P_VT_FHO_FR_avg_angles.m` returns P_VT molecule-atom collision, depending on *M<sub>1</sub>*, *Coll*, *i<sub>1</sub>*, *f<sub>1</sub>*, *E* and optional *method*. This fuction averages `P_VT_FHO_FR_MA_ij.m` over $\varepsilon_1$, *y*, $\vartheta_1$ and $\phi_1$ [1] using numeric multiple intergation `integralN` (for *method*='q' or by default) or using `trapz` integration (for *method*='t'). NOTE: before using be sure you downloaded [integralN](https://www.mathworks.com/matlabcentral/fileexchange/47919-integraln-m).
- `P_VT_FHO_FR_MM_ij.m` returns P_VT for molecule-molecule collision, depppending on $M_1$, *Coll*, *i<sub>1</sub>*, *f<sub>1</sub>*, *E*, $\varepsilon_1$, *y*, $\vartheta_1$, $\phi_1$, $\varepsilon_2$, $\vartheta_2$, $\phi_2$.
- `P_VT_FHO_FR_MM_ij_trapz.m` returns P_VT molecule-molecule collision, depending on *M<sub>1</sub>*, *Coll*, *i<sub>1</sub>*, *f<sub>1</sub>* and *E*. This fuction averages `P_VT_FHO_FR_MM_ij.m` over $\varepsilon_1$, *y*, $\vartheta_1$, $\phi_1$, $\varepsilon_2$, $\vartheta_2$ and $\phi_2$ using hybrid integration method with `integral2` for $\varepsilon_1$, $\varepsilon_2$ and `trapz` integration for other parameters.
- `particles_data_ini.m` initializes particles and collisions variables in the format prepared for present functions.  
- `par_data.mat` containes some particles and collisions variables.  
- `levels_e_ex.m` calculates vibrational energies of states.  
##### Examples
- `gbuild_p_VT_FHO_FR_avg_angles.m` script provides an example of `P_VT_FHO_FR_avg_angles.m` usage and plotting data vs [1].
- `gbuild_p_VT_FHO_FR_avg_angles_MM_trapz.m` script provides an example of `P_VT_FHO_FR_MM_ij_trapz.m` usage and plotting data vs [4].

Here $M_1$ is the molecule under consideration variable, *Coll* is the collision variable, *i<sub>1</sub>* is the initial vibrational level of *M1*, *f<sub>1</sub>* is the final vibrational level of *M1*, *E* is the collision energy, $\varepsilon_1$ is the fruction of vibrational energy for *M1*, *y* is the collision parameter, $\vartheta_1$ is the rotation angle of *M1*, $\phi_1$ is the angle between the plane of rotation and the radius vector for *M1*.

## to do
- add fast trapz integral for P_VT MA
- restruct repository


## References
[1]  [Adamovich I V and Rich J W 1998 *J. Chem. 
Phys.* **109** 18](https://doi.org/10.1063/1.477417).

[2]  [Adamovich I V 2001 *AIAA J.* **39** 10](https://doi.org/10.2514/2.1181).

[3]  [Mike Hosea (2022). integralN.m, MATLAB Central File Exchange.](https://www.mathworks.com/matlabcentral/fileexchange/47919-integraln-m)

[4]  [Gimelshein S F *et al* 2017 *J. Thermophys. Heat Transf.* **32** 4](https://doi.org/10.2514/1.T5228).