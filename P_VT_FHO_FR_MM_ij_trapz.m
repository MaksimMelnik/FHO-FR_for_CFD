function pVT=P_VT_FHO_FR_MM_ij_trapz(M1, Coll, i1, f1, E)
% FHO-FR probability averaged over collision parameters depending on 
% initial and final vibrational levels, collison molecules and collision 
% energy. The averaging uses hybrid trapz and integral2 integration. 
% Diatom-diatom collision. M1 is the first molecule, Coll is the 
% collision variable, i1 is the initial vibrational state of M1, 
% f1 is the final vibrational state of M1, E is the collision energy (J).
% src: eq.15, I. V. Adamovich et al., J. Chem. Phys. 109, 7711-7724 (1998) 
% https://doi.org/10.1063/1.477417
% 03.11.2022 Maksim Melnik

eps2max=@(eps1) 1-eps1;     % the integration limit function for eps2
error_v=1e-2;
pVT=integral2(@(eps1, eps2) ... integration over eps via integral2
        integrant_trapz(M1, Coll, i1, f1, E, eps1, eps2), ...
                    0,1, 0,eps2max, 'RelTol',error_v,'AbsTol',error_v);
pVT=pVT/pi^4;
end

function out=integrant_trapz(M1, Coll, i1, f1, E, eps1, eps2)
% integration over thetas, phis and y via trapz
maxdiv=10;  % 3 - 0.2с, 4 - wrong, 17 - 60с
ymax=1; theta_v1max=pi; theta_v2max=pi; phi1max=pi; phi2max=pi;
y=0:ymax/maxdiv:ymax;
yrs=reshape(y, 1, 1, []);                       % 3
theta_v1=0:theta_v1max/maxdiv:theta_v1max;
theta_v1rs=reshape(theta_v1, 1, 1, 1, []);      % 4
theta_v2=0:theta_v2max/maxdiv:theta_v2max;
theta_v2rs=reshape(theta_v2, 1, 1, 1, 1, []);   % 5
phi1=0:phi1max/maxdiv:phi1max;
phi1rs=reshape(phi1, 1, 1, 1, 1, 1, []);        % 6
phi2=0:phi2max/maxdiv:phi2max;
phi2rs=reshape(phi2, 1, 1, 1, 1, 1, 1, []);     % 7
out=P_VT_FHO_FR_MM_ij(M1, Coll, i1, f1, ...
            E, eps1, yrs, theta_v1rs, phi1rs, eps2, theta_v2rs, phi2rs);
out=trapz(phi2, trapz(phi1, trapz(theta_v2, trapz(theta_v1, ...
                                        trapz(y, out, 3), 4), 5), 6), 7);
end