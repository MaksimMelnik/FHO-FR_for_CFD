% The example of usage and plotting of P_VT_FHO_FR_MM_ij_trapz(T) vs data
% from fig. 1 (right) S. F. Gimelshein et al., J. Thermophys. Heat Transf.
% 2018 32:4, 882-891 https://doi.org/10.2514/1.T5228
% and vs Maksim's data obtained using pure integration with integral8 and
% without trapz.
% 06.09.2022 Maksim Melnik

tic
load par_data
k = 1.380649e-23;       % Boltzmann constant, J/K
h = 6.626070041e-34;    % Plank constant, J*sec
c = 299792458;          % speed of light, m/sec
Emin=4e3; Emax=2e5;     % min and max E values on the Gimelshein2017
                        % paper's plot, cm^-1
i1=1;                   % initial state, P_VT for the 1->0 transition
fr=[0 0.0005 0.001 0.002 0.003 0.005 ... just a fragmentation of E array
0.006 0.008 0.01 0.02 0.04 0.05 0.055 0.06 0.07 0.08 0.1 0.2 0.5 0.8 1];
fr=[0 0.0005 0.001 0.005 0.01 0.05 0.1 0.5 1]; % a light version of a net
E_arr=(fr*Emax+(1-fr)*Emin)*100*h*c;    % fragmentation
M1=N2;                                  % firs molecule
p_VT_10_tr=zeros(1, length(E_arr));
%% call of P_VT using hybrid trapz integration
for ind=1:length(E_arr)
p_VT_10_tr(ind)=...
            P_VT_FHO_FR_MM_ij_trapz(M1, Coll_N2_N2, i1, i1-1, E_arr(ind));
end
%%
E_cm=E_arr/h/c/100;             % E in cm-1
load data_Gimelshein2017_f1r    % loading of data from Gimelshein plot
loglog(p_VT_FHO_FR_N2N2_10_G17(:,1), p_VT_FHO_FR_N2N2_10_G17(:, 2), ...
            'linewidth', 1.5)   % plotting Gimelshein plot data
hold on
    % data obtained by Maksim via integral8
p_VT=[1.0042529537425651E-8 0.008 1.5393];
E_cm_int8=[4980 23600 1.0200e+05];     % E array for p_VT ^ for plotting v
loglog(E_cm_int8, p_VT, 'o', 'markerfacecolor', 'r', 'markersize', 10)
loglog(E_cm(1:ind), p_VT_10_tr, 'linewidth', 1.5) % plotting trapz results
legend('FHO-FR Gimelshein et al. 2017', "Maksim's integral8 code", ...
                                "Maksim's trapz code", 'location', 'best')
xlim([1e3 2e5])
ylim([1e-11 1e1])
title('N_2-N_2 p_{VT} FHO-FR 1->0')
%%
toc