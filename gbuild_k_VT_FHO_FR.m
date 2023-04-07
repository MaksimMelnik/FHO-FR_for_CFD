% The example of usage and plotting of k_VT_FHO_FR vs data from figs 12-15
% I. V. Adamovich et al., J. Chem. Phys. 109, 7711-7724
% (1998) https://doi.org/10.1063/1.477417
% 31.03.2023
% Maksim Melnik
% 1. k_VT calculation for comparison with Adamovich 1998 FR fig 12
% 2. plotting k_VT from Adamovich 1998 FR fig 12 vs Maksim's calculation
%% 1. k_VT calculation for comparison with Adamovich 1998 FR fig 12
T_num_parts=10;
Tmin=200;
Tmax=50200;
T_arr=10.^(log10(Tmin):(log10(Tmax)-log10(Tmin))/T_num_parts:log10(Tmax));
kO2Ar_10=T_arr*0;
kO2Ar_109=T_arr*0;
kO2Ar_2019=T_arr*0;
kO2Ar_3029=T_arr*0;
for ind=1:length(T_arr)
    kO2Ar_10(ind)=k_VT_FHO_FR(T_arr(ind), O2, Ar, 1, 0);
    kO2Ar_109(ind)=k_VT_FHO_FR(T_arr(ind), O2, Ar, 10, 9);
    kO2Ar_2019(ind)=k_VT_FHO_FR(T_arr(ind), O2, Ar, 20, 19);
    kO2Ar_3029(ind)=k_VT_FHO_FR(T_arr(ind), O2, Ar, 30, 29);
end
kO2Ar_10=kO2Ar_10*1e6;
kO2Ar_109=kO2Ar_109*1e6;
kO2Ar_2019=kO2Ar_2019*1e6;
kO2Ar_3029=kO2Ar_3029*1e6;
%% 2. plotting k_VT from Adamovich 1998 FR fig 12 vs Maksim's calculation
load data_Adamovich1998_FR_fig12
loglog(k_VT_FHO_FR_O2Ar_10_A98_FR(:, 1), ...
                                k_VT_FHO_FR_O2Ar_10_A98_FR(:, 2), 'k', ...
    k_VT_FHO_FR_O2Ar_109_A98_FR(:, 1), ...
                            k_VT_FHO_FR_O2Ar_109_A98_FR(:, 2), 'k', ...
	k_VT_FHO_FR_O2Ar_2019_A98_FR(:, 1), ...
                            k_VT_FHO_FR_O2Ar_2019_A98_FR(:, 2), 'k', ...
	k_VT_FHO_FR_O2Ar_3029_A98_FR(:, 1), ...
                                k_VT_FHO_FR_O2Ar_3029_A98_FR(:, 2), ...
                                                    'k', 'linewidth', 1.5)
hold on
loglog(T_arr, kO2Ar_10, ':', T_arr, kO2Ar_109, ':', ...
        T_arr, kO2Ar_2019, ':', T_arr, kO2Ar_3029, ':', 'linewidth', 1.5)

legend('ADIAV {\it k}_{10}', 'ADIAV {\it k}_{109}', ...
        'ADIAV {\it k}_{2019}', 'ADIAV {\it k}_{3029}', ...
        "Maksim's {\it k}_{10}", "Maksim's {\it k}_{109}", ...
        "Maksim's {\it k}_{2019}", "Maksim's {\it k}_{3029}", ...
                                                    'location', 'best')
xlim([200 50000])
ylim([1e-20 1e-9])
title('Adamovich 1998 FR fig. 12, {\it k}_{VT} O_2-Ar')
ylabel('{\it k}_{VT}, cm^3/s')
xlabel('{\it T}, K')