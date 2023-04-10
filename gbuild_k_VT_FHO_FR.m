% The example of usage and plotting of k_VT_FHO_FR vs data from figs 12-15
% I. V. Adamovich et al., J. Chem. Phys. 109, 7711-7724
% (1998) https://doi.org/10.1063/1.477417
% 31.03.2023
% Maksim Melnik
% 1. k_VT calculation for comparison with Adamovich 1998 FR fig 12
% 2. plotting k_VT from Adamovich 1998 FR fig 12 vs Maksim's calculation
%% 1. k_VT calculation for comparison with Adamovich 1998 FR fig 12
load particles O2 Ar
T_num_parts=10;
Tmin=200;
Tmax=50200;
T_arr=10.^(log10(Tmin):(log10(Tmax)-log10(Tmin))/T_num_parts:log10(Tmax));
kO2Ar_10=T_arr*0;
kO2Ar_109=T_arr*0;
kO2Ar_2019=T_arr*0;
kO2Ar_3029=T_arr*0;
tic
for ind=1:length(T_arr)
    kO2Ar_10(ind)=k_VT_FHO_FR(T_arr(ind), O2, Ar, 1, 0);
    kO2Ar_109(ind)=k_VT_FHO_FR(T_arr(ind), O2, Ar, 10, 9);
    kO2Ar_2019(ind)=k_VT_FHO_FR(T_arr(ind), O2, Ar, 20, 19);
    kO2Ar_3029(ind)=k_VT_FHO_FR(T_arr(ind), O2, Ar, 30, 29);
end
toc
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
%% optimization research (currently failed)
k = 1.380649e-23;                   % Boltzmann constant, J/K
h = 6.626070041e-34;                % Plank constant, J*sec
c = 299792458;                      % speed of light, m/sec
T=15000;
i1=1;
f1=i1-1;
T_cm1=T*k/h/c/100;                  % in cm-1
error_val=1e-3;
tic
Emin=1e0; Emax=1e7;     % standart min and max E value on the plot
fr_num=250;
E_arr=10.^(log10(Emin):(log10(Emax)-log10(Emin))/fr_num:log10(Emax));%...
%                                                                 *100*h*c;
temp=(E_arr./T_cm1).^2 .* P_VT_FHO_FR_avg_angles(O2, Ar, ...
            i1, f1, E_arr*h*c*100, 't') .* exp(-E_arr./T_cm1)/T_cm1;
int_t=trapz(E_arr, temp);
toc
tic
int=integral(@(E) (E./T_cm1).^2 .* P_VT_FHO_FR_avg_angles(O2, Ar, ...
            i1, f1, E*h*c*100, 't') .* exp(-E./T_cm1)/T_cm1, 0, Inf, ...
                            'RelTol', error_val, 'AbsTol', error_val);
toc
disp('delta, %')
(int_t-int)/int
%%
% какие параметры могут меняться? 
% могут меняться i1, f1, T, M1, M2
% Хочу придумать верхний и нижний пределы Emin, Emax для интегрирования 
% трапециями
% на что влияет изменение M2:
%     вероятность растёт с массой M2
% Потом посмотреть, на что влияет изменение i1
%     нижний предел падает Emin с ростом i1
%     верхний предел не понятно, вроде тоже падает Emax
% потом изменение f1
%     с ростом i1-f1 вероятность падает, так что оптимально для di=1
% потом изменение T
%     с ростом T растёт Emax
%     
% ну и в конце зависимость от изначальной молекулы M1