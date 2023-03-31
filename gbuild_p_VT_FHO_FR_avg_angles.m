% The example of usage and plotting of P_VT_FHO_FR_avg_angles(E) vs data
% from figs 8 and 9 I. V. Adamovich et al., J. Chem. Phys. 109, 7711-7724
% (1998) https://doi.org/10.1063/1.477417
% and vs data from fig. 1 (left) S. F. Gimelshein et al. J. Thermophys.
% Heat. Trans. 32 (4) (2017) https://doi.org/10.2514/1.T5228
% 1. Calculation of P_VT using integralN for comparison with Adamovich FR
%       1998 fig. 8.
% 2. Calculation of P_VT using trapz for comparison with Adamovich FR 1998
%       fig. 8.
% 3. Plotting the comparison of P_VT with the paper plot Adamovich FR 1998
%       fig. 8.
% 4. Calculation of P_VT using integralN for comparison with Adamovich FR
%       1998 fig. 9.
% 5. Calculation of P_VT using trapz for comparison with Adamovich FR 1998
%       fig. 9.
% 6. Plotting the comparison of P_VT with the paper plot Adamovich FR 1998
%       fig. 9.
% 7. Gimelshein 2017 fig. 1 left comparison, trapz
% 8. Gimelshein 2017 fig. 1 left comparison, plotting
% 30.08.2022 Maksim Melnik

tic
% clear fr E_arr E_cm %p4039 %p01

k = 1.380649e-23;       % Boltzmann constant, J/K
h = 6.626070041e-34;    % Plank constant, J*sec
c = 299792458;          % speed of light, m/sec
load particles N2       % loading particles
Emin=1e3; Emax=1e6;     % standart min and max E value on the plot
calc_integralN=false;   % use integralN method for calculation or 
                        %                                     only trapz?
%% 1. Adamovich fig 8 comparison (to 0th level), calculation
i1=[1 2 3 5];           % initial levels
f1=0;                   % final state
Emax=9e5;               % max E for fig 8
fr_num=50;
E_arr=10.^(log10(Emin):(log10(Emax)-log10(Emin))/fr_num:log10(Emax))...
                                                                *100*h*c;
% fr=[0 0.0005 0.001 0.002 0.003 0.005 ... the splitiing E interval array
%   0.006 0.008 0.01 0.02 0.04 0.05 0.055 0.06 0.07 0.08 0.1 0.2 0.5 0.8 1];
% fr=[0 0.0005 0.001 0.005 0.01 0.05 0.1 0.5 1];
% E_arr=(fr*Emax+(1-fr)*Emin)*100*h*c;    % splitting of E interval
M1=N2;                  % N2-N2a collision
M2=N2;                  % currently doesn't matter
disp(['Length of E is ' num2str(length(E_arr)) '.'])
disp('Collision of N2 with a fictitious N2a.')
if calc_integralN
 disp('IntegralN calculations started.')
 pif=zeros(length(i1), length(E_arr));
 for ind_i1=1:length(i1)
    disp(['i1=', num2str(i1(ind_i1)), ', f1=', num2str(f1)])
    for ind=1:length(E_arr)
        pif(ind_i1, ind)=P_VT_FHO_FR_avg_angles(M1, M2, ...
                                            i1(ind_i1), f1, E_arr(ind));
    end
 end
end
E_cm=E_arr/h/c/100;     % switching from J to cm-1
%% 2. Adamovich fig 8 comparison (to 0th level), calculation with trapz
pif_trapz=zeros(length(i1), length(E_arr));
disp('Trapz calculations started.')
for ind_i1=1:length(i1)
    disp(['i1=', num2str(i1(ind_i1)), ', f1=', num2str(f1)])
    for ind=1:length(E_arr)
        pif_trapz(ind_i1, ind)=P_VT_FHO_FR_avg_angles(M1, M2, ...
                                        i1(ind_i1), f1, E_arr(ind), 't');
    end
end
%% 3. Adamovich fig 8 comparison (to 0th level), plotting
load data_Adamovich98_2_figs8_9.mat
figure
if calc_integralN
 loglog(E_cm, pif(1,:), E_cm, pif(2,:), E_cm, pif(3,:), E_cm, pif(4,:), ...
                                                        'linewidth', 1.5)
 hold on
end
loglog(E_cm, pif_trapz(1,:), ':', E_cm, pif_trapz(2,:), ':',...
        E_cm, pif_trapz(3,:), ':', E_cm, pif_trapz(4,:), ':',...
                                                        'linewidth', 1.5)
if ~calc_integralN
    hold on
end
loglog(p_VT_FHO_FR_N2N2a_10_A98_2(:,1), ...
                                     p_VT_FHO_FR_N2N2a_10_A98_2(:,2), ...
    p_VT_FHO_FR_N2N2a_20_A98_2(:,1), p_VT_FHO_FR_N2N2a_20_A98_2(:,2), ...
    p_VT_FHO_FR_N2N2a_30_A98_2(:,1), p_VT_FHO_FR_N2N2a_30_A98_2(:,2), ...
    p_VT_FHO_FR_N2N2a_50_A98_2(:,1), p_VT_FHO_FR_N2N2a_50_A98_2(:,2), ...
                                    'color', [0,0,0], 'linewidth', 1.5)
legend_content=["Maksim's avg trapz 1->0", "Maksim's avg trapz 2->0", ...
    "Maksim's avg trapz 3->0", "Maksim's avg trapz 5->0",...
    "Adamovich98\_2 p10", "Adamovich98\_2 p20", "Adamovich98\_2 p30", ...
                            "Adamovich98\_2 p50"];
if calc_integralN
    legend_content=["Maksim's avg int 1->0", "Maksim's avg int 2->0", ...
    "Maksim's avg int 3->0", "Maksim's avg int 5->0", legend_content];
end
legend(legend_content, 'location', 'best')
xlim([1e3 1e6])     % M-A
ylim([1e-16 1])
title('N2-N2atom')
%% 4. Adamovich fig 9 comparison (from 40th level), calculation
i1=40;                                  % initial state
Emax=6e4;                               % max E value for fig 9
f1=[39 38 37 35];                       % final levels
% E_arr_40=(fr*Emax+(1-fr)*Emin)*100*h*c;    % splitting E interval
E_arr_40=10.^(log10(Emin):(log10(Emax)-log10(Emin))/fr_num:log10(Emax))...
                                                                *100*h*c;
M1=N2;                                  % the first molecule
M2=N2;
disp('Calculations for P_VT from 40th level.')
disp(['E array length is ', num2str(length(E_arr_40)) '.'])
disp('Collision with a fictitious N2a.')
if calc_integralN
 disp('IntegralN calculations started.')
 pif_40=zeros(length(f1), length(E_arr_40));
 for ind_f1=1:length(f1)
    disp(['i1=', num2str(i1), ', f1=', num2str(f1(ind_f1))])
    for ind=1:length(E_arr_40)
        pif_40(ind_f1, ind)=P_VT_FHO_FR_avg_angles(M1, M2, ...
                                        i1, f1(ind_f1), E_arr_40(ind));
    end
 end
end
E_cm_40=E_arr_40/h/c/100;     % switching from J to cm-1
%% 5. Adamovich fig 9 comparison (from 40th level), calculation trapz
pif_40_trapz=zeros(length(f1), length(E_arr_40));
disp('Trapz calculations started.')
for ind_f1=1:length(f1)
    disp(['i1=', num2str(i1), ', f1=', num2str(f1(ind_f1))])
    for ind=1:length(E_arr_40)
        pif_40_trapz(ind_f1, ind)=P_VT_FHO_FR_avg_angles(M1, M2, ...
                                    i1, f1(ind_f1), E_arr_40(ind), 't');
    end
end
%% 6. Adamovich fig 9 comparison (from 40th level), plotting
load data_Adamovich98_2_figs8_9.mat
figure
if calc_integralN
 loglog(E_cm_40, pif_40(1,:), E_cm_40, pif_40(2,:), ...
            E_cm_40, pif_40(3,:), E_cm_40, pif_40(4,:), 'linewidth', 1.5)
 hold on
end
loglog(E_cm_40, pif_40_trapz(1,:),':', E_cm_40,pif_40_trapz(2,:), ':',...
        E_cm_40, pif_40_trapz(3,:), ':', E_cm_40,pif_40_trapz(4,:),':',...
                                                        'linewidth', 1.5)
if ~calc_integralN
    hold on
end
loglog(p_VT_FHO_FR_N2N2a_4039_A98_2(:,1), ...
                                p_VT_FHO_FR_N2N2a_4039_A98_2(:,2), ...
        p_VT_FHO_FR_N2N2a_4038_A98_2(:,1), ...
                                p_VT_FHO_FR_N2N2a_4038_A98_2(:,2), ...
        p_VT_FHO_FR_N2N2a_4037_A98_2(:,1), ...
                                p_VT_FHO_FR_N2N2a_4037_A98_2(:,2), ...
        p_VT_FHO_FR_N2N2a_4035_A98_2(:,1), ...
                                p_VT_FHO_FR_N2N2a_4035_A98_2(:,2), ...
                                    'color', [0,0,0], 'linewidth', 1.5)
legend_content_40=["Maksim's avg trapz 40->39", ...
                                        "Maksim's avg trapz 40->38", ...
        "Maksim's avg trapz 40->37", "Maksim's avg trapz 40->35",...
        "Adamovich98\_2 p4039", "Adamovich98\_2 p4038", ...
        "Adamovich98\_2 p4037", "Adamovich98\_2 p4035"];

if calc_integralN
    legend_content_40=["Maksim's avg int 40->39",...
                                            "Maksim's avg int 40->38",...
                "Maksim's avg int 40->37", "Maksim's avg int 40->35",...
                                                    legend_content_40];
end
legend(legend_content_40, 'location', 'best')
xlim([1e3 6e4])     % M-A
ylim([1e-16 1])
title('N2-N2atom')
%% 7. Gimelshein 2017 fig. 1 left comparison, trapz
Emax=9e5;               % max E for fig 8
fr_num=50;
E_arr=10.^(log10(Emin):(log10(Emax)-log10(Emin))/fr_num:log10(Emax))...
                                                                *100*h*c;
M1=N2;                  % N2-N2a collision
M2=N2;                  % currently doesn't matter
i1=[1 5];               % initial levels
f1=0;                   % final state
pif_trapz=zeros(length(i1), length(E_arr));
disp('Trapz calculations started.')
for ind_i1=1:length(i1)
    disp(['i1=', num2str(i1(ind_i1)), ', f1=', num2str(f1)])
    for ind=1:length(E_arr)
        pif_trapz(ind_i1, ind)=P_VT_FHO_FR_avg_angles(M1, M2, ...
                                    i1(ind_i1), f1, E_arr(ind), 't');
    end
end
i1=40;
f1=[35 39];
pif_trapz_40=zeros(length(f1), length(E_arr));
for ind_f1=1:length(f1)
    disp(['i1=', num2str(i1), ', f1=', num2str(f1(ind_f1)), ])
    for ind=1:length(E_arr)
        pif_trapz_40(ind_f1, ind)=P_VT_FHO_FR_avg_angles(M1, M2, ...
                                    i1, f1(ind_f1), E_arr(ind), 't');
    end
end
%% 8. Gimelshein 2017 fig. 1 left comparison, plotting
E_cm=E_arr/h/c/100;             % E in cm-1
load data_Gimelshein2017_f1l    % loading of data from Gimelshein plot
figure
loglog(p_VT_FHO_FR_N2N2a_10_G17(:,1), p_VT_FHO_FR_N2N2a_10_G17(:, 2), ...
            'linewidth', 1.5)   % plotting Gimelshein plot data
hold on
loglog(p_VT_FHO_FR_N2N2a_50_G17(:,1), p_VT_FHO_FR_N2N2a_50_G17(:, 2),...
                                                        'linewidth', 1.5)
loglog(p_VT_FHO_FR_N2N2a_4039_G17(:,1), ...
                    p_VT_FHO_FR_N2N2a_4039_G17(:, 2), 'linewidth', 1.5)
loglog(p_VT_FHO_FR_N2N2a_4035_G17(:,1), ...
                    p_VT_FHO_FR_N2N2a_4035_G17(:, 2), 'linewidth', 1.5)
    % plotting trapz results
loglog(E_cm(1:ind), pif_trapz(1,:), '--', 'linewidth', 1.5)
loglog(E_cm(1:ind), pif_trapz(2,:), '--', 'linewidth', 1.5)
loglog(E_cm(1:ind), pif_trapz_40(1,:), '--', 'linewidth', 1.5)
loglog(E_cm(1:ind), pif_trapz_40(2,:), '--', 'linewidth', 1.5)
legend('FHO-FR Gimelshein et al. 2017, 1->0', ...
        'FHO-FR Gimelshein et al. 2017, 5->0', ...
        'FHO-FR Gimelshein et al. 2017, 40->39', ...
        'FHO-FR Gimelshein et al. 2017, 40->35', ...
        "Maksim's trapz code, 1->0", "Maksim's trapz code, 5->0", ...
        "Maksim's trapz code, 40->39", "Maksim's trapz code, 40->35", ...
                                                    'location', 'best')
xlim([1e3 9e5])
ylim([1e-15 1e-1])
title('N_2-N_2a p_{VT} FHO-FR')
%%
toc
