% The example of usage and plotting of P_VT_FHO_FR_avg_angles(T) vs data
% from I. V. Adamovich et al., J. Chem. Phys. 109, 7711-7724 (1998) 
% https://doi.org/10.1063/1.477417
% For figs 8 and 9
% 30.08.2022 Maksim Melnik

tic
% clear fr E_arr E_cm %p4039 %p01

k = 1.380649e-23;       % Boltzmann constant, J/K
h = 6.626070041e-34;    % Plank constant, J*sec
c = 299792458;          % speed of light, m/sec
load par_data.mat       % loading particles and collisions
Emin=1e3; Emax=1e6;     % standart min and max E value on the plot
calc_integralN=false;   % use integralN method for calculation or 
                        %                                     only trapz?
%% Adamovich fig 8 comparison (to 0th level), calculation
i1=[1 2 3 5];           % initial levels
f1=0;                   % final state
Emax=9e5;               % max E for fig 8
fr_num=50;
E_arr=10.^(log10(Emin):(log10(Emax)-log10(Emin))/fr_num:log10(Emax))...
                                                                *100*h*c;
fr=[0 0.0005 0.001 0.002 0.003 0.005 ... the splitiing E interval array
  0.006 0.008 0.01 0.02 0.04 0.05 0.055 0.06 0.07 0.08 0.1 0.2 0.5 0.8 1];
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
        pif(ind_i1, ind)=P_VT_FHO_FR_avg_angles(M1, M2, Coll_N2_N2a, ...
                                            i1(ind_i1), f1, E_arr(ind));
    end
 end
end
E_cm=E_arr/h/c/100;     % switching from J to cm-1
%% Adamovich fig 8 comparison (to 0th level), calculation with trapz
pif_trapz=zeros(length(i1), length(E_arr));
disp('Trapz calculations started.')
for ind_i1=1:length(i1)
    disp(['i1=', num2str(i1(ind_i1)), ', f1=', num2str(f1)])
    for ind=1:length(E_arr)
        pif_trapz(ind_i1, ind)=P_VT_FHO_FR_avg_angles(M1, M2, ...
                        Coll_N2_N2a, i1(ind_i1), f1, E_arr(ind), 't');
    end
end
%% Adamovich fig 8 comparison (to 0th level), plotting
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
%% Adamovich fig 9 comparison (from 40th level), calculation
i1=40;                                  % initial state
Emax=6e4;                               % max E value for fig 9
f1=[39 38 37 35];                       % final levels
% E_arr_40=(fr*Emax+(1-fr)*Emin)*100*h*c;    % splitting E interval
E_arr_40=10.^(log10(Emin):(log10(Emax)-log10(Emin))/fr_num:log10(Emax))...
                                                                *100*h*c;
M1=N2;                                  % the first molecule
M2=N;                                   % currently doesn't matter
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
                                Coll_N2_N2a, i1, f1(ind_f1), E_arr_40(ind));
    end
 end
end
E_cm_40=E_arr_40/h/c/100;     % switching from J to cm-1
%% Adamovich fig 9 comparison (from 40th level), calculation trapz
pif_40_trapz=zeros(length(f1), length(E_arr_40));
disp('Trapz calculations started.')
for ind_f1=1:length(f1)
    disp(['i1=', num2str(i1), ', f1=', num2str(f1(ind_f1))])
    for ind=1:length(E_arr_40)
        pif_40_trapz(ind_f1, ind)=P_VT_FHO_FR_avg_angles(M1, M2, ...
                        Coll_N2_N2a, i1, f1(ind_f1), E_arr_40(ind), 't');
    end
end
%% Adamovich fig 9 comparison (from 40th level), plotting
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
%%
toc
