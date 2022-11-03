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
%% Adamovich fig 8 comparison (to 0th level), calculation
i1=[1 2 3 5];           % initial levels
f1=0;                   % final state
Emax=9e5;               % max E for fig 8
fr=[0 0.0005 0.001 0.002 0.003 0.005 ... the splitiing E interval array
  0.006 0.008 0.01 0.02 0.04 0.05 0.055 0.06 0.07 0.08 0.1 0.2 0.5 0.8 1];
% fr=[0 0.0005 0.001 0.005 0.01 0.05 0.1 0.5 1];
E_arr=(fr*Emax+(1-fr)*Emin)*100*h*c;    % splitting of E interval
M1=N2;                  % N2-N2a collision
M2=N2;                  % currently doesn't matter
disp(['Length of E is ' num2str(length(E_arr))])
disp('collision of N2 with a fictitious N2a')
pif=zeros(length(i1), length(E_arr));
for ind_i1=1:length(i1)
    disp(['i1=', num2str(i1(ind_i1)), ', f1=', num2str(f1)])
    for ind=1:length(E_arr)
        disp(ind)
        pif(ind_i1, ind)=P_VT_FHO_FR_avg_angles(M1, M2, Coll_N2_N2a, ...
                                            i1(ind_i1), f1, E_arr(ind));
    end
end
E_cm=E_arr/h/c/100;     % switching from J to cm-1
%% Adamovich fig 8 comparison (to 0th level), plotting
load data_Adamovich98_2_figs8_9.mat
figure
loglog(E_cm, pif(1,:), E_cm, pif(2,:), E_cm, pif(3,:), E_cm, pif(4,:), ...
                                                        'linewidth', 1.5)
hold on
loglog(p_VT_FHO_FR_N2N2a_10_A98_2(:,1), ...
                                     p_VT_FHO_FR_N2N2a_10_A98_2(:,2), ...
    p_VT_FHO_FR_N2N2a_20_A98_2(:,1), p_VT_FHO_FR_N2N2a_20_A98_2(:,2), ...
    p_VT_FHO_FR_N2N2a_30_A98_2(:,1), p_VT_FHO_FR_N2N2a_30_A98_2(:,2), ...
    p_VT_FHO_FR_N2N2a_50_A98_2(:,1), p_VT_FHO_FR_N2N2a_50_A98_2(:,2), ...
                                    'color', [0,0,0], 'linewidth', 1.5)
legend("Maksim's avg int 1->0", "Maksim's avg int 2->0", ...
    "Maksim's avg int 3->0", "Maksim's avg int 5->0",...
    'Adamovich98\_2 p10', 'Adamovich98\_2 p20', 'Adamovich98\_2 p30', ...
                            'Adamovich98\_2 p50', 'location', 'best')
xlim([1e3 1e6])     % M-A
ylim([1e-16 1])
title('N2-N2atom')
%% Adamovich fig 9 comparison (from 40th level), calculation
i1=40;                                  % initial state
Emax=6e4;                               % max E value for fig 9
f1=[39 38 37 35];                       % final levels
E_arr=(fr*Emax+(1-fr)*Emin)*100*h*c;    % splitting E interval
M1=N2;                                  % the first molecule
M2=N;                                   % currently doesn't matter
disp('From 40th level')
disp(['E array length ', num2str(length(E_arr))])
disp('collision with a fictitious N2a')
pif_40=zeros(length(f1), length(E_arr));
for ind_f1=1:length(f1)
    disp(['i1=', num2str(i1), ', f1=', num2str(f1(ind_f1))])
    for ind=1:length(E_arr)
        disp(ind)
        pif_40(ind_f1, ind)=P_VT_FHO_FR_avg_angles(M1, M2, ...
                                Coll_N2_N2a, i1, f1(ind_f1), E_arr(ind));
    end
end
E_cm=E_arr/h/c/100;     % switching from J to cm-1
%% Adamovich fig 9 comparison (from 40th level), plotting
load data_Adamovich98_2_figs8_9.mat
figure
loglog(E_cm, pif_40(1,:), E_cm, pif_40(2,:), E_cm, pif_40(3,:), ...
                                    E_cm, pif_40(4,:), 'linewidth', 1.5)
hold on
loglog(p_VT_FHO_FR_N2N2a_4039_A98_2(:,1), ...
                                p_VT_FHO_FR_N2N2a_4039_A98_2(:,2), ...
        p_VT_FHO_FR_N2N2a_4038_A98_2(:,1), ...
                                p_VT_FHO_FR_N2N2a_4038_A98_2(:,2), ...
        p_VT_FHO_FR_N2N2a_4037_A98_2(:,1), ...
                                p_VT_FHO_FR_N2N2a_4037_A98_2(:,2), ...
        p_VT_FHO_FR_N2N2a_4035_A98_2(:,1), ...
                                p_VT_FHO_FR_N2N2a_4035_A98_2(:,2), ...
                                    'color', [0,0,0], 'linewidth', 1.5)
legend("Maksim's avg int 40->39", "Maksim's avg int 40->38", ...
        "Maksim's avg int 40->37", "Maksim's avg int 40->35",...
        'Adamovich98\_2 p4039', 'Adamovich98\_2 p4038', ...
        'Adamovich98\_2 p4037', 'Adamovich98\_2 p4035', ...
                                                    'location', 'best')
xlim([1e3 6e4])     % M-A
ylim([1e-16 1])
title('N2-N2atom')
%%
toc
