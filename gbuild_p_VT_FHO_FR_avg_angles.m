% Строим графики P_VT_FHO_FR_avg_angles(T)
% 30.08.2022 Maksim Melnik

tic
clear fr E_arr E_cm %p4039 %p01

k = 1.380649e-23;       % Boltzmann constant, J/K
h = 6.626070041e-34;    % Plank constant, J*sec
c = 299792458;          % speed of light, m/sec
Emin=1e3; Emax=1e6;     % минимальное и макс значения на графике 
                        % Адамович2017, см^-1
                        % нужно умножить на hc по идее
i1=40;
Emax=6e4;       % for 40->39
i1=[1 2 3 5];
f1=0;
Emax=9e5;
fr=[0 0.0005 0.001 0.002 0.003 0.005 0.006 0.008 0.01 0.02 0.04 0.05 ...
                                0.055 0.06 0.07 0.08 0.1 0.2 0.5 0.8 1];
% fr=[0 0.0005 0.001 0.005 0.01 0.05 0.1 0.5 1];
E_arr=(fr*Emax+(1-fr)*Emin)*100*h*c;
% T=
M1=N2;
M2=N2;
disp(length(E_arr))
% p01=zeros(1, length(E_arr));
% p20=zeros(1, length(E_arr));
p30=zeros(1, length(E_arr));
% p4039=zeros(1, length(E_arr));
disp('collision with a fictitious N2a')
pif=zeros(length(i1), length(E_arr));
for ind_i1=1:length(i1)
    disp(['i1=', num2str(i1(ind_i1)), ', f1=', num2str(f1)])
    for ind=1:length(E_arr)
        disp(ind)
        p_VT=P_VT_FHO_FR_avg_angles(M1, M2, Coll_N2_N2a, ...
                                            i1(ind_i1), f1, E_arr(ind));
        pif(ind_i1, ind)=p_VT;
    end
end
E_cm=E_arr/h/c/100;
%% gbuild
load data_Adamovich98_2_figs8_9.mat
figure
% loglog(E_cm, p01, E_cm, p4039, 'linewidth', 1.5)
loglog(E_cm, pif(1,:), E_cm, pif(2,:), E_cm, pif(3,:), E_cm, pif(4,:), ...
                                                        'linewidth', 1.5)
hold on
loglog(p_VT_FHO_FR_N2N2a_10_A98_2(:,1), ...
                p_VT_FHO_FR_N2N2a_10_A98_2(:,2), ...
    p_VT_FHO_FR_N2N2a_20_A98_2(:,1), p_VT_FHO_FR_N2N2a_20_A98_2(:,2), ...
    p_VT_FHO_FR_N2N2a_30_A98_2(:,1), p_VT_FHO_FR_N2N2a_30_A98_2(:,2), ...
    p_VT_FHO_FR_N2N2a_50_A98_2(:,1), p_VT_FHO_FR_N2N2a_50_A98_2(:,2), ...
                                    'color', [0,0,0], 'linewidth', 1.5)
% legend('1->0', '40->39', 'location', 'best')
legend("Maksim's avg int 1->0", "Maksim's avg int 2->0", ...
    "Maksim's avg int 3->0", "Maksim's avg int 5->0",...
    'Adamovich98\_2 p10', ...
    'Adamovich98\_2 p20', 'Adamovich98\_2 p30', 'Adamovich98\_2 p50', ...
    'location', 'best')
% xlim([1e3 2e5])     % M-M
% ylim([1e-11 1e-1])
xlim([1e3 1e6])     % M-A
ylim([1e-16 1])
title('N2-N2atom')
%% From 40th level calculating
i1=40;
Emax=6e4;       % for 40->39
f1=[39 38 37 35];
E_arr=(fr*Emax+(1-fr)*Emin)*100*h*c;
% T=
M1=N2;
M2=N;
disp('From 40th level')
disp(['E array length ', num2str(length(E_arr))])
disp('collision with a fictitious N2a')
pif_40=zeros(length(f1), length(E_arr));
for ind_f1=1:length(f1)
    disp(['i1=', num2str(i1), ', f1=', num2str(f1(ind_f1))])
    for ind=1:length(E_arr)
        disp(ind)
        p_VT=P_VT_FHO_FR_avg_angles(M1, M2, Coll_N2_N2a, ...
                                            i1, f1(ind_f1), E_arr(ind));
        pif_40(ind_f1, ind)=p_VT;
    end
end
E_cm=E_arr/h/c/100;
%% From 40th level gbuild
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