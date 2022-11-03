% The initialization of particles and collisions valiables.
    %% the constants
h = 6.626070041e-34;      % Plank constant, J*sec
N_a=6.02214076e23;        % Avogadro constant
c = 299792458;            % speed of light, m/sec
k = 1.380649e-23;         % Boltzmann constant, J/K

    %% particles data
C.name='C';
C.mass=1.994473440944e-26;          % kg
C.diameter=3.298e-10;               % m
C.s_e=1;
C.m_mass=12.01;                     % molar mass
C.form_e=4.97472e-19;               % online table for kappa
C.form_e=1.180914618115280e-18;     % by Alena Kosareva
C.fr_deg_c=3;                       % freedom degree for room temperature
C.num_atoms=1;                      % number of atoms in the particle

O.name='O';
O.mass=2.6567628316576e-26; O.diameter=2.75e-10;    O.s_e=9;
O.m_mass=16;                        % molar mass
O.form_e= 4.098045681049634e-19;    % online table for kappa
O.form_e= 4.098111014876693e-19;    % by Alena Kosareva
O.fr_deg_c=3;                       % freedom degree for room temperature
O.num_atoms=1;                      % number of atoms in the particle

N.name='N';
N.mass=2.32587E-26;                 % kg
N.diameter=3.29800E-10;             % m
N.num_atoms=1;                      % number of atoms in the particle

O2.name='O2';
O2.mass=5.31353E-26;                % kg
O2.num_elex_levels=1;               % number of electronic level
O2.num_vibr_levels=[36];  % numbers of vibrational states for each el lvl
O2.we=[158019];                     % Frequency of vibrations, 1/m
O2.wexe=[1198];
O2.weye=[4.747];
e_i=levels_e_ex(O2, 1);             
O2.ev_0(1)=e_i(1);  O2.ev_i{1}=e_i-e_i(1);  % vibrational energy
O2.num_atoms=2;                     % number of atoms in the particle

N2.name='N2';
N2.mass=4.65173E-26;                % kg
N2.num_elex_levels=1;
N2.num_vibr_levels=[47];
N2.we=[235857];                     % 1/m
N2.wexe=[1432.40];                  % 1/m
N2.weye=0;%[-0.226];                % 1/m
e_i=levels_e_ex(N2, 1);
N2.ev_0(1)=e_i(1);  N2.ev_i{1}=e_i-e_i(1);
N2.diss_e=[1.56362E-18];            % J
N2.diameter=3.4039E-10;             % m
N2.red_osc_mass=N.mass/2;           % kg
N2.num_atoms=2;                     % number of atoms in the particle

CO.name='CO';
CO.mass=4.651236272601599e-26;            % molecular mass, kg
CO.m_mass=C.m_mass+O.m_mass;              % molar mass
CO.red_osc_mass=C.mass*O.mass/(C.mass+O.mass);    
% CO.diss_e=1.77716869535000e-18;         % dissociation energy
CO.diss_e=[89460 41121 24848]*100*h*c; % dissociation energies, Krupenie
CO.form_e=-1.8349e-19;                    % online table for kappa
CO.form_e=-1.864225959573993e-19;         % by Alena Kosareva
CO.form_e_atoms_sum=C.form_e+O.form_e;
CO.num_elex_levels=3;                     % number of electronical levels
CO.num_vibr_levels=[68 34 18]; % numbers of vibr levels for each el. state
CO.diameter=3.65e-10;                     % molecular diameter
CO.sigma=1;     % symmetry factor (1 for heteronuclear, 2 for homonuclear)
CO.Be=[193.128087 169.124 161.15];        % spectroscopic constant
CO.moment_I=1.4560E-46;
CO.mltpl_atoms_mass=C.mass*O.mass;        % product of atoms masses
CO.mltpl_atoms_s_e=C.s_e*O.s_e;           % product of statistical weights
CO.s_e=[1 6 2];                           % electronic statistical weight
CO.we = [216981.358 174341 151824];       % spectroscopic constants, 1/m
CO.wexe=[1328.831 1436 1940]; 
CO.weye=[1.0511 -4.5 76.6];         
CO.e_E=[0  4868740*h*c  6507460*h*c]; % el. excitation energy, Krupenie
e_i=levels_e_ex(CO, 1);                   % vibrational levels energy
CO.ev_0(1)=e_i(1);  CO.ev_i{1}=e_i-e_i(1);  % the new e_vibr standart
e_i=levels_e_ex(CO, 2);  
CO.ev_0(2)=e_i(1);  CO.ev_i{2}=e_i-e_i(1);  % the new e_vibr standart
e_i=levels_e_ex(CO, 3);  
CO.ev_0(3)=e_i(1);  CO.ev_i{3}=e_i-e_i(1);  % the new e_vibr standart
CO.C_VE=[0 13.58e-2 5.49e-2];                   % VE-exchange parameters
CO.S_VE=[0 5e-13 9.94e-13]/1e6;                 % in SI
CO.mA_mAB=0.42880501528003884;  CO.mB_mAB=0.5711949847199612;  % FHO param
CO.fr_deg_c=5;                      % freedom degree for room temperature
CO.num_atoms=2;                     % number of atoms in the particle

C2.name='C2';
C2.mass=3.988946881888e-26;         C2.diameter=3.621e-10; 
C2.m_mass=2*C.m_mass;               % molar mass
C2.diss_e=9.949436628352846e-19;    C2.ev_i{1}=0;      C2.ev_0=0;
% C2.e_i=0; 
C2.e_E=0;
C2.Be=181.98399999999998;   C2.sigma=2;     C2.mltpl_atoms_mass=C.mass^2;
C2.s_e=1;   C2.e_E=0;   C2.num_elex_levels=1;   C2.num_vibr_levels=1;
C2.fr_deg_c=5;                      % freedom degree for room temperature
C2.mltpl_atoms_mass=C.mass*C.mass;
C2.mltpl_atoms_s_e=C.s_e*C.s_e;
C2.form_e=0;                        % Formation energy
C2.form_e_atoms_sum=C.form_e+C.form_e;
C2.CO_C_C2_O_e=58000*k;
C2.num_atoms=2;                     % number of atoms in the particle

Ar.name='Ar';
Ar.mass=6.633521356992e-26; Ar.diameter=3.33e-10;
Ar.m_mass=39.95;                    % molar mass
Ar.fr_deg_c=3;                      % freedom degree for room temperature
Ar.num_atoms=1;                     % number of atoms in the particle



    %% collisions data
% CO-M
Coll_CO_CO.coll_diameter=CO.diameter;
    % Park
Coll_CO_CO.ArrA(1)=2.3e20/N_a*1e-6;     Coll_CO_CO.ArrN(1)=-1;
Coll_CO_C.ArrA(1) =3.4e20/N_a*1e-6;     Coll_CO_C.ArrN(1) =-1;
Coll_CO_O.ArrA(1) =3.4e20/N_a*1e-6;     Coll_CO_O.ArrN(1) =-1;
Coll_CO_Ar.ArrA(1)=2.3e19/N_a*1e-6;     Coll_CO_Ar.ArrN(1)=-1;
Coll_CO_C2.ArrA(1)=2.3e20/N_a*1e-6;     Coll_CO_C2.ArrN(1)=-1;
Coll_C2.ArrA(1)   =3.7e14/N_a*1e-6;     Coll_C2.ArrN(1)   = 0;
Coll_CO_C__C2_O.ArrA(1)=2e17/N_a*1e-6;  Coll_CO_C__C2_O.ArrN(1)=-1;
    % Ibraguimova
Coll_CO_CO.ArrA(2)=2.8e21/N_a*1e-6;     Coll_CO_CO.ArrN(2)=-1.39;
Coll_CO_C.ArrA(2) =1.4e22/N_a*1e-6;     Coll_CO_C.ArrN(2) =-1.39;
Coll_CO_O.ArrA(2) =2.1e22/N_a*1e-6;     Coll_CO_O.ArrN(2) =-1.39;
Coll_CO_Ar.ArrA(2)=1.4e21/N_a*1e-6;     Coll_CO_Ar.ArrN(2)=-1.39;
    % McKenzie
Coll_CO_CO.ArrA(3)=0;	Coll_CO_CO.ArrN(3)=0;
    % Fairbairn diss coefficients
CO_CO_krec=2e-34/1e12;  % cm6 s-1, Fairbairn
Tfb=7600;               % Temperature in Fairbairn experiment
Zvibr=sum(exp(-CO.ev_i{1}'/(k*Tfb)));
Io2=1.95790474021232e-46;
Zrot=8*pi^2*Io2*k*Tfb/(2*h^2);
Zint=Zrot.*Zvibr;
kdr_fb=(CO.mass/CO.mltpl_atoms_mass)^1.5*h^3*...
            (2*pi*k*Tfb).^(-1.5).*Zint.*exp(CO.diss_e(1)./(k*Tfb));
Coll_CO_CO.ArrA(4)=CO_CO_krec/kdr_fb;    Coll_CO_CO.ArrN(4)=0;
Coll_CO_C.ArrA(4) =CO_CO_krec/kdr_fb;    Coll_CO_C.ArrN(4) =0;
Coll_CO_O.ArrA(4) =CO_CO_krec/kdr_fb;    Coll_CO_O.ArrN(4) =0;
Coll_CO_Ar.ArrA(4)=CO_CO_krec/kdr_fb;    Coll_CO_Ar.ArrN(4)=0;
Coll_CO_C2.ArrA(4)=CO_CO_krec/kdr_fb;    Coll_CO_C2.ArrN(4)=0;
Coll_C2.ArrA(4)   =6e-10/1e6;            Coll_C2.ArrN(4)   =0;
Coll_CO_C__C2_O.ArrA(4)=6e-10/1e6;       Coll_CO_C__C2_O.ArrN(4)=0;
    % Appleton diss coefficients with global diss energy
Coll_CO_CO.ArrA(5)=2.8e-3/1e6;      Coll_CO_CO.ArrN(5)=-2.86;
Coll_CO_C.ArrA(5) =2.8e-3/1e6;      Coll_CO_C.ArrN(5) =-2.86;
Coll_CO_O.ArrA(5) =2.8e-3/1e6;      Coll_CO_O.ArrN(5) =-2.86;
Coll_CO_Ar.ArrA(5)=2.8e-3/1e6;      Coll_CO_Ar.ArrN(5)=-2.86;
Coll_CO_C2.ArrA(5)=2.8e-3/1e6;      Coll_CO_C2.ArrN(5)=-2.86;
Coll_C2.ArrA(5)=3.7e14/N_a*1e-6;    Coll_C2.ArrN(5)   =0;
    % Appleton diss coefficients with own diss energy
% CO.diss_e(1)=98600*k;
%     Coll_CO_CO.ArrA(6)=4.4e-10/1e6;     Coll_CO_CO.ArrN(6)=0;
%     Coll_CO_C.ArrA(6) =4.4e-10/1e6;     Coll_CO_C.ArrN(6) =0;
%     Coll_CO_O.ArrA(6) =6.86e-9/1e6;     Coll_CO_O.ArrN(6) =0;
%     Coll_CO_Ar.ArrA(6)=4.4e-10/1e6;     Coll_CO_Ar.ArrN(6)=0;
%     Coll_CO_C2.ArrA(6)=4.4e-10/1e6;     Coll_CO_C2.ArrN(6)=0;
%     Coll_C2.ArrA(6)=3.7e14/N_a*1e-6;    Coll_C2.ArrN(6)=   0;
    % Mick diss coefficients=32 
Coll_CO_CO.ArrA(7)=4.3e27/N_a*1e-6;     Coll_CO_CO.ArrN(7)=-3.1;
Coll_CO_C.ArrA(7) =4.3e27/N_a*1e-6;     Coll_CO_C.ArrN(7) =-3.1;
Coll_CO_O.ArrA(7) =4.3e27/N_a*1e-6;     Coll_CO_O.ArrN(7) =-3.1;
Coll_CO_Ar.ArrA(7)=4.3e27/N_a*1e-6;     Coll_CO_Ar.ArrN(7)=-3.1;
Coll_CO_C2.ArrA(7)=4.3e27/N_a*1e-6;     Coll_CO_C2.ArrN(7)=-3.1;
Coll_C2.ArrA(7)   =3.72e14/N_a*1e-6;    Coll_C2.ArrN(7)   =0;
Coll_CO_C__C2_O.ArrA(7)=6e-10/1e6;       Coll_CO_C__C2_O.ArrN(7)=0; 

% N2-N2
Coll_N2_N2.red_mass=N2.mass/2;
Coll_N2_N2.alpha=3.7e10;

% N2-N
    Coll_N2_N.red_mass =N2.mass*N.mass/(N2.mass+N.mass);
    Coll_N2_N.alpha =4e10;
    
% N2-N2a (fictitious atom of N2's mass and diameter)
    Coll_N2_N2a.red_mass=N2.mass/2;
    Coll_N2_N2a.alpha=4e10;
    
save par_data O2 O N2 N Coll_N2_N2 Coll_N2_N Coll_N2_N2a