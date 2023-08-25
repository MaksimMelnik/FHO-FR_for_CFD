function k_VT=k_VT_FHO_FR(T, M1, M2, i1, f1)
% FHO-FR k_VT function. Averaging P_VT_FHO_FR_avg_angles over E, depending
% on initial and final states of M1. Currently may be used only for M-A
% collision. Also currently w/o calculation optimization.
% T is the gas temperature; M1 is the molecule under consideration; M2 is
% the collision partner particle; i1 is the initial vibrational level;
% f1 is the final vibrational level.

k = 1.380649e-23;                   % Boltzmann constant, J/K
h = 6.626070041e-34;                % Plank constant, J*sec
c = 299792458;                      % speed of light, m/sec
coll.alpha= 4.0000e+10;
coll.red_mass = M1.mass*M2.mass/(M1.mass+M2.mass);
R = (M1.diameter+M2.diameter)/2;    % collision diameter
u = sqrt(8*k*T/pi/coll.red_mass);
T_cm1=T*k/h/c/100;                  % in cm-1
error_val=1e-7;
max_eps1 = @(eps1) 1-eps1;
if M2.num_vibr_levels(1)==1
    int=integral(@(E) (E./T_cm1).^2 .* P_VT_FHO_FR_avg_angles(M1, M2, ...
            i1, f1, E*h*c*100, 't') .* exp(-E./T_cm1)/T_cm1, 0, Inf, ...
                            'RelTol', error_val, 'AbsTol', error_val);
else
    max_angle = pi;
    
    int=integral(@(E) (E./T_cm1).^3 .* P_VT_FHO_FR_MM_ij_trapz(M1, M2, ...
            i1, f1, E*h*c*100)/pi^4 .* exp(-E./T_cm1)/T_cm1, 0, Inf, ...
                            'RelTol', error_val, 'AbsTol', error_val);

    % pvt=integral_8(@(phi2, phi1, theta_v2, theta_v1, y, eps2, eps1) ...
    
    % int=integral_8(@(E, eps1, eps2, y, theta_v1, theta_v2, phi1, phi2) ...
	%        (E./T_cm1).^3 .* P_VT_FHO_FR_MM_ij(M1, M2, i1, f1, E*h*c*100, eps1, y, ...
    %                              theta_v1, phi1, eps2, theta_v2, phi2)/pi^4 .* exp(-E./T_cm1)/T_cm1 , ...
    %         0, Inf,  0, 1, 0, 1,  0, 1, 0, max_angle, 0, max_angle, 0, max_angle, 0, max_angle, ... limits
    %                                 'RelTol',error_val, 'AbsTol',error_val);
    
    % pvt=pvt/pi^4;
    % int=integral(@(E) (E./T_cm1).^3 .* pvt .* exp(-E./T_cm1)/T_cm1, 0, Inf, ...
    %                         'RelTol', error_val, 'AbsTol', error_val);


end
k_VT = pi*R^2*u * int;
display(T)
display(k_VT)
end