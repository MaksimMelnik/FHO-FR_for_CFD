function pvt=P_VT_FHO_FR_avg_angles(M1, M2, Coll, i1, f1, E)
% The FRO-FR probability calculation function averaging P_VT over all
% angles and parameters except of collision energy E.
% 24.08.2022 Maksim Melnik

% P_VT_FHO_FR_MM(M1, Coll, E, eps1, y, theta_v1, phi1, ...
%                                                    eps2, theta_v2, phi2) 
% int_1=@(y, eps1, theta_v1, phi1) ...
%     integral3(@(eps2, theta_v2, phi2) P_VT_FHO_FR_MM(M1, Coll, E, ...
%     eps1, y, theta_v1, phi1, eps2, theta_v2, phi2), 0, 1, 0, pi, 0, pi);
% int_2=@(y) integral3(@(eps1, theta_v1, phi1) int_1())

% RelTol=1e-1;
% switch M2.num_atoms
%     case 1
% %         P_VT_FHO_FR_MA(M1, Coll, E, eps1, y, theta_v1, phi1)
%         int_1=@(eps1, y, theta_v1) integral(@(phi1) ...
%             P_VT_FHO_FR_MA(M1, Coll, E, eps1, y, theta_v1, phi1), ...
%                         0, pi, 'ArrayValued', true, 'RelTol', RelTol);
%         int_2=@(eps1, y) ...
%             integral(@(theta_v1) int_1(eps1, y, theta_v1), 0, pi, ...
%                                 'ArrayValued', true, 'RelTol', RelTol);
%         int_3=@(y) ...
%              integral(@(eps1) int_2(eps1, y), 0, 1, ...
%                                 'ArrayValued', true, 'RelTol', RelTol);
%         int_4=integral(@(y) int_3(y), 0, 1, ...
%                                 'ArrayValued', true, 'RelTol', RelTol);
%         pvt=int_4;
%     case 2
% int_1=@(y, eps1, theta_v1, phi1, eps2, theta_v2) ...
%     integral(@(phi2) P_VT_FHO_FR_MM(M1, Coll, E, eps1, y, theta_v1, ...
%                 phi1, eps2, theta_v2, phi2), 0, pi, 'ArrayValued', true);
% int_2=@(y, eps1, theta_v1, phi1, eps2) ...
%     integral(@(theta_v2) int_1(y, eps1, theta_v1, phi1, eps2, ...
%                                   theta_v2), 0, pi, 'ArrayValued', true);
% int_3=@(y, eps1, theta_v1, phi1) integral(@(eps2) ...
%        int_2(y, eps1, theta_v1, phi1, eps2),  0, 1, 'ArrayValued', true);
% int_4=@(y, eps1, theta_v1) integral(@(phi1) ...
%              int_3(y, eps1, theta_v1, phi1), 0, pi, 'ArrayValued', true);
% int_5=@(y, eps1) integral(@(theta_v1) ...
%                    int_4(y, eps1, theta_v1), 0, pi, 'ArrayValued', true);
% int_6=@(y) integral(@(eps1) int_5(y, eps1), 0, 1, 'ArrayValued', true);
% int_7=integral(@(y) int_6(y), 0, 1, 'ArrayValued', true);
% pvt=int_7;
% end

% i1=1;
% f1=i1-1;

% my nested integrals
% int_1=@(y) integral3(@(eps1, theta_v1, phi1) ...
%     P_VT_FHO_FR_MA_ij(M1, Coll, i1, f1, E, eps1, y, theta_v1, phi1), ...
%                                                     0, 1, 0, pi, 0, pi);
% % pvt=int_1(0);
% pvt=integral(@(y) int_1(y), 0, 1, 'ArrayValued', true);

% % MATLAB integralN
max_angle=pi;
error_v=1e-2;
pvt=integralN(@(y, eps1, theta_v1, phi1) ...
	P_VT_FHO_FR_MA_ij(M1, Coll, i1, f1, E, eps1, y, theta_v1, phi1),...
        0, 1, 0, 1, 0, max_angle, 0, max_angle, ...
                                    'RelTol',error_v,'AbsTol',error_v);
end