function pvt=P_VT_FHO_FR_avg_angles(M1, M2, Coll, i1, f1, E, method)
% The FRO-FR probability calculation function averaging P_VT over all
% angles and parameters except of collision energy E, depending on initial 
% and final vibrational states of M1. Currently may be used only for M-A 
% collision.
% M1 is the molecule under consideration, M2 is a collision partner
% particle, Coll is the collision variable, i1 is the initial vibrational
% level, f1 is the final vibrational level, E is the collision energy (J),
% method is integration method ('q' for integralN, 't' for trapz).
% src: eq.15, I. V. Adamovich et al., J. Chem. Phys. 109, 7711-7724 (1998) 
% https://doi.org/10.1063/1.477417
% NOTE: to use present function multiple integration function integralN 
% is required:
% https://www.mathworks.com/matlabcentral/fileexchange/47919-integraln-m
% 24.08.2022 Maksim Melnik

if nargin<7
   method='q';
end
max_angle=pi; % the upper limit of theta_v1 and phi1 angles integration
switch method
    case 'q'
        error_v=1e-2; % integration error
        pvt=integralN(@(y, eps1, theta_v1, phi1) ...
	        P_VT_FHO_FR_MA_ij(M1, Coll, i1, f1, E, eps1, y, ...
                                                    theta_v1, phi1), ...
            0, 1, 0, 1, 0, max_angle, 0, max_angle, ... limits
                                    'RelTol',error_v, 'AbsTol',error_v);
    case 't'
        % integration over epsilon, theta, phi and y via trapz
        maxdiv=25;  % 25 seems fine, 32 times faster
        epsilonmax=1; ymax=1; theta_v1max=max_angle; phi1max=max_angle;
        eps1=0:epsilonmax/maxdiv:epsilonmax;    % 1
        eps1rs=eps1';
        y=0:ymax/maxdiv:ymax;                   % 2
        yrs=y;
        theta_v1=0:theta_v1max/maxdiv:theta_v1max;
        theta_v1rs=reshape(theta_v1, 1, 1, []); % 3
        phi1=0:phi1max/maxdiv:phi1max;
        phi1rs=reshape(phi1, 1, 1, 1, []);      % 4
        out=P_VT_FHO_FR_MA_ij(M1, Coll, i1, f1, E, eps1rs, yrs, ...
                                                theta_v1rs, phi1rs);
%         disp(length(out))
        pvt=trapz(phi1, trapz(theta_v1, trapz(y, ...
                                        trapz(eps1, out, 1), 2), 3), 4);
end
end
