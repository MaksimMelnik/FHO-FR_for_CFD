function pVT=P_VT_FHO_FR_MA_ij(M1, Coll, i1, f1, ...
                                            E, eps1, y, theta_v1, phi1)
% FHO-FR probability with a full array of arguments depending on the
% vibrational level. Diatom-atom collision.
% M1 is the first molecule, Coll is the collision, i1 is the initial
% vibrational state of M1, f1 is the final vibr state of M1, E is the
% collision energy (J), eps1 is the fraction of rotational energy, y is 
% the collision parameter, theta_v1 and phi1 are orientation angles.
% 30.08.2022 Maksim Melnik
% source: I. V. Adamovich et al., J. Chem. Phys. 109, 7711-7724 (1998) 
% https://doi.org/10.1063/1.477417

h = 6.626070041e-34;	% Plank constant, J*s
h_bar=h/(2*pi);         % reduced Planck constant, J*s
k = 1.380649e-23;       % Boltzmann constant, J/K
alpha=Coll.alpha;       % FHO parameter, 1/m
el_lvl=1;               % electronic level
e1=M1.ev_i{el_lvl}(i1+1); % initial state, J
e2=M1.ev_i{el_lvl}(f1+1); % final state, J
m0=M1.red_osc_mass;

s=abs(i1-f1);
m_r=Coll.red_mass;      % collision reduced mass
xi=m_r/m0/2;
omega=abs(e1-e2)./s/h_bar; % 1/с
theta_p=4*pi^2*omega.^2*m_r/(alpha^2*k); % К
theta=h_bar*omega/k;    % K
u=sqrt(2*E/m_r);        % m/s
gamma=max(0, -0.5*sin(2*theta_v1).*cos(phi1).*sqrt(xi*eps1) ...
                                                +sqrt((1-eps1).*(1-y)));
Q=theta_p*xi*cos(theta_v1).^2.*cos(phi1).^2 ...
    ./(4*theta.*sinh(pi*omega./(alpha*u*gamma)).^2);
ns=(factorial(max(i1,f1))./factorial(min(i1,f1))).^(1./s);
pVT=(ns.*Q).^s./factorial(s).^2 ...
                   .*exp(-2*ns.*Q./(s+1)-(ns.*Q).^2./((s+1).^2.*(s+2)));
end