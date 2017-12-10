function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4Example6Problem()

% disable warning of `Support of character vectors will be removed in a future release.` 
% which is produced by function `monomials`.
warning('off')

% independent variables
syms v a d vf af;
vars = [v a d vf af];

% Constructing the vector field dx/dt = f
f = [a;
     -3*a-3*(v-vf)+(d-(v+10));
     vf-v;
     af;
     0];

eps = [0.001,0.001];

% Constructing the theta constraint
theta1 = d-4;
theta2 = d-5;
theta3 = v-vf+1;
theta4 = v-vf;
theta5 = a+1;
theta6 = a;
g_theta = [theta1,theta2,theta3,theta4,theta5,theta6];

% Constructing the psy constraint
psy1 = d;
psy2 = v;
psy3 = vf;
g_psy = [psy1, psy2, psy3];

% Constructing the zeta constraint
zeta1 = (a+2)/7;
zeta2 = (af+2)/7;
zeta3 = d+1;
zeta4 = d;
g_zeta = [zeta1,zeta2,zeta3,zeta4];

end
