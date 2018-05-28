function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4EmsoftBiologicalModel1()

% form https://ths.rwth-aachen.de/research/projects/hypro/biological-model-i/

% independent variables
syms x1 x2 x3 x4 x5 x6 x7;
vars = [x1, x2, x3, x4, x5, x6, x7];

% Constructing the vector field dx/dt = f
f = [-0.4*x1 + 5*x3*x4;
    0.4*x1 - x2;
    x2 - 5*x3*x4;
    5*x5*x6 - 5*x3*x4;
    -5*x5*x6 + 5*x3*x4;
    0.5*x7 - 5*x5*x6;
    -0.5*x7 + 5*x5*x6];

eps = [0.001, 0.001];

% Constructing the theta constraint
g_theta = (vars - 0.99)*50;

% Constructing the psy constraint
g_psy = (vars+2)/4; 

% Constructing the zeta constraint
zeta1 = 5*(x1-0.6);
zeta2 = 5*(x2-0.6);
zeta3 = 5*(x3-1);
zeta4 = 5*(x4-0.6);
zeta5 = 5*(x5-1.2);
zeta6 = 5*(x6-1);
zeta7 = 5*(x7-0.6);
g_zeta = [zeta1, zeta2, zeta3, zeta4, zeta5, zeta6, zeta7];


end
