function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4EmsoftBiologicalModel2()

% form https://ths.rwth-aachen.de/research/projects/hypro/biological-model-ii/

% independent variables
syms x1 x2 x3 x4 x5 x6 x7 x8 x9;
vars = [x1, x2, x3, x4, x5, x6, x7, x8, x9];

% Constructing the vector field dx/dt = f
f = [3*x1 - x1*x6;
    x4 - x2*x6;
    x1*x6 - 3*x3;
    x2*x6 - x4;
    3*x3 + 5*x1 - x5;
    5*x5 + 3*x3 + x4 - x6*(x1 + x2 + 2*x8 + 1);
    5*x4 + x2 - 0.5*x7;
    5*x7 - 2*x6*x8 + x9 - 0.2*x8;
    2*x6*x8 - x9];

eps = [0.001, 0.001];

% Constructing the theta constraint
g_theta = (vars - 0.99)*50;

% Constructing the psy constraint
g_psy = (vars+5)/10; 

% Constructing the zeta constraint
zeta1 = 5*(x1-0.6);
zeta2 = 5*x2;
zeta3 = 5*(x3-1.6);
zeta4 = 5*(x4-1.8);
zeta5 = 5*x5;
zeta6 = 5*x6;
zeta7 = 5*x7;
zeta8 = 5*x8;
zeta9 = 5*x9;
g_zeta = [zeta1, zeta2, zeta3, zeta4, zeta5, zeta6, zeta7, zeta8, zeta9];


end
