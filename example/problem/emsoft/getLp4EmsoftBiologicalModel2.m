function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4EmsoftBiologicalModel2()

% form https://ths.rwth-aachen.de/research/projects/hypro/biological-model-ii/
% barrier with lambda = -1
% 0.4788*x1^2 + 0.5227*x1*x2 + 2.492*x1*x3 + 0.3089*x1*x4 + 1.391*x1*x5 - 0.114*x1*x6 - 0.8628*x1*x7 + 0.2594*x1*x8 + 0.1612*x1*x9 - 8.967*x1 + 0.5521*x2^2 + 0.8764*x2*x3 + 0.9352*x2*x4 - 0.3628*x2*x5 - 0.04419*x2*x6 + 0.4258*x2*x7 - 0.08721*x2*x8 - 0.004147*x2*x9 - 6.376*x2 - 0.1785*x3^2 + 1.022*x3*x4 - 0.05612*x3*x5 - 0.002406*x3*x6 - 1.058*x3*x7 + 0.197*x3*x8 + 0.2396*x3*x9 - 10.73*x3 + 0.2405*x4^2 - 0.1054*x4*x5 + 0.3233*x4*x6 + 1.207*x4*x7 - 0.4179*x4*x8 - 0.115*x4*x9 - 6.947*x4 - 0.4099*x5^2 + 0.5629*x5*x6 - 0.3395*x5*x7 - 0.3277*x5*x8 + 0.06521*x5*x9 + 1.357*x5 - 0.1549*x6^2 - 0.09829*x6*x7 - 0.06071*x6*x8 - 0.2566*x6*x9 - 0.2111*x6 - 0.04051*x7^2 + 0.4806*x7*x8 + 0.2038*x7*x9 + 0.6266*x7 - 0.02868*x8^2 + 0.002664*x8*x9 - 0.1777*x8 - 0.1492*x9^2 - 0.4161*x9 + 30.41

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
% g_psy = (vars+5)/10; % -5 <= x <= 5
g_psy = vars/5; % 0 <= x <= 5

% Constructing the zeta constraint
% zeta1 = 5*(x1-0.6);
% zeta2 = 5*x2;
% zeta3 = 5*(x3-1.6);
% zeta4 = 5*(x4-1.8);
% zeta5 = 5*x5;
% zeta6 = 5*x6;
% zeta7 = 5*x7;
% zeta8 = 5*x8;
% zeta9 = 5*x9;
% g_zeta = [zeta1, zeta2, zeta3, zeta4, zeta5, zeta6, zeta7, zeta8, zeta9];
g_zeta = (vars-1.8)*5;

end
