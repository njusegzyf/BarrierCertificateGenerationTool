function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4EmsoftBiologicalModel1()

% form https://ths.rwth-aachen.de/research/projects/hypro/biological-model-i/
% barrier for lambda = -1
% -2.549*x3-2.995*x6-3.576*x4-2.557*x1-3.389*x5-1.348*x2-2.615*x7+.754*x4^2+0.5552e-1*x4*x1+.7386*x5^2+0.4643e-1*x5*x7+.1847*x5*x2-.7219*x6*x3+.7259*x6^2+.3867*x1^2+1.299*x1*x2+.1235*x1*x3-.6037*x1*x5+.1977*x1*x6-.4632*x1*x7+0.4123e-1*x4*x2-.713*x4*x3+.7813*x4*x5+.215*x4*x6-.517*x4*x7-2.273*x2^2+1.446*x2*x3-.1724*x2*x6-0.3034e-1*x2*x7+0.9971e-1*x5*x3-.7006*x5*x6+0.9994e-1*x7*x3+.6089*x7*x6+.5814*x7^2+.4251*x3^2+22.71

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
% g_psy = (vars+2)/4; % -2 <= x <= 2
% g_psy = vars/2; % 0 <= x <= 2
g_psy = (vars + 2)/2; % -2 <= x <= 0

% Constructing the zeta constraint
% zeta1 = 5*(x1-0.6);
% zeta2 = 5*(x2-0.6);
% zeta3 = 5*(x3-1);
% zeta4 = 5*(x4-0.6);
% zeta5 = 5*(x5-1.2);
% zeta6 = 5*(x6-1);
% zeta7 = 5*(x7-0.6);
% g_zeta = [zeta1, zeta2, zeta3, zeta4, zeta5, zeta6, zeta7];
g_zeta = (vars-1.8)*5;

end
