function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4Emsoft6d2Problem()

%B=- 0.1182- 0.06094*x1- 0.00726*x2+ 0.3824*x3+ 0.3248*x4-0.1632*x5-0.1656*x6+ 0.03366*x1*x2- 0.006225*x1*x3- 0.004584*x1*x4- 0.009548*x1*x5+ 0.01588*x1*x6- 0.003326*x2*x3+ 0.005185*x2*x4+ 0.01407*x2*x5- 0.01693*x2*x6+ 0.04518*x3*x4  - 0.0167*x3*x5- 0.008248*x3*x6+ 0.004361*x4*x5+ 0.01167*x4*x6+ 0.0214*x5*x6  - 0.009488*x1^2- 0.009429*x2^2- 0.06402*x3^2- 0.08304*x4^2+0.005361*x5^2-0.007063*x6^2.

% independent variables
syms x1 x2 x3 x4 x5 x6;
vars = [x1; x2; x3; x4; x5; x6];

% Constructing the vector field dx/dt = f
f = [-x1^3+4*x2^3-6*x3*x4; -x1-x2+x5^3; x1*x4-x3+x4*x6; x1*x3+x3*x6-x4^3; -2*x2^3-x5+x6; -3*x3*x4-x5^3-x6];

eps = [0.0001,0.0001];

% Constructing the theta constraint
theta1 = 10*(x1-3);
theta2 = 10*(x2-3);
theta3 = 10*(x3-3);
theta4 = 10*(x4-3);
theta5 = 10*(x5-3);
theta6 = 10*(x6-3);
g_theta = [theta1,theta2,theta3,theta4,theta5,theta6];

% Constructing the psy constraint
psy1 = (x1)/10;
psy2 = (x2)/10;
psy3 = (x3-2)/8;
psy4 = (x4)/10;
psy5 = (x5)/10;
psy6 = (x6)/10;
g_psy = [psy1, psy2, psy3, psy4, psy5, psy6];

% Constructing the zeta constraint
zeta1 = 10*(x1-4);
zeta2 = 10*(x2-4.1);
zeta3 = 10*(x3-4.2);
zeta4 = 10*(x4-4.3);
zeta5 = 10*(x5-4.4);
zeta6 = 10*(x6-4.5);
g_zeta = [zeta1,zeta2,zeta3,zeta4,zeta5,zeta6];

end
