function [vars, stateNum, fs, eps, thetaStateIndex, g_theta, g_psys, g_zeta, g_guards] = getLp4HybridEx4Problem()

% independent variables
syms x1 x2 x3;
vars = [x1 x2 x3];

stateNum = 2;

% Constructing the vector field dx/dt = f
f1 = [x2;    -x1+x3;    x1+(2*x2+3*x3)*(1+x3^2)];
f2 = [x2;    -x1+x3;    -x1-2*x2-3*x3];
fs = [f1, f2];
% Note: use `fs(:, 1)` to get f1

import lp4.Lp4Config
% eps = [Lp4Config.DEFAULT_EPS, Lp4Config.DEFAULT_EPS];
eps = [0.01, 0.01];


% Constructing the theta constraint
theta1 = (x1^2+x2^2+x3^2)/0.01;
theta2 = (x1+0.1)/(0.2);
theta3 = (x2+0.1)/(0.2);
theta4 = (x3+0.1)/(0.2);
g_theta = [theta1, theta2, theta3, theta4];
thetaStateIndex = 1;



% Constructing the psy constraint
psy11 = (x1^2+0.01*x2^2+0.01*x3^2)/1.01;
psy12 = (x1+1.01)/(2.02);
psy13 = (0.1*x2+1.01)/(2.02);
psy14 = (0.1*x3+1.01)/(2.02);
g_psy1 = [psy11, psy12, psy13, psy14];

psy21 = (x1+5.1)/10.2;
psy22 = (x2+20)/(40);
psy23 = (x3+20)/(40);
psy24 = (x1^2+x2^2+x3^2-0.03)/800;
g_psy2 = [psy21, psy22, psy23, psy24];

g_psys = [g_psy1; g_psy2];



% Constructing the Guard constraint
import lp4.Guard

guard121=(x1^2+0.01*x2^2+0.01*x3^2-0.09)/0.02;
guard122=(x1+1.01)/(2.02);
guard123=(0.1*x2+1.01)/(2.02);
guard124=(0.1*x2+1.01)/(2.02); 
guard12 = Guard(1, 2, [guard121, guard122, guard123, guard124], [], []);

guard211=(x1^2+x2^2+x3^2-0.03)/0.02;
guard212=(x1+2.25)/4.5;
guard213=(x2+2.25)/4.5;
guard214=(x3+2.25)/4.5;
guard21 = Guard(2, 1, [guard211, guard212, guard213, guard214], [], []);

g_guards = [guard12, guard21];



% Constructing the zeta constraint
import lp4.UnsafeConstraint
zeta11 = (x2)/(12);
zeta12 = (x3)/(12);
zeta13 = (x1-5)/0.1;
g_zeta = [UnsafeConstraint(2, [zeta11, zeta12, zeta13])];

end
