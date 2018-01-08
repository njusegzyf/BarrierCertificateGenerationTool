function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4Example4Problem()

% From FM2016 ref 18 Safety Verification of Hybrid Systems by Constraint Propagation Based Abstraction Refinement
% Example FOCUS:
% Flow: (x1, x2) = (x1 − x2, x1 + x2)
% Empty jump relation
% Init: 2.5 ≤ x1 ≤ 3 ∧ x2 = 0
% Unsafe: x1 ≤ 2
% The state space: [0,4] × [0,4]

% independent variables
syms x1 x2;
vars = [x1 x2];

% Constructing the vector field dx/dt = f
f = [x1-x2;
     x1+x2];

eps = [0.1,0.1];

% Constructing the theta constraint
theta1 = 2*x1-5;
theta2 = -x2;
theta3 = x2;
g_theta = [theta1,theta2,theta3];

% Constructing the psy constraint
psy1 = x1/3;
psy2 = x2/3;
g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = x1/2;
zeta2 = x2/3;
g_zeta = [zeta1, zeta2];

end
