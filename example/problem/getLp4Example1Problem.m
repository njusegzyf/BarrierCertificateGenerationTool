function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4Example1Problem()

% From FM2016 ref 18 Safety Verification of Hybrid Systems by Constraint Propagation Based Abstraction Refinement
% Example CLOCK: A simple example with a clock variable.
% Flow: (x, y, t) = (−5.5y + y^2, 6x − x^2, 1)
% Empty jump relation
% Init: 4 ≤ x ≤ 4.5 ∧ y = 1 ∧ t = 0
% Unsafe: (1 ≤ x < 2 ∧ 2 < y < 3 ∧ 2 ≤ t ≤ 4)
% The state space: [1, 5] × [1, 5] × [0, 4]

% independent variables
syms x1 x2;
vars = [x1 x2];

% Constructing the vector field dx/dt = f
f = [-(11/2)*x2+x2^2;
    6*x1-x1^2];

eps = [0.001, 0.001];

% Constructing the theta constraint
theta1 = 2*x1-8;
theta2 = x2;
theta3 = x2-1;

g_theta = [theta1, theta2, theta3];

% Constructing the psy constraint
psy1 = (x1-1)/4;
psy2 = (x2-1)/4;

g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = x1-1;
zeta2 = x2-2;

g_zeta = [zeta1, zeta2];

end
