function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4Example5Problem()

% From FM2016 ref 18 Safety Verification of Hybrid Systems by Constraint Propagation Based Abstraction Refinement
% Example ECO: A predator-prey example of ecosystem problems.
% Flow: ¡s = 1 → ¡x x ˙ ˙ 1 2¢ = ¡− xx 21 − + xx 11 xx 22¢¢ ∧ ¡s = 2 → ¡x x ˙ ˙ 1 2¢ = ¡− xx 21 − + xx 11 xx 22¢¢
% Jump: ¡(s = 1 ∧ 0.875 ≤ x2 ≤ 0.9) → (s0 = 2 ∧ (x0 1 − 1.2)2 + (x0 2 − 1.8)2 ≤ 0.01¢
% ∨ ¡(s = 2∧1.1 ≤ x2 ≤ 1.125) → (s0 = 1∧(x0 1 −0.7)2 +(x0 2 −0.7)2 ≤ 0.01)¢
% Init: s = 1 ∧ (x1 − 0.8)2 + (x2 − 0.2)2 ≤ 0.01
% Unsafe: ¡s = 1 ∧ x1 > 0.8 ∧ x2 > 0.8 ∧ x1 <= 0.9 ∧ x2 ≤ 0.9¢
% State space: (1,[0.1,0.9] × [0.1,0.9]) ∪ (2,[1.1,1.9] × [1.1,1.9])

% independent variables
syms x1 x2;
vars = [x1 x2];

% Constructing the vector field dx/dt = f
f = [-x1+x1*x2;
     x2-x1*x2];

eps = [0.001,0.001];

% Constructing the theta constraint
theta1 = 5*(x1-0.6);
theta2 = 5*(x2-0.1);
g_theta = [theta1,theta2];

% Constructing the psy constraint
psy1 = 5/4*x1-1/8;
psy2 = 5/4*x2-1/8;
g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = 5*(x1-0.5);
zeta2 = 5*(x2-0.5);
g_zeta = [zeta1,zeta2];

end
