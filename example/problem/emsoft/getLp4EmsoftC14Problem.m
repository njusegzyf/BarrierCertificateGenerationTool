function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4EmsoftC14Problem()

%====
%I= (x1-3/2)^2+(x2-1/2)^2+(x3-3/2)^2-1/4;
%U=(x1-1/2)^2+(x2-3/2)^2+(x3-3/2)^2-1/4;
%Inv1=((x1-2)*(x1+2));
%Inv2=((x2-2)*(x2+2));
%Inv3=((x3-2)*(x3+2));
%B=17.9226844-5.4569*x1-0.7915*x2-19.4404*x3+50.2870*x1^2-30.3634*x1*x2+4.8555*x1*x3+12.9877*x2^2-19.7978*x2*x3+0.7835*x3^2;
%====

% independent variables
syms x1 x2 x3;
vars = [x1, x2, x3];

% Constructing the vector field dx/dt = f
f =  [x1^2+x1*x2-x1*x3; 2*x1*x2+x2^2; x2*x3-2*x3^2];

eps = [0.00015, 0.00015];

% Constructing the theta constraint
theta1 = 4*((x1-3/2)^2+(x2-1/2)^2+(x3-3/2)^2-1/4);
theta2 = x1-1;
theta3 = x2;
theta4 = x3-1;
g_theta = [theta1, theta2, theta3, theta4];

% Constructing the psy constraint
psy1 = (x1+2)/4;
psy2 = (x2+2)/4;
psy3 = (x3+2)/4;
g_psy = [psy1, psy2, psy3];

% Constructing the zeta constraint
zeta1 = 4*((x1-1/2)^2+(x2-3/2)^2+(x3-3/2)^2-1/4);
zeta2 = x1;
zeta3 = x2-1;
zeta4 = x3-1;
g_zeta = [zeta1, zeta2, zeta3, zeta4];

end
