function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4EmsoftC14_2Problem()

%====
%Init1=(x1+1)*(x1+2);
%Init2=(x2+2)*(x2+1);
%Init3=(x3+2)*(x3+1);
%Unsafe1=(x1-1)*(x1-2);
%Unsafe2=(x2-2)*(x2-1);
%Unsafe3=(x3-2)*(x3-1);
%Inv1=((x1-2)*(x1+2));
%Inv2=((x2-2)*(x2+2));
%Inv3=((x3-2)*(x3+2));
%B=165.1787867-115.3431*x1+16.1618*x2-210.1866*x3+51.0304*x1^2-23.1069*x1*x2+4.9340*x1*x3+9.9151*x2^2+31.0738*x2*x3-11.6747*x3^2;
%====

% independent variables
syms x1 x2 x3;
vars = [x1, x2, x3];

% Constructing the vector field dx/dt = f
f =  [x1^2+x1*x2-x1*x3; 2*x1*x2+x2^2; x2*x3-2*x3^2];

eps = [0.0001, 0.0001];

% Constructing the theta constraint
theta1 = x1+2;
theta2 = x2+2;
theta3 = x3+2;
g_theta = [theta1, theta2, theta3];

% Constructing the psy constraint
psy1 = (x1+2)/4;
psy2 = (x2+2)/4;
psy3 = (x3+2)/4;
g_psy = [psy1, psy2, psy3];

% Constructing the zeta constraint
zeta1 = x1-1;
zeta2 = x2-1;
zeta3 = x3-1;
g_zeta = [zeta1, zeta2, zeta3];

end
