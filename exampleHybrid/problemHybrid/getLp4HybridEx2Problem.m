function [vars, stateNum, fs, eps, thetaStateIndex, g_theta, g_psys, g_zeta, g_guards] = getLp4HybridEx2Problem()

% independent variables
syms x y;
vars = [x y];

stateNum = 2;

% Constructing the vector field dx/dt = f
f1 = [y^2+10*y+25; 2*x*y+10*x-40*y-200];
f2 = [-y^2-10*y-25; 8*x*y+40*x-160*y-800];
fs = [f1, f2];
% Note: use `fs(:, 1)` to get f1

import lp4.Lp4Config
% eps = [Lp4Config.DEFAULT_EPS, Lp4Config.DEFAULT_EPS];
eps = [0.01, 0.01];


% Constructing the theta constraint
theta1 = (x-7)/4;
theta2 = (y-18)/4;
theta3 = 1/4*((x-9)^2+(y-20)^2);
g_theta = [theta1, theta2, theta3];
thetaStateIndex = 1;



% Constructing the psy constraint
psy11 = (x-5)/30;
psy12 = y/60;
g_psy1 = [psy11, psy12];

psy21 = (x-5)/30;
psy22 = y/60;
g_psy2 = [psy21, psy22];

g_psys = [g_psy1; g_psy2];



% Constructing the Guard constraint
import lp4.Guard

guard121=x-35;
guard122=35-x;
guard123=y/60;
guard12 = Guard(1, 2, [guard121, guard122, guard123], [], []);

guard211=x-5;
guard212=5-x;
guard213=y/60;
guard21 = Guard(2, 1, [guard211, guard212, guard213], [], []);

g_guards = [guard12, guard21]; % [guard12, guard21, guard23, guard32];



% Constructing the zeta constraint
import lp4.UnsafeConstraint
zeta11 = (x-5)/30;
zeta12 = (y-48)/12;
g_zeta = [UnsafeConstraint(1, [zeta11, zeta12]),UnsafeConstraint(2, [zeta11, zeta12])];

end

% --------------------------------------------------------------
% Begin with random start 21 :
% --------------------------------------------------------------
% --------------------------------------------------------------
% Iteration 1 :
% --------------------------------------------------------------
% constraint theta is processed: 2018-01-30 15:07:23
% constraint psy is processed: 2018-01-30 15:07:35
% constraint psy is processed: 2018-01-30 15:07:48
% constraint guard is processed: 2018-01-30 15:08:08
% constraint guard is processed: 2018-01-30 15:08:28
% constraint re is processed: 2018-01-30 15:08:35
% constraint re is processed: 2018-01-30 15:08:43
% constraint zeta is processed: 2018-01-30 15:08:51
% constraint zeta is processed: 2018-01-30 15:08:59
% 
% Optimal solution found.
% 
% --------------------------------------------------------------
% The parameter setting:
% degree: 2; eps1: 0.01; eps2: 0.01
% --------------------------------------------------------------
% The coefficients of function phy1 is:
%    1.0e-03 *
% 
%   1 至 5 列
% 
%   -0.370433625372584  -0.079913144496633  -0.093770752382241   0.002255756952604   0.001269015390544
% 
%   6 列
% 
%   -0.001711894756446
% 
% The function phy1 is:
% (93636808539863*x*y)/73786976294838206464 - (1729765070794385*y)/18446744073709551616 - (5896549298619021*x)/73786976294838206464 + (332890969577445*x^2)/147573952589676412928 - (8084194420039151*y^2)/4722366482869645213696 - 1708323570886115/4611686018427387904
%  
% The coefficients of function phy2 is:
%   1 至 5 列
% 
%   -0.003402012456214   0.000264142515237  -0.000014752428952  -0.000016901100867   0.000001769064518
% 
%   6 列
% 
%   -0.000002996628532
% 
% The function phy2 is:
% (304535586097563*x)/1152921504606846976 - (4354148501446257*y)/295147905179352825856 + (1044271373024283*x*y)/590295810358705651712 - (1247081129027395*x^2)/73786976294838206464 - (7075589070384567*y^2)/2361183241434822606848 - 3922253319709225/1152921504606846976
%  
% --------------------------------------------------------------
% The computation time is:
%    0.028935946960509
% 
% --------------------------------------------------------------
% The rou is: 0.00054711
% --------------------------------------------------------------
% constraint theta is processed: 2018-01-30 15:09:12
% constraint psy is processed: 2018-01-30 15:09:19
% constraint psy is processed: 2018-01-30 15:09:26
% constraint guard is processed: 2018-01-30 15:09:33
% constraint guard is processed: 2018-01-30 15:09:40
% constraint re is processed: 2018-01-30 15:09:42
% constraint re is processed: 2018-01-30 15:09:44
% constraint zeta is processed: 2018-01-30 15:09:48
% constraint zeta is processed: 2018-01-30 15:09:52
% 
% Optimal solution found.
% 
% --------------------------------------------------------------
% The parameter setting:
% ; lambda degree: 1; re degree: 1; eps1: 0.01; eps2: 0.01
% --------------------------------------------------------------
% The coefficients of lambda1 is:
%  -25.394555461694051   3.317696046506840  -0.741292877827920
% 
% The function lambda1 is:
% (3735393669694193*x)/1125899906842624 - (3338486328358221*y)/4503599627370496 - 7147931907157795/281474976710656
%  
% The coefficients of lambda2 is:
%   -0.031548804939597   0.274933375531270  -0.793168254270356
% 
% The function lambda2 is:
% (4952759390377357*x)/18014398509481984 - (55814253974595*y)/70368744177664 - 4546661957438497/144115188075855872
%  
% --------------------------------------------------------------
% The coefficients of re1 is:
%   -0.002188435204734                   0                   0
% 
% The function re1 is:
% -5046188017953363/2305843009213693952
%  
% The coefficients of re2 is:
%   -0.002188435204734                   0                   0
% 
% The function re2 is:
% -5046188017953361/2305843009213693952
%  
% --------------------------------------------------------------
% The computation time is:
%    0.015518444276060
% 
% --------------------------------------------------------------
% The rou is: 0.00054711
% --------------------------------------------------------------
