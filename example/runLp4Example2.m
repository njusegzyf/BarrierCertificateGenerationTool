function [lp, solveRes] = runLp4Example2() 

% Example 2

clear; 
echo on;

% disable warning of `Support of character vectors will be removed in a future release.` 
% which is produced by function `monomials`.
warning('off')

% independent variables
syms x1 x2;
vars = [x1 x2];

% Constructing the vector field dx/dt = f
f = [x2;
    -x1+(1/3)*x1^3-x2];

eps = [0.00001,0.00001];

% phy = 151/99+152/99*x1+62/33*x2+106/99*x1*x2+4/9*x1^2
% Constructing the theta constraint
theta1 = 4*(x1-1.5)^2+4*x2^2;
theta2 = x2+1/2;
theta3 = x1-1;

g_theta = [theta1,theta2,theta3];

% Constructing the psy constraint
psy1 = (x1+2)/4;
psy2 = (x2+2)/4;

g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = 6.25*(x1+1)^2+6.25*(x2+1)^2;
zeta2 = x1+7/5;
zeta3 = x2+7/5;
g_zeta = [zeta1,zeta2,zeta3];

import lp4.LinearProgram4_v2
lp = LinearProgram4_v2(vars);

lp.f = f;
lp.eps = eps;

% Set the degree of phy
degree = 3;
pLambdaDegree = 3;
lp = lp.setDegreeAndInit(degree, pLambdaDegree);

% Note: degree 4 is OK.

lp = lp.setThetaConstraint(g_theta);
lp = lp.setPsyConstraint(g_psy);
lp = lp.setZetaConstraint(g_zeta);
lp = lp.generateEqsForConstraint1To3();

import lp4util.Partition
lp.pPartitions = repmat(Partition(-1, 1), 1000, 1);
lp.pLambdaPartitions = repmat(Partition(-1, 1), 1000, 1);
lp = lp.setWConstraint();

lp = lp.setDevVarsConstraint();

lp = lp.setLinprogF();

% see http://blog.sina.com.cn/s/blog_68b0c65f0100mq5m.html
% options = optimset(��MaxIter��, 2000);

% solve the lp problem
[lp, solveRes] = lp.solve();
% [lp, solveRes] = lp.solve(options);

% [lp, solveRes] = lp.solve1And3();

% --------------------------------------------------------------
% The parameter setting:
% degree: 4; r: 0; M: 0; eps1: 0.0001; eps2: 0.0001
% --------------------------------------------------------------
% The coefficients of function phy is:
%    1.0e+03 *
% 
%    -2.1335
%    -0.0647
%     0.0590
%     0.2532
%    -0.0026
%     0.2172
%    -0.0276
%     0.0010
%     0.0004
%    -0.0257
%    -0.0001
%    -0.0001
%     0.0000
%    -0.0001
%    -0.0000
% 
% --------------------------------------------------------------
% The function phy is:
% - (12205204661007*x1^4)/140737488355328 - (4000727484283*x1^3*x2)/70368744177664 - (3886697156809601*x1^3)/140737488355328 + (3374997630337*x1^2*x2^2)/140737488355328 + (141737353587345*x1^2*x2)/140737488355328 + (278401594010247*x1^2)/1099511627776 - (13981847497209*x1*x2^3)/140737488355328 + (57641084288791*x1*x2^2)/140737488355328 - (372890363379441*x1*x2)/140737488355328 - (2275945995120445*x1)/35184372088832 - (5710497927511*x2^4)/281474976710656 - (7233419319987091*x2^3)/281474976710656 + (1910264361808549*x2^2)/8796093022208 + (2077126693041003*x2)/35184372088832 - 4691662008282723/2199023255552
%  
% --------------------------------------------------------------
% The computation time is:
%     0.0249
% 
% --------------------------------------------------------------

% --------------------------------------------------------------
% The parameter setting:
% degree: 4; r: 0; M: 0; eps1: 0.01; eps2: 0.01
% --------------------------------------------------------------
% The coefficients of function phy is:
%    1.0e+03 *
% 
%    -2.1336
%    -0.0647
%     0.0590
%     0.2532
%    -0.0026
%     0.2172
%    -0.0276
%     0.0010
%     0.0004
%    -0.0257
%    -0.0001
%    -0.0001
%     0.0000
%    -0.0001
%    -0.0000
% 
% --------------------------------------------------------------
% The function phy is:
% - (12205373839385*x1^4)/140737488355328 - (4000800761423*x1^3*x2)/70368744177664 - (7773722672459515*x1^3)/281474976710656 + (3375044523105*x1^2*x2^2)/140737488355328 + (141739594308091*x1^2*x2)/140737488355328 + (2227305669736395*x1^2)/8796093022208 - (1747757236599*x1*x2^3)/17592186044416 + (7205269964973*x1*x2^2)/17592186044416 - (372898337587405*x1*x2)/140737488355328 - (4551983968877367*x1)/70368744177664 - (5710465947261*x2^4)/281474976710656 - (7233739129416035*x2^3)/281474976710656 + (7462295001795*x2^2)/34359738368 + (4154318668019867*x2)/70368744177664 - 4691877763811201/2199023255552
%  
% --------------------------------------------------------------
% The computation time is:
%     0.0250
% 
% --------------------------------------------------------------

% --------------------------------------------------------------
% The parameter setting:
% degree: 4; r: -1; M: 0; eps1: 0.01; eps2: 0.01
% --------------------------------------------------------------
% The coefficients of function phy is:
%    1.0e+03 *
% 
%    -1.7477
%    -0.3068
%     0.2563
%     0.2477
%    -0.0258
%     0.1290
%    -0.0134
%    -0.0028
%     0.0183
%    -0.0183
%    -0.0013
%     0.0003
%     0.0000
%    -0.0022
%     0.0001
% 
% --------------------------------------------------------------
% The function phy is:
% - (91196474856733*x1^4)/70368744177664 + (21082238871629*x1^3*x2)/70368744177664 - (235976748503395*x1^3)/17592186044416 + (3112517737129*x1^2*x2^2)/140737488355328 - (393249026201687*x1^2*x2)/140737488355328 + (1089236740472215*x1^2)/4398046511104 - (19524536876495*x1*x2^3)/8796093022208 + (1288729701754775*x1*x2^2)/70368744177664 - (906004446802527*x1*x2)/35184372088832 - (674623777254217*x1)/2199023255552 + (12895035011651*x2^4)/140737488355328 - (2579600017717517*x2^3)/140737488355328 + (4538725214328039*x2^2)/35184372088832 + (281841654922081*x2)/1099511627776 - 7686253856605263/4398046511104
%  
% --------------------------------------------------------------
% The computation time is:
%     0.0266
% 
% --------------------------------------------------------------

% --------------------------------------------------------------
% The parameter setting:
% degree: 4; r: 1; M: 0; eps1: 0.01; eps2: 0.01
% --------------------------------------------------------------
% The coefficients of function phy is:
%    1.0e+03 *
% 
%    -1.5098
%     0.0369
%    -0.0622
%     0.1838
%     0.0084
%     0.1862
%    -0.0300
%     0.0036
%    -0.0115
%    -0.0219
%     0.0010
%    -0.0004
%    -0.0000
%     0.0014
%     0.0006
% 
% --------------------------------------------------------------
% The function phy is:
% (136865841786809*x1^4)/140737488355328 - (53233726009715*x1^3*x2)/140737488355328 - (8434546009386931*x1^3)/281474976710656 - (2302117655549*x1^2*x2^2)/70368744177664 + (252390933583385*x1^2*x2)/70368744177664 + (6466364674180581*x1^2)/35184372088832 + (197124203941161*x1*x2^3)/140737488355328 - (1611704810074159*x1*x2^2)/140737488355328 + (1181244140810231*x1*x2)/140737488355328 + (5196648642626241*x1)/140737488355328 + (41314830665085*x2^4)/70368744177664 - (1538203060113773*x2^3)/70368744177664 + (6552744574288613*x2^2)/35184372088832 - (8756596348656219*x2)/140737488355328 - 6640143366182191/4398046511104
%  
% --------------------------------------------------------------
% The computation time is:
%     0.0277
% 
% --------------------------------------------------------------

warning('on')

echo off;
end
