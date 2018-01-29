function [lpVer, solveResVer, resNorms] = runLp4HybridEx1WithIterations1_3()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

[vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards] = getLp4HybridEx1Problem();
T = vars(1);
t = vars(2);

% Set the degree
degree = 2;
pLambdaDegree = 1;
pReDegree = 1;

maxIterations = 20;

initPhys = [(1299*T)/9740 - (8657902762668003*T*t)/2361183241434822606848 - 1299/1948,...
    (1299*T)/9740 - (2398308818870159*t)/144115188075855872 - (8657902762670907*T*t)/2361183241434822606848 - 1299/1948,...
    (3637*T)/27272 - 3153064689246259/4503599627370496];



import lp4.runAndVerifyHLPWithIterations1
[lpVer, solveResVer, resNorms] = runAndVerifyHLPWithIterations1(...
    vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards, degree, pLambdaDegree, pReDegree,...
    maxIterations, initPhys);

warning('on')
echo off;

end

% Optimal solution found.
% 
% --------------------------------------------------------------
% The parameter setting:
% degree: 2; eps1: 0.1; eps2: 0.1
% --------------------------------------------------------------
% The coefficients of function phy1 is:
%   1 至 5 列
% 
%    2.719325779659566                   0                   0                   0                   0
% 
%   6 列
% 
%    0.000048928655807
% 
% The function phy1 is:
% (1805148783081117*t^2)/36893488147419103232 + 6123377283986903/2251799813685248
%  
% The coefficients of function phy2 is:
%   1 至 5 列
% 
%   -0.067293323092515   0.024062506974034   0.031041186913244  -0.000000000000000                   0
% 
%   6 列
% 
%   -0.000113240177380
% 
% The function phy2 is:
% (1733886359069771*T)/72057594037927936 + (4473506490099937*t)/144115188075855872 - (1413697157682783*T^2)/20769187434139310514121985316880384 - (8355650283930911*t^2)/73786976294838206464 - 4848994956863559/72057594037927936
%  
% The coefficients of function phy3 is:
%   1 至 5 列
% 
%   -0.097691035112777  -0.000078356377265                   0   0.000000000000000                   0
% 
%   6 列
% 
%                    0
% 
% The function phy3 is:
% (8645256203897683*T^2)/83076749736557242056487941267521536 - (5781680151835247*T)/73786976294838206464 - 7039380949301427/72057594037927936
%  
% --------------------------------------------------------------
% The computation time is:
%    0.013282391359828
% 
% --------------------------------------------------------------
% The rou is: 0.0097857
% --------------------------------------------------------------
% --------------------------------------------------------------
% Iteration 9 :
% --------------------------------------------------------------
% constraint theta is processed: 2018-01-30 01:24:39
% constraint psy is processed: 2018-01-30 01:24:49
% constraint psy is processed: 2018-01-30 01:24:56
% constraint psy is processed: 2018-01-30 01:25:04
% constraint guard is processed: 2018-01-30 01:25:10
% constraint guard is processed: 2018-01-30 01:25:16
% constraint guard is processed: 2018-01-30 01:25:22
% constraint guard is processed: 2018-01-30 01:25:28
% constraint re is processed: 2018-01-30 01:25:30
% constraint re is processed: 2018-01-30 01:25:32
% constraint re is processed: 2018-01-30 01:25:34
% constraint re is processed: 2018-01-30 01:25:36
% constraint zeta is processed: 2018-01-30 01:25:40
% 
% Optimal solution found.
% 
% --------------------------------------------------------------
% The parameter setting:
% ; lambda degree: 1; re degree: 1; eps1: 0.1; eps2: 0.1
% --------------------------------------------------------------
% The coefficients of lambda1 is:
%   -0.015194642928874   0.000000000000000  -0.000107875315658
% 
% The function lambda1 is:
% (6644932866285051*T)/83076749736557242056487941267521536 - (7959793359266645*t)/73786976294838206464 - 8759115293760403/576460752303423488
%  
% The coefficients of lambda2 is:
%      0     0     0
% 
% The function lambda2 is:
% 0
%  
% The coefficients of lambda3 is:
%    0.022506060862957   0.003585423739641   0.400451728267932
% 
% The function lambda3 is:
% (4133712132560379*T)/1152921504606846976 + (56358570443979*t)/140737488355328 + 810866298527939/36028797018963968
%  
% --------------------------------------------------------------
% The coefficients of re1 is:
%   -0.010360207466149                   0  -0.000035965297115
% 
% The function re1 is:
% - (5307541051325415*t)/147573952589676412928 - 2986126494977875/288230376151711744
%  
% The coefficients of re2 is:
%      0     0     0
% 
% The function re2 is:
% 0
%  
% The coefficients of re3 is:
%   -0.019277242834532   0.000151860818885  -0.000100631650509
% 
% The function re3 is:
% (2801337660796025*T)/18446744073709551616 - (7425305210591519*t)/73786976294838206464 - 2778143476682501/144115188075855872
%  
% The coefficients of re4 is:
%    0.088160672182258  -0.000000000000000                   0
% 
% The function re4 is:
% 198520185194375/2251799813685248 - (4812363774998949*T)/332306998946228968225951765070086144
%  
% --------------------------------------------------------------
% The computation time is:
%    0.012631328602593
% 
% --------------------------------------------------------------
% The rou is: 0.0097801
% --------------------------------------------------------------
% constraint theta is processed: 2018-01-30 01:25:50
% constraint psy is processed: 2018-01-30 01:25:55
% constraint psy is processed: 2018-01-30 01:25:56
% constraint psy is processed: 2018-01-30 01:26:01
% constraint guard is processed: 2018-01-30 01:26:06
% constraint guard is processed: 2018-01-30 01:26:09
% constraint guard is processed: 2018-01-30 01:26:14
% constraint guard is processed: 2018-01-30 01:26:20
% constraint re is processed: 2018-01-30 01:26:21
% constraint re is processed: 2018-01-30 01:26:22
% constraint re is processed: 2018-01-30 01:26:23
% constraint re is processed: 2018-01-30 01:26:25
% constraint zeta is processed: 2018-01-30 01:26:28
% 
% Problem is unbounded.
% 
% --------------------------------------------------------------
% The parameter setting:
% degree: 2; eps1: 0.1; eps2: 0.1
% --------------------------------------------------------------
% The problem with degree 2 maybe have no solution.
% --------------------------------------------------------------
% Unable to find lambda and re for next iterations.
% 
% ans = 
% 
%   HybridLinearProgramVerificationWithGivenLambdaAndRe - 属性:
% 
%             indvars: [1×2 sym]
%              degree: 2
%            stateNum: 3
%     thetaStateIndex: 2
%            guardNum: 4
%                phys: [1×3 sym]
%      phyPolynomials: [1×3 lp4util.SymbolicPolynomial]
%                 eps: [0.100000000000000 0.100000000000000]
%             lambdas: [1×3 sym]
%                 res: [1×4 sym]
%                  fs: [2×3 sym]
%               theta: [1×3 sym]
%                psys: [3×2 sym]
%               zetas: [1×1 lp4.UnsafeConstraint]
%              guards: [1×4 lp4.Guard]
%             decvars: [1×273 sym]
%      decvarsIndexes: [1×1 lp4.HybridLinearProgramDecVarIndexes]
%               exprs: [1×14 Constraint]
%            linprogF: [1×273 double]
%         isAttachRou: 1
%              rouVar: [1×1 sym]
