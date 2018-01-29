function [lpVer, solveResVer, resNorms] = runLp4HybridEx1WithIterations1()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

[vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards] = getLp4HybridEx1Problem();
T = vars(1);
t = vars(2);

% Set the degree
degree = 3;
pLambdaDegree = 1;
pReDegree = 1;

maxIterations = 20;

initPhys = [(5809433765776575*T^3)/75557863725914323419136 + (6773964550827919*T^2*t)/151115727451828646838272 - (7511444873090635*T^2)/590295810358705651712 - (6826529722291971*T*t^2)/1208925819614629174706176 - (3548774222853607*T*t)/590295810358705651712 + (3098470933514395*T)/4611686018427387904 + (570205939113153*t^3)/604462909807314587353088 - (5131637099961173*t^2)/18889465931478580854784 + (1270854446229395*t)/9223372036854775808 - 1786964784349139/72057594037927936,...
    (1860778171644435*T)/1152921504606846976 - (8393104122009939*t)/2305843009213693952 - (5433625772041543*T*t)/147573952589676412928 - (1604740725452595*T*t^2)/2361183241434822606848 - (4929237851062029*T^2*t)/1180591620717411303424 - (604418235700997*T^2)/18446744073709551616 - (8834635609939629*T^3)/1180591620717411303424 + (5181608096556247*t^2)/73786976294838206464 - (6531363188139113*t^3)/18889465931478580854784 - 4482866174013461/144115188075855872,...
    (6272931833902507*T)/4611686018427387904 - (1741427425521125*t)/36893488147419103232 - (8922114595439145*T*t)/4722366482869645213696 + (2847605919070353*T*t^2)/302231454903657293676544 - (3490611404985401*T^2*t)/75557863725914323419136 - (2498013317334515*T^2)/147573952589676412928 - (1530364546513995*T^3)/73786976294838206464 + (2229552391879787*t^2)/4722366482869645213696 - 3344134566277361/36028797018963968];



import lp4.runAndVerifyHLPWithIterations1
[lpVer, solveResVer, resNorms] = runAndVerifyHLPWithIterations1(...
    vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards, degree, pLambdaDegree, pReDegree,...
    maxIterations, initPhys);

warning('on')
echo off;

end



% --------------------------------------------------------------
% Iteration 1 :
% --------------------------------------------------------------
% constraint theta is processed: 2018-01-30 01:18:10
% constraint psy is processed: 2018-01-30 01:18:32
% constraint psy is processed: 2018-01-30 01:18:53
% constraint psy is processed: 2018-01-30 01:19:14
% constraint guard is processed: 2018-01-30 01:19:35
% constraint guard is processed: 2018-01-30 01:19:57
% constraint guard is processed: 2018-01-30 01:20:19
% constraint guard is processed: 2018-01-30 01:20:41
% constraint re is processed: 2018-01-30 01:20:45
% constraint re is processed: 2018-01-30 01:20:50
% constraint re is processed: 2018-01-30 01:20:54
% constraint re is processed: 2018-01-30 01:20:58
% constraint zeta is processed: 2018-01-30 01:21:13
% 
% Optimal solution found.
% 
% --------------------------------------------------------------
% The parameter setting:
% ; lambda degree: 1; re degree: 1; eps1: 0.1; eps2: 0.1
% --------------------------------------------------------------
% The coefficients of lambda1 is:
%    3.438836645327228   0.125306851070823   0.022448916535679
% 
% The function lambda1 is:
% (2257327551158111*T)/18014398509481984 + (6470459657277231*t)/288230376151711744 + 7743571717241855/2251799813685248
%  
% The coefficients of lambda2 is:
%    2.979930734713809   0.190516237375207                   0
% 
% The function lambda2 is:
% (6864070845208083*T)/36028797018963968 + 1677551868305875/562949953421312
%  
% The coefficients of lambda3 is:
%    1.822645333297169   0.000037173232411                   0
% 
% The function lambda3 is:
% (5485800837455593*T)/147573952589676412928 + 8208464843865703/4503599627370496
%  
% --------------------------------------------------------------
% The coefficients of re1 is:
%    0.759190166559266                   0   0.007377809860181
% 
% The function re1 is:
% (8506035644703529*t)/1152921504606846976 + 6838177102439307/9007199254740992
%  
% The coefficients of re2 is:
%    0.018516769465381  -0.000019646863763   0.000048362780940
% 
% The function re2 is:
% (1784271685395123*t)/36893488147419103232 - (5798730683115349*T)/295147905179352825856 + 1334273857030325/72057594037927936
%  
% The coefficients of re3 is:
%    2.499146049338394   0.045727605577685                   0
% 
% The function re3 is:
% (6590042478086637*T)/144115188075855872 + 1406894152068105/562949953421312
%  
% The coefficients of re4 is:
%    0.373305064036888                   0                   0
% 
% The function re4 is:
% 105076034205753/281474976710656
%  
% --------------------------------------------------------------
% The computation time is:
%    0.021670355232166
% 
% --------------------------------------------------------------
% The rou is: 0.0012863
% --------------------------------------------------------------
% constraint theta is processed: 2018-01-30 01:21:40
% constraint psy is processed: 2018-01-30 01:21:59
% constraint psy is processed: 2018-01-30 01:22:19
% constraint psy is processed: 2018-01-30 01:22:39
% constraint guard is processed: 2018-01-30 01:22:58
% constraint guard is processed: 2018-01-30 01:23:18
% constraint guard is processed: 2018-01-30 01:23:38
% constraint guard is processed: 2018-01-30 01:23:52
% constraint re is processed: 2018-01-30 01:23:56
% constraint re is processed: 2018-01-30 01:24:00
% constraint re is processed: 2018-01-30 01:24:05
% constraint re is processed: 2018-01-30 01:24:06
% constraint zeta is processed: 2018-01-30 01:24:21
% 错误使用 linprog (line 336)
% LINPROG encountered an internal error that caused the solution to lose feasibility. We are sorry for the
% inconvenience.
% 
% Please contact technical support for assistance with your problem, quoting the code "-2_4".
% 
% 出错 lp4.HybridLinearProgramVerificationWithGivenLambdaAndRe/solve (line 439)
%             [x, fval, flag, ~] = linprog(this.linprogF, Aie, bie, Aeq, beq);
% 
% 出错 lp4.runAndVerifyHLPWithIterations1 (line 42)
%     [lpVer, solveResVer, resNorms] = lpVer.solve();
% 
% 出错 runLp4HybridEx1WithIterations1 (line 29)
%     vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards, degree, pLambdaDegree,
%     pReDegree,...
