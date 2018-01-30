function [lpVer, solveResVer, resNorms] = runLp4HybridEx1WithRounds()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

[vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards] = getLp4HybridEx1Problem();

% Set the degree
degree = 1;
pLambdaDegree = 1;
pReDegree = 1;

randomStartCount = 100;
randomStartBeginIndex = 1;
maxIterations = 20;

randomSeed = 0;
rng(randomSeed);
randomSeeds = randi(10000000, 2, randomStartCount, 'int32');

lambdaStart = -1;
lambdaEnd = 1;
resStart = 0;
resEnd = 1;

for i = randomStartBeginIndex : randomStartCount
    lp4.Lp4Config.displayDelimiterLine();
    disp(['Begin with random start ', num2str(i), ' :']);
    lp4.Lp4Config.displayDelimiterLine();
    
    import lp4util.generateRamdomExprs
    initLambdas = generateRamdomExprs(vars, pLambdaDegree, stateNum, lambdaStart, lambdaEnd, randomSeeds(1, i));
    initRes = generateRamdomExprs(vars, degree, length(guards), resStart, resEnd, randomSeeds(2, i));
    
    try
        [lpVer, solveResVer, resNorms] = lp4.runAndVerifyHLPWithIterations2(...
            vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards, degree, pLambdaDegree, pReDegree,...
            maxIterations, initLambdas, initRes);
        
        if solveResVer.hasSolutionWithRou()
            return;
        end
    catch err
        % continue if an error occurs in one iteration
        disp(err);
    end
end

warning('on')
echo off;

end



% degree = 1;
% pLambdaDegree = 1;
% pReDegree = 1;
% 
% rou 5.658246025804554e-06
% 
% 
% --------------------------------------------------------------
% The coefficients of lambda1 is:
%    1.0e+02 *
% 
%    7.432734789827285  -0.029128404484393  -0.001968549722421
% 
% The function lambda1 is:
% 6537902662072243/8796093022208 - (221638994908859*t)/1125899906842624 - (6559133579090519*T)/2251799813685248
%  
% The coefficients of lambda2 is:
%    1.0e+03 *
% 
%    4.917667308774452  -0.048650576370670  -0.530306607985331
% 
% The function lambda2 is:
% 5407032387531419/1099511627776 - (4664626254130565*t)/8796093022208 - (1711739981361793*T)/35184372088832
%  
% The coefficients of lambda3 is:
%    1.0e+02 *
% 
%    2.051167093887285  -0.000926840206275  -0.981485431746266
% 
% The function lambda3 is:
% 3608451312384921/17592186044416 - (431661857879107*t)/4398046511104 - (417411720761115*T)/4503599627370496
%  
% --------------------------------------------------------------
% The coefficients of re1 is:
%    1.0e+02 *
% 
%    8.164019663658847  -0.102033324290714  -0.004255764720636
% 
% The function re1 is:
% 897643454958481/1099511627776 - (3833252082006455*t)/9007199254740992 - (1435991378921983*T)/140737488355328
%  
% The coefficients of re2 is:
%    1.0e+02 *
% 
%    1.790641499352800  -0.001607742997633  -0.003558930124422
% 
% The function re2 is:
% 3150129839546647/17592186044416 - (100174977388657*t)/281474976710656 - (2896252306019613*T)/18014398509481984
%  
% The coefficients of re3 is:
%    1.0e+02 *
% 
%    1.700159972538914  -0.000905607866830  -0.003570661370910
% 
% The function re3 is:
% 5981906108434755/35184372088832 - (6432331687798537*t)/18014398509481984 - (1631398100639305*T)/18014398509481984
%  
% The coefficients of re4 is:
%    0.677285320643057   0.000000066810689  -0.000616534220059
% 
% The function re4 is:
% (154054960997*T)/2305843009213693952 - (1421631121263547*t)/2305843009213693952 + 3050221917671577/4503599627370496
%  
% 
% 
% 
% The coefficients of function phy1 is:
%   -0.003098969269797  -0.000001313083145  -0.005597630966492
% 
% The function phy1 is:
% - (6200859831147461*T)/4722366482869645213696 - (6453629116121589*t)/1152921504606846976 - 3572868313264945/1152921504606846976
%  
% The coefficients of function phy2 is:
%    1.0e-03 *
% 
%   -0.021417800404425  -0.000121456044324  -0.317274346346662
% 
% The function phy2 is:
% - (4588479622878917*T)/37778931862957161709568 - (2926339334105175*t)/9223372036854775808 - 3160709461457793/147573952589676412928
%  
% The coefficients of function phy3 is:
%   -0.107245548519106   0.000287580026298  -5.104440677342149
% 
% The function phy3 is:
% (5304915145821303*T)/18446744073709551616 - (2873544641551613*t)/562949953421312 - 1931964049391163/18014398509481984
%  
% 
% 
% Norms
% 
%    1.0e-08 *
% 
%   1 至 4 列
% 
%    0.000000000000305   0.000712035467816   0.744869031713575   0.000015007075534
% 
%   5 至 8 列
% 
%    0.000027394403205   0.018318616994791   0.000907361564893   0.000000630233324
% 
%   9 至 12 列
% 
%    0.000000533022688   0.000000011102230   0.000000022377260   0.000000000092713
% 
%   13 至 14 列
% 
%    0.000000201747005                   0
