function results = verifyLp4EmsoftBiologicalModel2WithDifferentRegions2(combinationsPrefix)
% This function can specify to run just a part of the problems,
% that is, we only run combinations whose prefix matches `combinationsPrefix`
%
% For example, run verifyLp4EmsoftBiologicalModel2WithDifferentRegions2([2, 2, 2]) to run 1/8 of all problems,
% or just run verifyLp4EmsoftBiologicalModel2WithDifferentRegions2() to run all problems.
% 
% Note: In Matlab2017b, the problems can be verified with builtin linprog.
% However, the verification failed with builtin linprog, but succeeded with cvx.

% clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% set dfault argument
if nargin < 1
    combinationsPrefix = [];
end

% since the process prdouce a lot of  text output in command window, we creat a log for the text output
% @see https://www.mathworks.com/help/matlab/ref/diary.html 
logFile = strcat([lp4.Lp4Config.LOG_DIR, 'BiologicalModel2WithDifferentRegions2-', num2str(combinationsPrefix), '.txt']);
lp4.Lp4Config.beginLogging(logFile, true);



% get the problem
[vars, f, eps, theta, ~, zeta] = getLp4EmsoftBiologicalModel2();

% get the test problem for ts
[varsTest, fTest, epsTest, ~, ~, ~] = getLp4EmsoftGeneralSubItemsTestProblem(9);

x1 = vars(1);
x2 = vars(2);
x3 = vars(3);
x4 = vars(4);
x5 = vars(5);
x6 = vars(6);
x7 = vars(7);
x8 = vars(8);
x9 = vars(9);

t1 = -0.4105e-1*x1-0.8081e-1*x2-0.6141e-1*x3-0.7016e-1*x4-0.9917e-2*x5-0.3012e-1*x6-0.115e-1*x7-0.7954e-2*x8-0.6927e-2*x9+.2-0.5727e-1*x1*x2-0.1009e-1*x1*x3-.1057*x1*x4-.1554*x1*x5+0.2475e-1*x1*x6-0.3015e-1*x1*x7-0.4277e-1*x1*x8+0.201e-1*x1*x9-0.5481e-1*x2*x3+0.9027e-1*x2*x4-0.1572e-1*x2*x5-0.5454e-2*x2*x6+0.5935e-2*x2*x7+0.1227e-1*x2*x8+0.3387e-2*x2*x9-.1465*x3*x4+0.7915e-1*x3*x5+0.9248e-2*x3*x6-0.8794e-2*x3*x7-0.2567e-1*x3*x8-0.1254e-1*x3*x9+0.4661e-4*x4*x5-0.3312e-1*x4*x6-0.4957e-1*x4*x7+0.5629e-1*x4*x8+0.1284e-1*x4*x9-0.4994e-1*x5*x6-0.3897e-1*x5*x7-0.5854e-2*x5*x8-0.3928e-1*x5*x9+0.4684e-2*x6*x7-0.6267e-2*x6*x8+0.563e-2*x6*x9-0.652e-2*x7*x8+0.3872e-2*x7*x9-0.4024e-2*x8*x9+.2485*x1^2+0.5471e-1*x2^2+.198*x3^2+.1736*x4^2+.102*x5^2+0.5103e-1*x6^2+0.6661e-1*x7^2+0.3871e-1*x8^2+0.3678e-1*x9^2;

t2 = -0.5468e-1*x1-0.5941e-1*x2-0.916e-1*x3-.1761*x4-0.2401e-1*x5+0.1958e-1*x6-0.6482e-1*x7-0.1513e-1*x8-0.6664e-2*x9+.5497-0.5027e-1*x1*x2+.1077*x1*x3-.125*x1*x4-0.6964e-1*x1*x5+0.2012e-1*x1*x6-0.1612e-1*x1*x7-0.4179e-1*x1*x8+0.6631e-2*x1*x9-0.5041e-1*x2*x3+0.9864e-1*x2*x4-0.8405e-2*x2*x5-0.1266e-1*x2*x6-0.4374e-2*x2*x7+0.1661e-1*x2*x8+0.4552e-2*x2*x9-.1261*x3*x4+0.134e-1*x3*x5+0.9303e-2*x3*x6-0.1179e-1*x3*x7-0.3964e-1*x3*x8-0.2456e-3*x3*x9+0.2498e-1*x4*x5-0.2873e-1*x4*x6-0.8832e-1*x4*x7+0.6213e-1*x4*x8+0.1028e-1*x4*x9-0.4948e-1*x5*x6-0.8193e-1*x5*x7-0.1396e-1*x5*x8-0.4032e-1*x5*x9-0.2326e-1*x6*x7-0.1391e-1*x6*x8+0.4225e-3*x6*x9+0.5716e-3*x7*x8+0.5171e-2*x7*x9-0.2752e-2*x8*x9+.1503*x1^2+.1738*x2^2+.1901*x3^2+.2352*x4^2+.1244*x5^2+0.7041e-1*x6^2+.127*x7^2+0.8006e-1*x8^2+0.7896e-1*x9^2;

t3 = -0.4251e-1*x1-0.6791e-1*x2-0.5742e-1*x3-0.7355e-1*x4-0.1494e-1*x5+0.764e-2*x6-0.4868e-2*x7-0.8122e-2*x8-0.6796e-2*x9+.1564-0.5197e-1*x1*x2+0.9937e-1*x1*x3-.1117*x1*x4-0.4775e-1*x1*x5+0.2676e-1*x1*x6-0.1921e-1*x1*x7-0.3518e-1*x1*x8+0.5821e-2*x1*x9-0.6383e-1*x2*x3+0.9536e-1*x2*x4-0.1535e-1*x2*x5-0.9895e-2*x2*x6-0.2189e-2*x2*x7+0.1296e-1*x2*x8+0.2964e-2*x2*x9-.1532*x3*x4+0.3634e-1*x3*x5-0.9179e-2*x3*x6-0.1243e-1*x3*x7-0.3414e-1*x3*x8-0.181e-2*x3*x9-0.4067e-2*x4*x5-0.3123e-1*x4*x6-0.4116e-1*x4*x7+0.5302e-1*x4*x8+0.1283e-1*x4*x9-0.4551e-1*x5*x6-0.378e-1*x5*x7-0.7297e-2*x5*x8-0.2769e-1*x5*x9+0.4619e-2*x6*x7-0.9466e-2*x6*x8+0.5761e-2*x6*x9-0.3414e-2*x7*x8+0.3302e-2*x7*x9-0.202e-2*x8*x9+.12*x1^2+0.5026e-1*x2^2+.2302*x3^2+.1632*x4^2+0.7081e-1*x5^2+0.3975e-1*x6^2+0.5674e-1*x7^2+0.3266e-1*x8^2+0.3001e-1*x9^2;

t4 = -0.1829e-1*x1-0.8866e-1*x2-0.3888e-1*x3-.101*x4-0.1311e-1*x5+0.9231e-2*x6-0.7562e-2*x7-0.1389e-1*x8-0.9402e-2*x9-0.5658e-1*x1*x2+0.9595e-1*x1*x3-.1201*x1*x4-0.4986e-1*x1*x5+0.2365e-1*x1*x6-0.1423e-1*x1*x7-0.3526e-1*x1*x8+0.4848e-2*x1*x9-0.5454e-1*x2*x3+.1105*x2*x4-0.1277e-1*x2*x5-0.1209e-1*x2*x6-0.8841e-2*x2*x7+0.1666e-1*x2*x8+0.3076e-2*x2*x9-.1496*x3*x4+0.1716e-1*x3*x5+0.8221e-2*x3*x6-0.1356e-3*x3*x7-0.3318e-1*x3*x8-0.1955e-2*x3*x9-0.6768e-2*x4*x5-0.4722e-1*x4*x6-0.587e-1*x4*x7+0.6346e-1*x4*x8+0.1843e-1*x4*x9-0.4411e-1*x5*x6-0.3678e-1*x5*x7-0.8449e-2*x5*x8-0.2811e-1*x5*x9+0.9517e-2*x6*x7-0.1125e-1*x6*x8+0.4639e-2*x6*x9-0.6487e-2*x7*x8+0.1872e-2*x7*x9-0.9355e-3*x8*x9+.1149*x1^2+0.5584e-1*x2^2+.1491*x3^2+.2382*x4^2+0.7146e-1*x5^2+0.4092e-1*x6^2+0.6e-1*x7^2+0.3482e-1*x8^2+0.3213e-1*x9^2+.1716;

t5 = -0.2748e-1*x1-0.9823e-1*x2-0.1508e-1*x3-0.9888e-1*x4+0.124e-1*x5+.1252*x6+0.5161e-2*x7+0.576e-1*x8+0.9084e-2*x9-0.5731e-1*x1*x2+.134*x1*x3-.1136*x1*x4-0.4268e-1*x1*x5+0.2417e-1*x1*x6-0.1965e-1*x1*x7-0.3357e-1*x1*x8+0.5787e-2*x1*x9-0.5398e-1*x2*x3+.121*x2*x4-0.4048e-2*x2*x5-0.1283e-1*x2*x6-0.1145e-1*x2*x7+0.2073e-1*x2*x8+0.492e-2*x2*x9-.1259*x3*x4+0.5185e-1*x3*x5+0.1862e-1*x3*x6-0.8806e-2*x3*x7-0.1378e-1*x3*x8-0.4562e-2*x3*x9+0.5995e-2*x4*x5-0.2923e-1*x4*x6-0.4259e-1*x4*x7+0.672e-1*x4*x8+0.196e-1*x4*x9-0.3809e-1*x5*x6-0.456e-1*x5*x7+0.1167e-1*x5*x8-0.3725e-1*x5*x9+0.1476e-1*x6*x7+0.1116e-1*x6*x8+0.2008e-1*x6*x9+0.1174e-1*x7*x8+0.9659e-2*x7*x9+0.6629e-2*x8*x9+.1306*x1^2+0.7432e-1*x2^2+.196*x3^2+.1898*x4^2+.1636*x5^2+0.6906e-1*x6^2+0.8035e-1*x7^2+0.6204e-1*x8^2+0.5676e-1*x9^2+.3485;

t6 = -.2205*x1-.651*x2-.1132*x3-.1443*x4-0.9088e-1*x5-0.7105e-1*x6-0.5853e-1*x7-0.6263e-1*x8-0.2471e-1*x9-0.3612e-1*x1*x2-0.9382e-1*x1*x3-0.9492e-1*x1*x4-.2645*x1*x5+0.6884e-1*x1*x6-0.3579e-1*x1*x7-0.4276e-1*x1*x8+0.325e-1*x1*x9-0.6821e-1*x2*x3+.1347*x2*x4-0.1355e-1*x2*x5-0.1346e-1*x2*x6-0.1042e-1*x2*x7+0.1367e-1*x2*x8+0.2879e-2*x2*x9-.1738*x3*x4+.2108*x3*x5-0.4336e-1*x3*x6-0.4217e-2*x3*x7-0.2702e-1*x3*x8-0.3639e-1*x3*x9-0.6773e-2*x4*x5-0.5917e-1*x4*x6-0.8911e-1*x4*x7+0.6937e-1*x4*x8+0.1708e-1*x4*x9-.1066*x5*x6-0.5809e-1*x5*x7-0.101e-1*x5*x8-0.7111e-1*x5*x9+0.3087e-1*x6*x7-0.2194e-1*x6*x8+0.1676e-1*x6*x9-0.4054e-2*x7*x8+0.8176e-2*x7*x9-0.4548e-2*x8*x9+.3822*x1^2+.1677*x2^2+.3376*x3^2+.2352*x4^2+.1886*x5^2+.1853*x6^2+.1226*x7^2+0.7065e-1*x8^2+0.7063e-1*x9^2+1.246;

t7 = -0.2855e-1*x1-.1494*x2-0.5272e-1*x3-.1266*x4-0.2941e-1*x5+0.4639e-1*x6+0.8627e-3*x7-0.436e-2*x8-0.5805e-2*x9-0.6345e-1*x1*x2+.1078*x1*x3-.1165*x1*x4-0.5353e-1*x1*x5+0.2367e-1*x1*x6-0.2228e-1*x1*x7-0.3625e-1*x1*x8+0.6584e-2*x1*x9-0.7037e-1*x2*x3+.123*x2*x4-0.1549e-1*x2*x5-0.3178e-1*x2*x6-0.1707e-1*x2*x7+0.1464e-2*x2*x8+0.1289e-2*x2*x9-.1396*x3*x4+0.2707e-1*x3*x5+0.7915e-2*x3*x6-0.7204e-2*x3*x7-0.2341e-1*x3*x8-0.1289e-2*x3*x9+0.323e-3*x4*x5-0.4445e-1*x4*x6-0.5325e-1*x4*x7+0.5885e-1*x4*x8+0.1342e-1*x4*x9-0.5397e-1*x5*x6-0.4806e-1*x5*x7-0.5521e-2*x5*x8-0.3732e-1*x5*x9+0.4263e-2*x6*x7-0.4964e-2*x6*x8+0.5473e-2*x6*x9-0.7364e-2*x7*x8+0.3818e-3*x7*x9-0.5158e-2*x8*x9+.1323*x1^2+0.8317e-1*x2^2+.1752*x3^2+.1888*x4^2+0.8751e-1*x5^2+0.5576e-1*x6^2+.1295*x7^2+0.4673e-1*x8^2+0.4443e-1*x9^2+.27;

t8 = -0.7322e-1*x1-.1111*x2-0.9771e-1*x3-.1444*x4-0.671e-1*x5-0.7584e-1*x6-0.3818e-1*x7-0.1084e-1*x8-0.4522e-2*x9-0.5167e-1*x1*x2+0.8965e-1*x1*x3-.1077*x1*x4-0.7657e-1*x1*x5+0.2118e-1*x1*x6-0.3023e-1*x1*x7-0.3273e-1*x1*x8+0.6553e-2*x1*x9-0.5304e-1*x2*x3+0.9126e-1*x2*x4-0.2582e-1*x2*x5-0.1097e-1*x2*x6+0.9238e-2*x2*x7+0.1703e-1*x2*x8+0.5887e-2*x2*x9-.1267*x3*x4+0.8396e-2*x3*x5+0.101e-1*x3*x6-0.7431e-2*x3*x7-0.3667e-1*x3*x8+0.3566e-2*x3*x9-0.6793e-3*x4*x5-0.3048e-1*x4*x6-0.3869e-1*x4*x7+0.5133e-1*x4*x8+0.1985e-1*x4*x9-0.633e-1*x5*x6-0.3448e-1*x5*x7-0.2191e-1*x5*x8-0.2918e-1*x5*x9-0.1338e-1*x6*x7-0.2509e-2*x6*x8+0.2907e-2*x6*x9-0.1792e-1*x7*x8+0.1016e-2*x7*x9-0.5874e-2*x8*x9+.1612*x1^2+0.8855e-1*x2^2+.1868*x3^2+.2076*x4^2+.1192*x5^2+0.8511e-1*x6^2+.1037*x7^2+.1595*x8^2+0.7292e-1*x9^2+.5776;

t9 = -0.6047e-1*x1-.1122*x2-0.8302e-1*x3-.1356*x4-0.3278e-1*x5+0.273e-1*x6-0.2965e-1*x7-0.1311e-1*x8-0.4877e-2*x9-0.5373e-1*x1*x2+.1069*x1*x3-.1076*x1*x4-0.591e-1*x1*x5+0.1876e-1*x1*x6-0.2452e-1*x1*x7-0.3878e-1*x1*x8+0.3252e-2*x1*x9-0.5214e-1*x2*x3+.1029*x2*x4-0.1588e-1*x2*x5-0.1631e-1*x2*x6-0.2393e-2*x2*x7+0.128e-1*x2*x8+0.4855e-2*x2*x9-.1254*x3*x4+0.106e-1*x3*x5+0.5424e-2*x3*x6-0.7674e-2*x3*x7-0.3645e-1*x3*x8+0.4552e-2*x3*x9-0.6577e-2*x4*x5-0.3422e-1*x4*x6-0.3974e-1*x4*x7+0.5908e-1*x4*x8+0.198e-1*x4*x9-0.5384e-1*x5*x6-0.4895e-1*x5*x7-0.132e-1*x5*x8-0.2606e-1*x5*x9-0.3986e-2*x6*x7-0.1079e-1*x6*x8+0.3282e-2*x6*x9-0.5153e-2*x7*x8-0.7041e-2*x7*x9-0.1772e-2*x8*x9+.1424*x1^2+0.8334e-1*x2^2+.179*x3^2+.1999*x4^2+.106*x5^2+0.6938e-1*x6^2+0.9436e-1*x7^2+0.6909e-1*x8^2+.1411*x9^2+.4235;

ts = [t1, t2, t3, t4, t5, t6, t7, t8, t9];

initPhy = -8.967*x1-6.376*x2-10.73*x3-6.947*x4+1.357*x5-.2111*x6+.6266*x7-.1777*x8-.4161*x9+30.41+.5227*x1*x2+2.492*x1*x3+.3089*x1*x4+1.391*x1*x5-.114*x1*x6-.8628*x1*x7+.2594*x1*x8+.1612*x1*x9+.8764*x2*x3+.9352*x2*x4-.3628*x2*x5-0.4419e-1*x2*x6+.4258*x2*x7-0.8721e-1*x2*x8-0.4147e-2*x2*x9+1.022*x3*x4-0.5612e-1*x3*x5-0.2406e-2*x3*x6-1.058*x3*x7+.197*x3*x8+.2396*x3*x9-.1054*x4*x5+.3233*x4*x6+1.207*x4*x7-.4179*x4*x8-.115*x4*x9+.5629*x5*x6-.3395*x5*x7-.3277*x5*x8+0.6521e-1*x5*x9-0.9829e-1*x6*x7-0.6071e-1*x6*x8-.2566*x6*x9+.4806*x7*x8+.2038*x7*x9+0.2664e-2*x8*x9+.4788*x1^2+.5521*x2^2-.1785*x3^2+.2405*x4^2-.4099*x5^2-.1549*x6^2-0.4051e-1*x7^2-0.2868e-1*x8^2-.1492*x9^2;
pLambdaDegree = 0;

% x = [-K, K] => [-K, 0], [0, K]
K = 2;

varsLen = 9;
% pre allocate struct array
currentSelectedRegionIndexes(varsLen) = struct('pos', 1, 'regionIndex', 1);

results = [];

% regionCombinationsCount = 2^varsLen;
combinationsPrefixLen = length(combinationsPrefix);
    function processRegionCombinations(currentPos)
        if (currentPos == varsLen + 1)
            for i = 1 : combinationsPrefixLen
                % skip a combination if its prefix does not match
                if currentSelectedRegionIndexes(i).regionIndex ~= combinationsPrefix(i)
                    return
                end
            end
            
            % all region index is enumerated, process it
            
            % get psy from current regions
            verstring = char(version);
            if verstring(1) == '9' % for new matlab version like 2017b
                psy = arrayfun(@(x) getRegion(vars(x.pos), x.regionIndex), currentSelectedRegionIndexes);
            else
                psy(varsLen) = x1; % assign the last element to pre allocate memory
                for x = currentSelectedRegionIndexes
                    pos = x.pos;
                    psy(pos) = getRegion(vars(pos), x.regionIndex);
                end
            end
            
            psyAsStr = char(psy);
            lp4.Lp4Config.displayTextBetweenDelimiterLines(['Process psy: ', psyAsStr]);
            
            % verify each psy combination using ts
            for t = ts
                lpVerTest = lp4.LinearProgram4Verification3.createWithRou(varsTest, fTest, epsTest, psy, [], [], pLambdaDegree, t);
                [~, solveResVerTest, resNormsTest] = lpVerTest.solve();
                lp4.Lp4Config.printVerifyWithOnlyThetaResult(solveResVerTest, resNormsTest);
                
                % store result (append to the end of results)
                results = [results struct('psy', psyAsStr, 'ts', ts, 'result', solveResVerTest)];            
            end

        else
            % for each region for the var at `currentPos`, set the region,
            % and then call `processRegionCombinations` with next position to process remain vars
            for i = 1 : 2
                % set region for the var at `currentPos`
                currentSelectedRegionIndexes(currentPos) = struct('pos', currentPos, 'regionIndex', i);
                processRegionCombinations(currentPos + 1);
            end
        end
    end

    function res = getRegion(x, i)
        if i == 1
            res = x/K + 1; % for x = [-K, 0]
        else % i == 2
            res = x/K; % for x = [0, K]
        end
    end

% verify theta and zeta using phy, which is same to all combinations
lpVer = lp4.LinearProgram4Verification3.createWithRou(vars, f, eps, theta, [], zeta, pLambdaDegree, initPhy);
[lpVer, solveResVer, resNorms] = lpVer.solve();

lp4.Lp4Config.printVerifyWithOnlyThetaResult(solveResVer, resNorms);
lp4.Lp4Config.printVerifyWithOnlyZetaResult(solveResVer, resNorms);

% begin process from first psy combination
processRegionCombinations(1);



lp4.Lp4Config.endLogging(logFile, true);

warning('on')

echo off;

end
