function [lpVer, solveResVer, resNorms] = runLp4HybridEx1State3Verification()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, g_theta, g_psy, g_zeta] = getLp4HybridEx1State3Problem();
T = vars(1);
t = vars(2);

lambda = (4983313996117259*t)/9007199254740992;
import lp4.verifyWithGivenLambda
[lpVer, solveResVer, resNorms] = verifyWithGivenLambda(vars, f, eps, g_theta, g_psy, g_zeta, lambda, 5);

% The function phy coffs is:
% 
%    1.0e-03 *
% 
%   -0.210000000000000
%    0.000000004812633
%    0.010000000000000
%    0.000000001143596
%                    0
%   -0.000200000000000
%   -0.000000013087820
%    0.000000001143596
%   -0.000000000075077
%    0.000002000000000
%                    0
%                    0
%   -0.000000167219287
%                    0
%   -0.000000010000000
%   -0.000000000000006
%                    0
%    0.000000037155992
%                    0
%                    0
%                    0

% The function phy is:
% - (1024499878798137*T^5)/162259276829213363391578010288128 + (5749611197591933*T^3*t^2)/154742504910672534362390528 - (8100968501419161*T^3)/618970019642690137449562112 - (1617245711072605*T^2*t^2)/9671406556917033397649408 + (2831406466515435*T^2*t)/2475880078570760549798248448 + (2831406466515435*T^2)/2475880078570760549798248448 - (2974109352427815*T*t^2)/39614081257132168796771975168 + (5957751106626237*T)/1237940039285380274899124224 - t^4/100000000000 + (4835703278458515*t^3)/2417851639229258349412352 - (7555786372591427*t^2)/37778931862957161709568 + t/100000 - 21/100000
 

% degrees = [1, 2, 3, 4];
% lambdas = [-1, 0, 1];
% 
% for degree = degrees
%     for lambda = lambdas
%         import lp4.verifyWithGivenLambda
%         [lpVer, solveResVer, resNorms] = verifyWithGivenLambda(vars, f, eps, g_theta, g_psy, g_zeta, lambda, degree);
%         
%         if solveResVer.hasSolution()
%             return
%         end
%     end
% end

warning('on')

echo off;
end
