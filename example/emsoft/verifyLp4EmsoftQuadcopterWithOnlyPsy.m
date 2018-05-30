function [lpVer, solveResVer, resNorms] = verifyLp4EmsoftQuadcopterWithOnlyPsy()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, ~, psy, ~] = getLp4EmsoftQuadcopterProblem();

x1 = vars(1);
x2 = vars(2);
x3 = vars(3);
x4 = vars(4);
x5 = vars(5);
x6 = vars(6);
x7 = vars(7);
x8 = vars(8);
x9 = vars(9);
x10 = vars(10);
x11 = vars(11);
x12 = vars(12);

initPhy = (6125460164004971*x1^2)/1125899906842624 + (3658882615012259*x1*x2)/9903520314283042199192993792 - (2385093345331*x1*x3)/9671406556917033397649408 + (1173700150021249*x1*x4)/288230376151711744 + (16903800292119161*x1*x5)/5070602400912917605986812821504 - (10573170396772063*x1*x6)/20282409603651670423947251286016 + (3125104702180361*x1*x7)/2475880078570760549798248448 - (4512769157139313*x1*x8)/2251799813685248 + (399969536153861*x1*x9)/9903520314283042199192993792 + (6436898763663335*x1*x10)/9903520314283042199192993792 - (13100782248304567*x1*x11)/4503599627370496 + (17264417196630933*x1*x12)/2535301200456458802993406410752 + (6125460164005083*x2^2)/1125899906842624 - (142806315609*x2*x3)/75557863725914323419136 - (6513014589398907*x2*x4)/1267650600228229401496703205376 + (4694800600101863*x2*x5)/1152921504606846976 - (15885801231950227*x2*x6)/5070602400912917605986812821504 + (564096144642305*x2*x7)/281474976710656 + (5660816936221259*x2*x8)/4951760157141521099596496896 + (567458303955*x2*x9)/9671406556917033397649408 + (13100782248303719*x2*x10)/4503599627370496 + (852105435113047*x2*x11)/4951760157141521099596496896 + (8369090473342829*x2*x12)/158456325028528675187087900672 + (1126631741165689*x3^2)/1125899906842624 - (2468276917894969*x3*x4)/40564819207303340847894502572032 - (5650665813690413*x3*x5)/10141204801825835211973625643008 + (2997595911977533*x3*x6)/2305843009213693952 - (2797598577*x3*x7)/9444732965739290427392 + (785228817*x3*x8)/9444732965739290427392 - (2414730501763*x3*x9)/618970019642690137449562112 - (3727197428913271*x3*x10)/19807040628566084398385987584 + (5546841058673875*x3*x11)/79228162514264337593543950336 - (6214872250259387*x3*x12)/2535301200456458802993406410752 + (6917531029210581*x4^2)/18446744073709551616 + (3293716064124701*x4*x5)/81129638414606681695789005144064 - (2374304892175305*x4*x6)/2596148429267413814265248164610048 - (11012115521398003*x4*x7)/40564819207303340847894502572032 - (13835058055039847*x4*x8)/18446744073709551616 - (5402535322100043*x4*x9)/5070602400912917605986812821504 - (16598254578052589*x4*x10)/40564819207303340847894502572032 - (5020480818534999*x4*x11)/4611686018427387904 - (11039278094641599*x4*x12)/162259276829213363391578010288128 + (1729382757303101*x5^2)/4611686018427387904 - (70544152146249*x5*x6)/40564819207303340847894502572032 + (864691128457485*x5*x7)/1152921504606846976 - (6260299247955117*x5*x8)/5070602400912917605986812821504 - (9597626794418619*x5*x9)/1298074214633706907132624082305024 + (10040961637139021*x5*x10)/9223372036854775808 - (9640620747106241*x5*x11)/10141204801825835211973625643008 + (1452715009595571*x5*x12)/40564819207303340847894502572032 + (5995196876753173*x6^2)/9223372036854775808 - (1137698469066805*x6*x7)/5070602400912917605986812821504 - (11336567300837669*x6*x8)/162259276829213363391578010288128 + (9868741422229981*x6*x9)/162259276829213363391578010288128 - (142228739727579*x6*x10)/316912650057057350374175801344 + (5308959724894443*x6*x11)/40564819207303340847894502572032 + (7722837321408001*x6*x12)/2596148429267413814265248164610048 + (3269385089018433*x7^2)/2251799813685248 + (3451921197863257*x7*x8)/4951760157141521099596496896 + (10684630992076553*x7*x9)/633825300114114700748351602688 + (623753909514429*x7*x10)/562949953421312 + (1637574801149411*x7*x11)/9903520314283042199192993792 + (10503027216931833*x7*x12)/633825300114114700748351602688 + (6538770178036453*x8^2)/4503599627370496 + (12474783014528933*x8*x9)/79228162514264337593543950336 + (8002821255025837*x8*x10)/19807040628566084398385987584 + (2495015638057487*x8*x11)/2251799813685248 + (14981779813089907*x8*x12)/2535301200456458802993406410752 + (4649701592162465*x9^2)/4503599627370496 + (9768488797463847*x9*x10)/633825300114114700748351602688 + (6784730002109699*x9*x11)/79228162514264337593543950336 + (9367487224931077*x9*x12)/144115188075855872 + (6324983308235203*x10^2)/9007199254740992 + (9520490369815455*x10*x11)/158456325028528675187087900672 + (17047789820904799*x10*x12)/1267650600228229401496703205376 + (6324983308235001*x11^2)/9007199254740992 + (9742891414709511*x11*x12)/2535301200456458802993406410752 + (1173122889175201*x12^2)/36028797018963968 - 1;

pLambdaDegree = 0;



% run and verify
lpVer = lp4.LinearProgram4Verification3.createWithRou(vars, f, eps, [], psy, [], pLambdaDegree, initPhy);

% set phy by hand
constraint2 = -(x1^2 + x2^2 + x3^2 + x4^2 + x5^2 + x6^2 + x7^2 + x8^2 + x9^2 + x10^2 + x11^2 + x12^2) ; % + initPhy + lpVer.eps(1);
leftDegree = computeDegree(constraint2, lpVer.indvars);
de = lp4.Lp4Config.getVerificationCDegree(leftDegree);

c_gama_delta = sym('c_gama_delta', [1, lp4.Lp4Config.DEFAULT_DEC_VAR_SIZE]);
[~, expression] = constraintExpression(de, psy, c_gama_delta);
constraint2 = constraint2 + expression;

lpVer.exprs(2).polyexpr = constraint2;

[lpVer, solveResVer, resNorms] = lpVer.solve();

lp4.Lp4Config.printVerifyWithOnlyPsyResult(solveResVer, resNorms);

warning('on')

echo off;
end
