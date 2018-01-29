classdef LinearProgram4_v3
    %LINEARPROGRAM A linear program.
    % In this one, we only add rou to all ies.
    % And for resue a lp for different ranges, we let the constraints of
    % (p, c, w) to be the last expr (index 5).
    % We also allow theta, psy, zeta to be empty.
    
    properties
        indvars % 问题1中的独立变量，例如x = [x1, x2, ..., xn]，类型为符号化变量构成的行向量
        degree % 问题1中的带求多项式函数φ的次数，类型为正整数
        phy % 问题1中的带求多项式函数φ，当lp.degree初始化后，自动赋值为系数为决策变量的多项式
        % r % 问题1中的λ, which is replaced by pLambdaPolynomial
        eps % 问题1中的?1和?1构成的向量
        f % 问题1中的f
        decvars % 问题1中的决策变量，包括P, Cα,β, Cγ,δ, Cu,v
        exprs
        
        % new in lp4
        phyPolynomial
        
        pLambdaDegree
        pLambdaPolynomial % 修改后的问题1中的λ
        
        wSymbolicVars
        wExpression
        rouVar
        
        pPartitions
        pLambdaPartitions
        
        linprogF
        
        theta
        psy
        zeta
        
        rouInc  =  0
    end % properties
    
    methods
        
        function this = LinearProgram4_v3(indvarsArg)
            % vars can only be a vector of a matrix of symbolic variables
            if ~isa(indvarsArg, 'sym')
                error('vars is not a vector of a matrix of symbolic variables');
            end
            
            this.indvars = reshape(indvarsArg, 1, size(indvarsArg, 1) * size(indvarsArg, 2));
            
            this.degree = 0;
            this.pLambdaDegree = 0;
            
            this.eps = [];
            
            this.f = [];
        end
        
        function this = initLpWithoutRanges(this, f, eps, theta, psy, zeta, degree, pLambdaDegree)
            this.f = f;
            this.eps = eps;
            
            this = this.setDegreeAndInit(degree, pLambdaDegree);
            
            this = this.setThetaConstraint(theta);
            this = this.setPsyConstraint(psy);
            this = this.setZetaConstraint(zeta);
            this = this.generateEqsForConstraint1To3();
            
            this = this.setDevVarsConstraint();
            
            this = this.setLinprogF();
        end
        
        function this = initLp(this, f, eps, theta, psy, zeta, degree, pLambdaDegree, phyRange, pLambdaRange)
            this = this.initLpWithoutRanges(f, eps, theta, psy, zeta, degree, pLambdaDegree);
            
            import lp4util.Partition
            this.pPartitions = repmat(phyRange, 1024, 1);
            this.pLambdaPartitions = repmat(pLambdaRange, 1024, 1);
            this = this.setWConstraint();
        end
        
        function lp = set.f(lp, f)
            lp.f = f;
        end
        
        function lp = set.eps(lp, eps)
            lp.eps = eps;
        end
        
        function this = addDecisionVars(this, decvars)
            % addDecisionVars add decision variables.
            %
            % lp is a linear program.
            % decvars are decsion variables and should be symbolic.
            
            if isa(decvars, 'sym')
                % 将新的决策变量变为行向量的形式添加到lp的lp.decvars属性中
                % decision vars can only be a vector of a matrix of symbolic variables
                % reshape `decvars` to of dim `[1, size(decvars, 1) * size(decvars, 2)]`
                this.decvars = [this.decvars reshape(decvars, 1, size(decvars, 1) * size(decvars, 2))];
            end
        end
        
        function this = setDegreeAndInit(this, degree, pLambdaDegree)
            this.degree = degree;
            this.pLambdaDegree = pLambdaDegree;
            p = sym('p', [1, monomialNumber(length(this.indvars), degree)]);
            this = this.addDecisionVars(p);
            this.phy = p * monomials(this.indvars, 0 : degree);
            
            % new in lp4
            import lp4util.SymbolicPolynomial
            
            this.rouVar = sym('rou');
            this = this.addDecisionVars(this.rouVar);
            
            this.phyPolynomial = SymbolicPolynomial(this.indvars, degree, p, this.phy);
            
            % set PLambda
            pLambda = sym('pLambda', [1, monomialNumber(length(this.indvars), pLambdaDegree)]);
            this = this.addDecisionVars(pLambda);
            pLambdaExpr = pLambda * monomials(this.indvars, 0 : pLambdaDegree);
            this.pLambdaPolynomial = SymbolicPolynomial(this.indvars, pLambdaDegree, pLambda, pLambdaExpr);
            
            % introduce W = P * PLambda
            import lp4.LinearProgram4
            [this.wSymbolicVars, this.wExpression] = LinearProgram4.createWExpression(this.phyPolynomial, this.pLambdaPolynomial);
            this = this.addDecisionVars(this.wSymbolicVars);
        end
        
        function this = setThetaConstraint(this, theta)
            this.theta = theta;
            
            % for an empty constrain
            if isempty(theta)
                expr = Constraint.createEmptyConstraint();
                expr.num = 1;
                this.exprs = [this.exprs, expr];
                return;
            end
            
            expr = Constraint();
            expr.num = 1;
            expr.name = 'theta';
            expr.type = 'eq';
            expr.A = [];
            expr.b = [];
            
            constraint1 = -this.phy;
            import lp4.Lp4Config
            de = computeDegree(constraint1, this.indvars) + Lp4Config.C_DEGREE_INC;
            
            c_alpha_beta = sym('c_alpha_beta', [1,10000]); % pre-defined varibales, only a few of them are the actual variables
            [constraintDecvars, expression] = constraintExpression(de, theta, c_alpha_beta);
            this = this.addDecisionVars(constraintDecvars);
            
            constraint1 = constraint1 + expression;
            
            constraint1 = expand(constraint1);
            expr.polyexpr = constraint1;
            
            % Note: This is the first expression, and use `this.exprs(expr.num) = expr;` here
            % will cause an error in Matlab 2014b but is Ok in Matlab 2017b.
            this.exprs = [this.exprs, expr];
        end
        
        function this = setPsyConstraint(this, psy)
            this.psy = psy;
            
            % for an empty constrain
            if isempty(psy)
                expr = Constraint.createEmptyConstraint();
                expr.num = 2;
                this.exprs(expr.num) = expr;
                return;
            end
            
            expr = Constraint();
            expr.num = 2;
            expr.name = 'psy';
            expr.type = 'eq';
            expr.A = [];
            expr.b = [];
            
            phy_d = 0;
            for k = 1 : 1 : length(this.indvars)
                phy_d = phy_d + diff(this.phy, this.indvars(k)) * this.f(k);
            end
            
            constraint2 = -phy_d + this.wExpression + this.eps(1);
            import lp4.Lp4Config
            de = computeDegree(constraint2, this.indvars) + Lp4Config.C_DEGREE_INC;
            
            c_gama_delta = sym('c_gama_delta',[1, 1000 * 1000]);
            [constraintDecvars, expression] = constraintExpression(de, psy, c_gama_delta);
            this = this.addDecisionVars(constraintDecvars);
            
            constraint2 = constraint2 + expression;
            
            constraint2 = expand(constraint2);
            expr.polyexpr = constraint2;
            
            this.exprs(expr.num) = expr;
        end
        
        function this = setZetaConstraint(this, zeta)
            this.zeta = zeta;
            
            % for an empty constrain
            if isempty(zeta)
                expr = Constraint.createEmptyConstraint();
                expr.num = 3;
                this.exprs(expr.num) = expr;
                return;
            end
            
            expr = Constraint();
            expr.num = 3;
            expr.name = 'zeta';
            expr.type = 'eq';
            expr.A = [];
            expr.b = [];
            
            constraint3 = this.phy + this.eps(2); % different from lp2
            import lp4.Lp4Config
            de = computeDegree(constraint3, this.indvars) + Lp4Config.C_DEGREE_INC;
            
            c_u_v = sym('c_u_v',[1,10000]);
            [constraintDecvars, expression] = constraintExpression(de, zeta, c_u_v);
            this = this.addDecisionVars(constraintDecvars);
            
            constraint3 = constraint3 + expression;
            
            constraint3 = expand(constraint3);
            expr.polyexpr = constraint3;
            
            this.exprs(expr.num) = expr;
        end
        
        function this = generateEqsForConstraint1To3(this)
            for k = 1 : 1 : 3
                if this.exprs(k).isEmptyConstraint()
                    continue;
                end
                
                [ this.exprs(k).A, this.exprs(k).b ] = eqgenerate( this.indvars, this.decvars, this.exprs(k).polyexpr);
                disp(['constraint ',this.exprs(k).name,' is processed: ',datestr(now,'yyyy-mm-dd HH:MM:SS')]);
            end
        end
        
        function this = setDevVarsConstraint(this)
            expr = Constraint();
            expr.num = 4;
            expr.name = 'decvarconstraints';
            expr.type = 'ie';
            expr.polyexpr = [];
            
            cStart = this.getWStart() + this.getWLength() - 1;
            cLength = length(this.decvars) - cStart;
            expr.A = zeros(cLength, length(this.decvars));
            rouIndex = this.getRouIndex();
            for k = 1 : 1 : cLength
                expr.A(k, cStart + k) = -1;
                % - rou
                expr.A(k, rouIndex) = -1;
            end
            bc = repmat(this.rouInc, cLength, 1);
            expr.b = bc;
            
            this.exprs(expr.num) = expr;
        end % function setDevVarsConstraint
        
        function this = setWConstraint(this)
            expr = Constraint();
            expr.num = 5;
            expr.name = 'w';
            expr.type = 'ie';
            expr.A = [];
            expr.b = [];
            
            rouIndex = this.getRouIndex();
            for i = 1 : 1 : length(this.phyPolynomial.coefficientVars)
                p = this.phyPolynomial.coefficientVars(i);
                pPartition = this.pPartitions(i);
                
                expr = pPartition.createConstraintsWithRouAndAddToExpr(p, this.decvars, rouIndex, expr);
            end
            
            for i = 1 : 1 : length(this.pLambdaPolynomial.coefficientVars)
                p = this.pLambdaPolynomial.coefficientVars(i);
                pPartition = this.pLambdaPartitions(i);
                
                expr = pPartition.createConstraintsWithRouAndAddToExpr(p, this.decvars, rouIndex, expr);
            end
            
            for i1 = 1 : 1 : length(this.phyPolynomial.coefficientVars)
                for i2 = 1 : 1 : length(this.pLambdaPolynomial.coefficientVars)
                    p = this.phyPolynomial.coefficientVars(i1);
                    pPartition = this.pPartitions(i1);
                    pLambda = this.pLambdaPolynomial.coefficientVars(i2);
                    pLambdaPartition = this.pLambdaPartitions(i2);
                    w = this.wSymbolicVars(i1, i2);
                    
                    import lp4.createOverApproximationAsConstraintsWithRou
                    constraints = createOverApproximationAsConstraintsWithRou(...
                        p, pPartition, pLambda, pLambdaPartition, w, this.rouVar);
                    
                    % save expressions for debug
                    import lp4.Lp4Config
                    if Lp4Config.isDebug()
                        expr.polyexpr = [expr.polyexpr constraints];
                    end
                    
                    % constraints = arrayfun(@(x) x == 0, constraints);
                    [ A, b ] = equationsToMatrix(constraints, this.decvars);
                    A = double(A);
                    b = double(b);
                    expr.A = [expr.A; A];
                    expr.b = [expr.b; b];
                end
            end
            
            this.exprs(expr.num) = expr;
        end
        
        function this = setLinprogF(this)
            % init f with all zeros
            this.linprogF = zeros(1, length(this.decvars));
            
            % chagend in lp4, in which f = rou
            % rouVarIndex = find(this.decvars == this.rouVar);
            rouVarIndex = this.getRouIndex();
            this.linprogF(rouVarIndex) = 1;
        end
        
        function [this, solveRes] = solve(this)
            Aeq = [];
            beq = [];
            Aie = [];
            bie = [];
            
            for k = 1 : 1 : length(this.exprs)
                if this.exprs(k).isEmptyConstraint()
                    continue;
                end
                
                if strcmp(this.exprs(k).type, 'eq')
                    Aeq = [Aeq; this.exprs(k).A];
                    beq = [beq; this.exprs(k).b];
                else
                    Aie = [Aie; this.exprs(k).A];
                    bie = [bie; this.exprs(k).b];
                end
            end
            
            tic;
            [x, fval, flag, ~] = linprog(this.linprogF, Aie, bie, Aeq, beq);
            time = toc;
            
            import lp4.LinearProgram4SolveResult
            solveRes = LinearProgram4SolveResult(this, x, fval, flag, time);
            
            import lp4.Lp4Config
            if Lp4Config.isDebug()
                solveRes.printSolution();
            end
            
        end % function solve
        
        function [this, solveRes] = solve1And3(this)
            Aeq = [];
            beq = [];
            Aie = [];
            bie = [];
            
            for k = [1, 3, 5]
                if strcmp(this.exprs(k).name, 'empty')
                    continue;
                end
                
                if strcmp(this.exprs(k).type, 'eq')
                    Aeq = [Aeq; this.exprs(k).A];
                    beq = [beq; this.exprs(k).b];
                else
                    Aie = [Aie; this.exprs(k).A];
                    bie = [bie; this.exprs(k).b];
                end
            end
            
            % set rou, pLambda, w to 0
            unusedVarStart = this.getRouIndex() - 1;
            unusedVarLength = 1 + this.getPLambdaCoefficientLength() + this.getWLength();
            unusedVarAeq = zeros(unusedVarLength, length(this.decvars));
            for k = 1 : 1 : unusedVarLength
                unusedVarAeq(k, unusedVarStart + k) = 1;
            end
            unusedVarbeq = zeros(unusedVarLength, 1);
            
            Aeq = [Aeq; unusedVarAeq];
            beq = [beq; unusedVarbeq];
            
            % solve
            tic;
            [x, fval, flag, ~] = linprog(this.linprogF, Aie, bie, Aeq, beq);
            time = toc;
            
            import lp4.LinearProgram4SolveResult
            solveRes = LinearProgram4SolveResult(this, x, fval, flag, time);
            
            import lp4.Lp4Config
            if Lp4Config.isDebug()
                solveRes.printSolution();
            end
            
        end % function solve
        
        % `verify` is forward to `verifyWithLambda`
        function [lpVer, solveResVer, resNorms] = verify(lp, solveRes, phyRangeInVerify)
            [lpVer, solveResVer, resNorms] = verifyWithLambda(lp, solveRes, phyRangeInVerify);
        end
        
        function [lpVer, solveResVer, resNorms] = verifyWithLambda(lp, solveRes, phyRangeInVerify)
            if ~(solveRes.hasSolution())
                lpVer = 0;
                solveResVer = 0;
                resNorms = [];
                disp('Not find feasible solution to verify.');
                return;
            end
            
            % verify
            
            import lp4.LinearProgram4Verification2
            lpVer = LinearProgram4Verification2(lp.indvars);
            
            lpVer.f = lp.f;
            lpVer.eps = lp.eps;
            
            % set the degree of phy
            import lp4.Lp4Config
            lpVer = lpVer.setDegreeAndInit(lp.degree + Lp4Config.VERIFICATION_PHY_DEGREE_INC);
            
            % set lambda expression
            lpVer.lambda = solveRes.getPLmabdaExpression();
            
            lpVer = lpVer.setThetaConstraint(lp.theta);
            lpVer = lpVer.setPsyConstraint(lp.psy);
            lpVer = lpVer.setZetaConstraint(lp.zeta);
            lpVer = lpVer.generateEqsForConstraint1To3();
            
            if isa(phyRangeInVerify, 'lp4util.Partition')
                lpVer.pPartitions = repmat(phyRangeInVerify, 1024, 1);
                lpVer = lpVer.setPhyConstraint();
            end
            
            lpVer = lpVer.setDevVarsConstraint();
            
            % solve the lp problem
            [lpVer, solveResVer, resNorms] = lpVer.solve();
            
            if ~(solveResVer.hasSolution())
                disp('Verify feasible solution failed.');
            else
                disp('Verify feasible solution succeed, norms ;');
                disp(resNorms);
            end
            
            import lp4.Lp4Config
            Lp4Config.displayDelimiterLine();
        end
        
        function [lpVer, solveResVer, resNorms] = verifyWithPhy(lp, solveRes)
            if ~(solveRes.hasSolution())
                lpVer = 0;
                solveResVer = 0;
                resNorms = [];
                disp('Not find feasible solution to verify.');
                return;
            end
            
            % verify
            
            import lp4.LinearProgram4Verification3
            lpVer = LinearProgram4Verification3(lp.indvars);
            
            lpVer.f = lp.f;
            lpVer.eps = lp.eps;
            
            % set the degree of lambda
            import lp4.Lp4Config
            lpVer = lpVer.setDegreeAndInit(lp.pLambdaDegree + Lp4Config.VERIFICATION_LAMBDA_DEGREE_INC);
            
            % set phy expression
            lpVer.phy = solveRes.getPhyExpression();
            
            lpVer = lpVer.setThetaConstraint(lp.theta);
            lpVer = lpVer.setPsyConstraint(lp.psy);
            lpVer = lpVer.setZetaConstraint(lp.zeta);
            lpVer = lpVer.generateEqsForConstraint1To3();
            
            lpVer = lpVer.setDevVarsConstraint();
            
            [lpVer, solveResVer, resNorms] = lpVer.solve();
            
            if ~(solveResVer.hasSolution())
                disp('Verify feasible solution failed.');
            else
                disp('Verify feasible solution succeed, norms ;');
                disp(resNorms);
            end
            
            import lp4.Lp4Config
            Lp4Config.displayDelimiterLine();
        end
        
        function res = getPhyCoefficientStart(this)
            res = 1;
        end
        
        function res = getPhyCoefficientLength(this)
            res = length(this.phyPolynomial.coefficientVars);
        end
        
        function res = getRouIndex(this)
            res = this.getPhyCoefficientStart() + this.getPhyCoefficientLength();
        end
        
        function res = getPLambdaCoefficientStart(this)
            res = this.getRouIndex() + 1;
        end
        
        function res = getPLambdaCoefficientLength(this)
            res = length(this.pLambdaPolynomial.coefficientVars);
        end
        
        function res = getWStart(this)
            res = this.getPLambdaCoefficientStart() + this.getPLambdaCoefficientLength();
        end
        
        function res = getWLength(this)
            res = size(this.wSymbolicVars, 1) * size(this.wSymbolicVars, 2);
        end
        
    end % methods
    
    methods (Static)
        
        function lp = createLp(vars, f, eps, theta, psy, zeta, degree, pLambdaDegree, phyRange, pLambdaRange)
            import lp4.LinearProgram4_v3
            lp = LinearProgram4_v3(vars);
            lp = lp.initLp(f, eps, theta, psy, zeta, degree, pLambdaDegree, phyRange, pLambdaRange);
        end
        
        function lp = createLpWithoutRanges(vars, f, eps, theta, psy, zeta, degree, pLambdaDegree)
            import lp4.LinearProgram4_v3
            lp = LinearProgram4_v3(vars);
            lp = lp.initLpWithoutRanges(f, eps, theta, psy, zeta, degree, pLambdaDegree);
        end
        
        function [wSymbolicVars, wExpression] = createWExpression(pPolynomial, pLambdaPolynomial)
            wSymbolicVars = sym('w', [length(pPolynomial.coefficientVars), length(pLambdaPolynomial.coefficientVars)]);
            wExpression = pPolynomial.expression * pLambdaPolynomial.expression;
            
            wExpression = expand(wExpression); % feval(symengine, 'expand', wExpression);
            
            % wExpression = collect(wExpression);
            % wExpression = feval(symengine, 'collect', wExpression,...
            % converttochar([pPolynomial.symbolicVars, pLambdaPolynomial.symbolicVars]));
            
            for i = 1 : 1 : length(pPolynomial.coefficientVars)
                for j = 1 : 1 : length(pLambdaPolynomial.coefficientVars)
                    % replace `p_i * pLambda_j` with `w_i_j`
                    iSymbol = pPolynomial.coefficientVars(i);
                    jSymbol = pLambdaPolynomial.coefficientVars(j);
                    wSymbol = wSymbolicVars(i, j);
                    wExpression = subs(wExpression, iSymbol * jSymbol, wSymbol);
                end
            end
        end % function createWExpression
        
    end % methods (Static)
    
end % classdef
