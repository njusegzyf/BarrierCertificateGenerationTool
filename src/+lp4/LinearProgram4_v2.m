classdef LinearProgram4_v2
    %LINEARPROGRAM A linear program.
    
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
        pLambdaNormalizedSymbolicVars
        wSymbolicVars
        wExpression
        rouVar
        
        rouInc = 0.1
        
        pPartitions
        pLambdaPartitions
        
        linprogF
        
        theta
        psy
        zeta
    end % properties
    
    methods
        
        function this = LinearProgram4_v2(indvarsArg)
            % vars can only be a vector of a matrix of symbolic variables
            if ~isa(indvarsArg, 'sym')
                error('vars is not a vector of a matrix of symbolic variables');
            end
            
            this.indvars = reshape(indvarsArg, 1, size(indvarsArg, 1) * size(indvarsArg, 2));
            
            this.degree = 0;
            this.pLambdaDegree = 0;
            this.phy = 0;
            
            this.eps = [];
            
            this.f = [];
            this.decvars = [];
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
            
            expr = Constraint();
            expr.num = length(this.exprs) + 1;
            expr.name = 'theta';
            expr.type = 'ie';
            expr.A = [];
            expr.b = [];
            
            constraint1 = -this.phy;
            de = computeDegree(constraint1, this.indvars);
            
            c_alpha_beta = sym('c_alpha_beta', [1,10000]); % pre-defined varibales, only a few of them are the actual variables
            [constraintDecvars, expression] = constraintExpression(de, theta, c_alpha_beta);
            this = this.addDecisionVars(constraintDecvars);
            
            constraint1 = constraint1 + expression;
            
            constraint1 = expand(constraint1);
            expr.polyexpr = constraint1;
            
            this.exprs = [this.exprs expr];
        end
        
        function this = setPsyConstraint(this, psy)
            this.psy = psy;
            
            expr = Constraint();
            expr.num = length(this.exprs) + 1;
            expr.name = 'psy';
            expr.type = 'ie';
            expr.A = [];
            expr.b = [];
            
            phy_d = 0;
            for k = 1 : 1 : length(this.indvars)
                phy_d = phy_d + diff(this.phy, this.indvars(k)) * this.f(k);
            end
            
            constraint2 = -phy_d + this.wExpression + this.eps(1);
            de = computeDegree(constraint2, this.indvars);
            % make de even
            
            c_gama_delta = sym('c_gama_delta',[1,10000]);
            [constraintDecvars, expression] = constraintExpression(de, psy, c_gama_delta);
            this = this.addDecisionVars(constraintDecvars);
            
            constraint2 = constraint2 + expression;
            
            constraint2 = expand(constraint2);
            expr.polyexpr = constraint2;
            
            this.exprs = [this.exprs expr];
        end
        
        function this = setZetaConstraint(this, zeta)
            this.zeta = zeta;
            
            expr = Constraint();
            expr.num = length(this.exprs) + 1;
            expr.name = 'zeta';
            expr.type = 'ie';
            expr.A = [];
            expr.b = [];
            
            constraint3 = this.phy + this.eps(2); % different from lp2
            de = computeDegree(constraint3, this.indvars);
            
            c_u_v = sym('c_u_v',[1,10000]);
            [constraintDecvars, expression] = constraintExpression(de, zeta, c_u_v);
            this = this.addDecisionVars(constraintDecvars);
            
            constraint3 = constraint3 + expression;
            
            constraint3 = expand(constraint3);
            expr.polyexpr = constraint3;
            
            this.exprs = [this.exprs expr];
        end
        
        function this = generateEqsForConstraint1To3(this)
            rouIndex = this.getRouIndex();
            for k = 1 : 1 : 3
                [ A, b ] = eqgenerate( this.indvars, this.decvars, this.exprs(k).polyexpr);
                % convert to two ieq : -A*x <= -b, A*x <= b
                aLength = size(A, 1);
                
                A1 = -A;
                for i = 1 : aLength
                    A1(i, rouIndex) = -1;
                end
                b1 = -b + this.rouInc;
                
                A2 = A;
                for i = 1 : aLength
                    A2(i, rouIndex) = -1;
                end
                b2 = b + this.rouInc;
                
                this.exprs(k).A = [A1; A2];
                this.exprs(k).b = [b1; b2];
                
                disp(['constraint ',this.exprs(k).name,' is processed: ',datestr(now,'yyyy-mm-dd HH:MM:SS')]);
            end
        end
        
        function this = setWConstraint(this)
            expr = Constraint();
            expr.num = length(this.exprs) + 1;
            expr.name = 'w';
            expr.type = 'ie';
            expr.A = [];
            expr.b = [];
            
            for i = 1 : 1 : length(this.phyPolynomial.coefficientVars)
                p = this.phyPolynomial.coefficientVars(i);
                pPartition = this.pPartitions(i);
                
                expr = pPartition.createConstraintsAndAddToExpr(p, this.decvars, expr);
            end
            
            for i = 1 : 1 : length(this.pLambdaPolynomial.coefficientVars)
                p = this.pLambdaPolynomial.coefficientVars(i);
                pPartition = this.pLambdaPartitions(i);
                
                expr = pPartition.createConstraintsAndAddToExpr(p, this.decvars, expr);
            end
            
            for i1 = 1 : 1 : length(this.phyPolynomial.coefficientVars)
                for i2 = 1 : 1 : length(this.pLambdaPolynomial.coefficientVars)
                    p = this.phyPolynomial.coefficientVars(i1);
                    pPartition = this.pPartitions(i1);
                    pLambda = this.pLambdaPolynomial.coefficientVars(i2);
                    pLambdaPartition = this.pLambdaPartitions(i2);
                    w = this.wSymbolicVars(i1, i2);
                    
                    import lp4.createOverApproximationAsConstraints
                    constraints = createOverApproximationAsConstraints(p, pPartition, pLambda, pLambdaPartition, w);
                    
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
            
            this.exprs = [this.exprs expr];
        end
        
        function this = setDevVarsConstraint(this)
            decexpr = Constraint();
            decexpr.num = length(this.exprs) + 1;
            decexpr.name = 'decvarconstraints';
            decexpr.type = 'ie';
            decexpr.polyexpr = [];
            
            cStart = this.getWStart() + this.getWLength() - 1;
            cLength = length(this.decvars) - cStart;
            decexpr.A = zeros(cLength, length(this.decvars));
            rouIndex = this.getRouIndex();
            for k = 1 : 1 : cLength
                decexpr.A(k, cStart + k) = -1;
                % - rou
                decexpr.A(k, rouIndex) = -1;
            end
            bc = repmat(this.rouInc, cLength, 1);
            decexpr.b = bc;
            this.exprs = [this.exprs decexpr];
        end % function setDevVarsConstraint
        
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
        
        function [verification, res] = verify(this, solveRes)
            import lp4.LinearProgram4Verification2
            verification = LinearProgram4Verification2(this.indvars);
            
            verification.f = this.f;
            verification.eps = this.eps;
            
            % Set the degree of phy
            verification = verification.setDegreeAndInit(this.degree);
            
            verification.lambda = solveRes.getPLmabdaExpression();
            
            verification = verification.setThetaConstraint(this.theta);
            verification = verification.setPsyConstraint(this.psy);
            verification = verification.setZetaConstraint(this.zeta);
            verification = verification.generateEqsForConstraint1To3();
            
            verification = verification.setDevVarsConstraint();
            
            % solve the lp problem
            [verification, res] = verification.solve();
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
