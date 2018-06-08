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
        decvars % 问题1中的决策变量，包括 p, lambda, w, rou, c
        exprs
        
        % new in lp4
        phyPolynomial
        
        pLambdaDegree
        pLambdaPolynomial % 修改后的问题1中的λ
        
        wSymbolicVars
        wExpression
        
        rouVarIndex
        rouVar
        
        pPartitions
        pLambdaPartitions
        
        linprogF
        
        theta
        psy
        zeta
        
        rouInc  =  0
        
        cStart = -1
        lambdaEnd
        
        fDecvars
        nonlconConstraints % nonlconConstraints(i) <= 0
        nonlconConstraintsFunctions
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
            
            % set Phy
            p = sym('p', [1, monomialNumber(length(this.indvars), degree)]);
            this = this.addDecisionVars(p);
            this.phy = p * monomials(this.indvars, 0 : degree);
            
            % new in lp4
            import lp4util.SymbolicPolynomial
            this.phyPolynomial = SymbolicPolynomial(this.indvars, degree, p, this.phy);
            
            % set PLambda
            pLambda = sym('pLambda', [1, monomialNumber(length(this.indvars), pLambdaDegree)]);
            this = this.addDecisionVars(pLambda);
            pLambdaExpr = pLambda * monomials(this.indvars, 0 : pLambdaDegree);
            this.pLambdaPolynomial = SymbolicPolynomial(this.indvars, pLambdaDegree, pLambda, pLambdaExpr);
            this.lambdaEnd = length(this.decvars);
            
            % introduce W = P * PLambda
            [this.wSymbolicVars, this.wExpression] = lp4.LinearProgram4_v3.createWExpression(this.phyPolynomial, this.pLambdaPolynomial);
            this = this.addDecisionVars(this.wSymbolicVars);
            
            this.rouVar = sym('rou');
            this = this.addDecisionVars(this.rouVar);
            this.rouVarIndex = length(this.decvars);
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
            de = computeDegree(constraint1, this.indvars) + lp4.Lp4Config.C_DEGREE_INC;
            
            c_alpha_beta = sym('c_alpha_beta', [1, lp4.Lp4Config.DEFAULT_DEC_VAR_SIZE]); % pre-defined varibales, only a few of them are the actual variables
            [constraintDecvars, expression] = constraintExpression(de, theta, c_alpha_beta);
            this.cStart = length(this.decvars) + 1;
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
            de = computeDegree(constraint2, this.indvars) + lp4.Lp4Config.C_DEGREE_INC;
            
            c_gama_delta = sym('c_gama_delta', [1, lp4.Lp4Config.DEFAULT_DEC_VAR_SIZE]);
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
            de = computeDegree(constraint3, this.indvars) + lp4.Lp4Config.C_DEGREE_INC;
            
            c_u_v = sym('c_u_v', [1, lp4.Lp4Config.DEFAULT_DEC_VAR_SIZE]);
            [constraintDecvars, expression] = constraintExpression(de, zeta, c_u_v);
            this = this.addDecisionVars(constraintDecvars);
            
            constraint3 = constraint3 + expression;
            
            constraint3 = expand(constraint3);
            expr.polyexpr = constraint3;
            
            this.exprs(expr.num) = expr;
        end
        
        function this = generateEqsForConstraint1To3(this)
            this = lp4.Lp4AndHlpVerificationBase.generateConstraintEqsParallelly(this);
        end
        
        function this = setDevVarsConstraint(this)
            expr = Constraint();
            expr.num = 4;
            expr.name = 'decvarconstraints';
            expr.type = 'ie';
            expr.polyexpr = [];
            
            cStartMinus1 = this.cStart - 1;
            cLength = length(this.decvars) - cStartMinus1;
            expr.A = zeros(cLength, length(this.decvars));
            rouIndex = this.getRouIndex();
            for k = 1 : 1 : cLength
                expr.A(k, cStartMinus1 + k) = -1;
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
                    if lp4.Lp4Config.isDebug()
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
            this.linprogF(this.getRouIndex()) = 1;
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
            
            if lp4.Lp4Config.isDebug()
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
            
            solveRes = lp4.LinearProgram4SolveResult(this, x, fval, flag, time);
            
            if lp4.Lp4Config.isDebug()
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
            
            lpVer = lp4.LinearProgram4Verification2.create(lp.indvars,...
                lp.f, lp.eps, lp.theta, lp.psy, lp.zeta, lp.degree, solveRes.getPLmabdaExpression(), phyRangeInVerify);
            
            % solve the lp problem
            [lpVer, solveResVer, resNorms] = lpVer.solve();
            
            if ~(solveResVer.hasSolution())
                disp('Verify feasible solution failed.');
            else
                disp('Verify feasible solution succeed, norms ;');
                disp(resNorms);
            end
            
            lp4.Lp4Config.displayDelimiterLine();
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
            
            lpVer = lp4.LinearProgram4Verification3.create(lp.indvars,...
                lp.f, lp.eps, lp.theta, lp.psy, lp.zeta, lp.pLambdaDegree, solveRes.getPhyExpression());
            
            [lpVer, solveResVer, resNorms] = lpVer.solve();
            
            if ~(solveResVer.hasSolution())
                disp('Verify feasible solution failed.');
            else
                disp('Verify feasible solution succeed, norms ;');
                disp(resNorms);
            end
            
            lp4.Lp4Config.displayDelimiterLine();
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
        
        function res = getWEnd(this)
            res = this.rouVarIndex - 1;
        end
        
        function res = getAllCDecvars(this)
            res = this.decvars(this.cStart : end);
        end
        
        function res = solveDirectlyWithFmincon(this)
            phyAndLambdaVars = this.decvars(1 : this.lambdaEnd);
            phyAndLambdaVarsLen = length(phyAndLambdaVars);
            
            cVars = this.getAllCDecvars();
            cVarsLen = length(cVars);
            
            this.fDecvars = [phyAndLambdaVars, cVars];
            fDecvarsLen = phyAndLambdaVarsLen + cVarsLen;
            
            % optimtool
            % http://ww2.mathworks.cn/help/optim/ug/fmincon.html
            % http://www.mathworks.com/help/optim/ug/optimoptions.html
            
            % all available algorithms: 'interior-point', 'active-set', 'sqp', 'sqp-legacy', 'trust-region-reflective'
            problem.solver = 'fmincon';
            options = optimoptions('fmincon', 'Display', 'iter',...
                                   'Algorithm', 'interior-point',...
                                   'MaxFunctionEvaluations', 50000, 'MaxIterations', 200);
            problem.options = options;
            
            % set the objective function to maximize the sum of all Cs (to be the negative of the sum of all C)
            problem.objective = @(x) -sum(x(phyAndLambdaVarsLen + 1 : end));
            
            % init point
            problem.x0 = zeros([1, fDecvarsLen]); % ones([1, fDecvarsLen]);
            
            % no lb and up
            %             indvarsLbs = arrayfun(@(p) p.boundLow, indvarsBounds);
            %             lb = [indvarsLbs, zeros([1, cVarsLen])];
            %             indvarsUbs = arrayfun(@(p) p.boundHigh, indvarsBounds);
            %             ub = [indvarsLbs, repmat(10000000, [1, cVarsLen])];
            %             problem.lb = lb;
            %             problem.ub = ub;
            
            % use ineq for c >= 0
            aineq = zeros([cVarsLen, fDecvarsLen]);
            for k = 1 : 1 : cVarsLen
                aineq(k, phyAndLambdaVarsLen + k) = -1;
            end
            bineq = zeros([cVarsLen, 1]);
            problem.Aineq = aineq;
            problem.bineq = bineq;
            
            this = this.setNonlconConstraints();
            this.displayNonlcon();
            
            problem.nonlcon = @this.getNonlcon;
            
            res = fmincon(problem);
            % res = []; % if we just want to display the constraints
        end
        
        function this = setNonlconConstraints(this)
            % collect eqs for non linear constraints
            Aeq = [];
            beq = [];
            
            for k = 1 : length(this.exprs)
                if this.exprs(k).isEmptyConstraint()
                    continue;
                end
                
                if strcmp(this.exprs(k).type, 'eq')
                    Aeq = [Aeq; this.exprs(k).A];
                    beq = [beq; this.exprs(k).b];
                end
            end
            
            decvarsLen = length(this.decvars);
            eqLength = size(Aeq, 1);
            
            % pre allocate memory
            fakeConstraint = this.decvars(1);
            this.nonlconConstraints = [fakeConstraint];
            this.nonlconConstraints(eqLength) = fakeConstraint;
            
            this.nonlconConstraintsFunctions = cell(eqLength);
            %             fakeFunction = matlabFunction(fakeConstraint);
            %             % Note: we must use cell arrays for function handlers
            %             % this.nonlconConstraintsFunctions = {fakeFunction};
            %             this.nonlconConstraintsFunctions{eqLength} = fakeFunction;
            
            for k = 1 : eqLength
                expr = Aeq(k, :) * reshape(this.decvars, [decvarsLen, 1]) - beq(k);
                
                for i1 = 1 : length(this.phyPolynomial.coefficientVars)
                    for i2 = 1 : length(this.pLambdaPolynomial.coefficientVars)
                        p = this.phyPolynomial.coefficientVars(i1);
                        pLambda = this.pLambdaPolynomial.coefficientVars(i2);
                        w = this.wSymbolicVars(i1, i2);
                        
                        expr = subs(expr, w, p*pLambda);
                    end
                end
                
                this.nonlconConstraints(k) = expr;
                this.nonlconConstraintsFunctions{k} = matlabFunction(expr, 'vars', {this.fDecvars});
            end
        end
        
        function [c, ceq] = getNonlcon(this, x)
            
            nonlconConstraintsLen = length(this.nonlconConstraints);
            
            % pre alocate memory
            ceq(nonlconConstraintsLen) = 0;
            
            for k = 1 : nonlconConstraintsLen
                
                % since we have converted symbol expressions to function handlers, we can directly call the functions.
                ceq(k) = this.nonlconConstraintsFunctions{k}(x);
                
                %                 expr = this.nonlconConstraints(k);
                %
                %                 % eval the constraint expr with its actual value
                %                 expr = subs(expr, this.fDecvars, x);
                %                 %                 fDecvarsLen = length(this.fDecvars);
                %                 %                 for i = 1 : fDecvarsLen
                %                 %                      expr = subs(expr, this.fDecvars(i), x(i));
                %                 %                 end
                %
                %                 ceq(k) = double(expr);
            end
            
            c = [];
        end
        
        function displayNonlcon(this)
            lp4.Lp4Config.displayDelimiterLine();
            disp(['Phy degree: ', num2str(this.degree), ', Lambda degree:', num2str(this.pLambdaDegree)]);
            
            lp4.Lp4Config.displayDelimiterLine();
            disp('Dec vars:');
            disp(this.fDecvars);
            
            lp4.Lp4Config.displayDelimiterLine();
            disp('Nonlinear constraints:');
            
            disp('[');
            for expr = this.nonlconConstraints
                disp([char(expr), ',']);
            end
            disp(']');
            
            lp4.Lp4Config.displayDelimiterLine();
            disp('C constraints:');
            
            disp('[');
            for expr = this.decvars(this.cStart : end)
                disp([char(expr), ',']);
            end
            disp(']');
            
            lp4.Lp4Config.displayDelimiterLine();
        end
        
    end % methods
    
    methods (Static)
        
        function lp = createLp(vars, f, eps, theta, psy, zeta, degree, pLambdaDegree, phyRange, pLambdaRange)
            lp =  lp4.LinearProgram4_v3(vars);
            lp = lp.initLp(f, eps, theta, psy, zeta, degree, pLambdaDegree, phyRange, pLambdaRange);
        end
        
        function lp = createLpWithoutRanges(vars, f, eps, theta, psy, zeta, degree, pLambdaDegree)
            lp = lp4.LinearProgram4_v3(vars);
            lp = lp.initLpWithoutRanges(f, eps, theta, psy, zeta, degree, pLambdaDegree);
        end
        
        function [wSymbolicVars, wExpression] = createWExpression(pPolynomial, pLambdaPolynomial)
            [wSymbolicVars, wExpression] = lp4.HybridLinearProgram.createWExpression('w', pPolynomial, pLambdaPolynomial);
        end % function createWExpression
        
    end % methods (Static)
    
end % classdef
