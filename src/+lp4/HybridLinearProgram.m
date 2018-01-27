classdef HybridLinearProgram
    
    properties
        indvars % 问题1中的独立变量，例如x = [x1, x2, ..., xn]，类型为符号化变量构成的行向量
        degree % 问题1中的带求多项式函数φ的次数，类型为正整数
        
        stateNum
        thetaStateIndex
        zetaStateIndex
        
        guardNum
        
        phys % 不同状态下的 phy
        phyPolynomials
        
        eps
        
        decvars % 问题1中的决策变量，包括 不同状态下Phy的决策变量, rou, lambda, re, wLambda, wRe, Cαβ, Cγδ, Cu,v
        decvarsIndexes
        
        exprs
        
        fs % 在不同状态下的 f
        
        guards
        
        theta
        psys
        zeta
        
        pLambdaDegree
        pLambdaPolynomials % 修改后的问题1中的λ
        
        pReDegree
        pRePolynomials
        
        wLambdas
        
        wRes
        
        rouVar
        
        pPartitions
        pLambdaPartitions
        pRePartitions
        
        linprogF
        
        rouInc  =  0
    end % properties
    
    methods
        
        function this = HybridLinearProgram(indvarsArg, stateNum, guardNum)
            % vars can only be a vector of a matrix of symbolic variables
            if ~isa(indvarsArg, 'sym')
                error('vars is not a vector of a matrix of symbolic variables');
            end
            if stateNum <= 0
                error('Wrong state num.');
            end
            if guardNum <= 0
                error('Wrong guard num.');
            end
            
            import lp4util.reshapeToVector
            this.indvars = reshapeToVector(indvarsArg);
            
            this.stateNum = stateNum;
            this.guardNum = guardNum;
            
            this.degree = 0;
            this.pLambdaDegree = 0;
            this.pReDegree = 0;
            
            this.eps = [];
            this.fs = [];
            
            this.exprs = [];
            
            import lp4.HybridLinearProgramDecVarIndexes
            this.decvarsIndexes = HybridLinearProgramDecVarIndexes();
        end
        
        function this = initWithoutRanges(this, fs, eps, thetaStateIndex, zetaStateIndex, theta, psys, zeta, guards, degree, pLambdaDegree, pReDegree)
            
            this.fs = fs;
            this.eps = eps;
            
            this.thetaStateIndex = thetaStateIndex;
            this.zetaStateIndex = zetaStateIndex;
            
            % Note: we needs guards when genereting dec vars
            this.guards = guards;
            this = this.setDegreeAndInit(degree, pLambdaDegree, pReDegree);
            
            this = this.setThetaConstraint(theta);
            this = this.setPsyConstraint(psys);
            this = this.setGuardConstraint(guards);
            this = this.setReConstraint();
            this = this.setZetaConstraint(zeta);
            
            this = this.generateEqsForSafetyConstraints();
            
            this = this.setDecVarsConstraint();
            
            this = this.setLinprogF();
        end
        
        function this = setRangeAndWConstraint(this, phyRange, pLambdaRange, pReRange)
            
            import lp4util.Partition
            this.pPartitions = repmat(phyRange, this.stateNum, 1024);
            this.pLambdaPartitions = repmat(pLambdaRange, this.stateNum, 1024);
            this.pRePartitions = repmat(pReRange, this.guardNum, 1024);
            this = this.setWConstraint();
        end
        
        function this = init(this, fs, eps, thetaStateIndex, zetaStateIndex, theta, psys, zeta, guards, degree, pLambdaDegree, pReDegree,...
                phyRange, pLambdaRange, pReRange)
            
            this = this.initWithoutRanges(fs, eps, thetaStateIndex, zetaStateIndex, theta, psys, zeta, guards, degree, pLambdaDegree, pReDegree);
            this = this.setRangeAndWConstraint(phyRange, pLambdaRange, pReRange);
        end
        

        
        function this = set.fs(this, fs)
            this.fs = fs;
        end
        
        function this = set.eps(this, eps)
            this.eps = eps;
        end
        
        function [this, startIndex, endIndex] = addDecisionVars(this, decvars)
            % addDecisionVars add decision variables.
            %
            % lp is a linear program.
            % decvars are decsion variables and should be symbolic.
            
            if isa(decvars, 'sym')
                % 将新的决策变量变为行向量的形式添加到lp的lp.decvars属性中
                % decision vars can only be a vector of a matrix of symbolic variables
                % reshape `decvars` to of dim `[1, size(decvars, 1) * size(decvars, 2)]`
                startIndex = size(this.decvars, 2) + 1;
                this.decvars = [this.decvars reshape(decvars, 1, size(decvars, 1) * size(decvars, 2))];
                endIndex = size(this.decvars, 2);
            else
                error("Wrong decvars.");
            end
        end
        
        function this = setDegreeAndInit(this, degree, pLambdaDegree, pReDegree)
            this.degree = degree;
            this.pLambdaDegree = pLambdaDegree;
            this.pReDegree = pReDegree;
            
            import lp4util.SymbolicPolynomial
            
            % set phy for each state
            this.phys = [];
            this.phyPolynomials = [];
            for i = 1 : this.stateNum
                currentStateP = sym(strcat('p', num2str(i), '_'), [1, monomialNumber(length(this.indvars), degree)]);
                [this, this.decvarsIndexes.phyStarts(i), this.decvarsIndexes.phyEnds(i)] = this.addDecisionVars(currentStateP);
                
                currentStatePhy = currentStateP * monomials(this.indvars, 0 : degree);
                this.phys = [this.phys, currentStatePhy];
                this.phyPolynomials = [this.phyPolynomials, SymbolicPolynomial(this.indvars, degree, currentStateP, currentStatePhy)];
            end
            
            this.rouVar = sym('rou');
            [this, this.decvarsIndexes.rouIndex, ~] = this.addDecisionVars(this.rouVar);
            
            % set lambda for each state
            this.pLambdaPolynomials = [];
            for i = 1 : this.stateNum
                currentStatePLambda = sym(strcat('pLmabda', num2str(i), '_'), [1, monomialNumber(length(this.indvars), pLambdaDegree)]);
                [this, this.decvarsIndexes.lambdaStarts(i), this.decvarsIndexes.lambdaEnds(i)] = this.addDecisionVars(currentStatePLambda);
                
                currentStatePLambdaExpr = currentStatePLambda * monomials(this.indvars, 0 : pLambdaDegree);
                this.pLambdaPolynomials = [this.pLambdaPolynomials, SymbolicPolynomial(this.indvars, pLambdaDegree, currentStatePLambda, currentStatePLambdaExpr)];
            end
            
            % set re for each guard
            this.pRePolynomials = [];
            for i = 1 : this.guardNum
                currentStatePRe = sym(strcat('pRe', num2str(i), '_'), [1, monomialNumber(length(this.indvars), pReDegree)]);
                [this, this.decvarsIndexes.reStarts(i), this.decvarsIndexes.reEnds(i)] = this.addDecisionVars(currentStatePRe);
                
                currentStatePReExpr = currentStatePRe * monomials(this.indvars, 0 : pReDegree);
                this.pRePolynomials = [this.pRePolynomials, SymbolicPolynomial(this.indvars, pReDegree, currentStatePRe, currentStatePReExpr)];
            end
            
            % introduce WLambda = P * PLambda for each state in formula 2
            this.wLambdas = [];
            import lp4.HybridLinearProgram
            for i = 1 : this.stateNum
                name = strcat('wLambda', num2str(i), '_');
                [currentWLambdaSymbolicVars, currentWLambdaExpression] = HybridLinearProgram.createWExpression(name, this.phyPolynomials(i), this.pLambdaPolynomials(i));
                
                import lp4.DecVarsAndExpr
                this.wLambdas = [this.wLambdas, DecVarsAndExpr(currentWLambdaSymbolicVars, currentWLambdaExpression)];
                
                [this, this.decvarsIndexes.wLambdaStarts(i), this.decvarsIndexes.wLambdaEnds(i)] = this.addDecisionVars(currentWLambdaSymbolicVars);
            end
            
            % introduce WRe = P * PRe for each guard in formula 3
            this.wRes = [];
            for i = 1 : this.guardNum
                guard = this.guards(i);
                guardFrom = guard.fromStateIndex;
                guardDest = guard.destStateIndex;
                name = strcat('wRe', num2str(i), '_');
                
                [currentWReSymbolicVars, currentWReExpression] = HybridLinearProgram.createWExpression(name, this.phyPolynomials(guardFrom), this.pRePolynomials(i));
                
                import lp4.DecVarsAndExpr
                this.wRes = [this.wRes, DecVarsAndExpr(currentWReSymbolicVars, currentWReExpression)];
                
                [this, this.decvarsIndexes.wReStarts(i), this.decvarsIndexes.wReEnds(i)] = this.addDecisionVars(currentWReSymbolicVars);
            end
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
            
            % different from LP, the phy is the phy of theta state
            thetaStatePhy = this.phys(this.thetaStateIndex);
            constraint1 = -thetaStatePhy;
            import lp4.Lp4Config
            de = computeDegree(constraint1, this.indvars) + Lp4Config.C_DEGREE_INC;
            
            c_alpha_beta = sym('c_alpha_beta', [1,10000]); % pre-defined varibales, only a few of them are the actual variables
            [constraintDecvars, expression] = constraintExpression(de, theta, c_alpha_beta);
            [this, this.decvarsIndexes.cThetaStart, this.decvarsIndexes.cThetaEnd] = this.addDecisionVars(constraintDecvars);
            
            constraint1 = constraint1 + expression;
            
            constraint1 = expand(constraint1);
            expr.polyexpr = constraint1;
            
            % Note: This is the first expression, and use `this.exprs(expr.num) = expr;` here
            % will cause an error in Matlab 2014b but is Ok in Matlab 2017b.
            this.exprs = [this.exprs, expr];
        end
        
        function this = setPsyConstraint(this, psys)
            if this.stateNum ~= size(psys, 1)
                error("Psys are of wrong number.")
            end
            
            this.psys = psys;
            
            for i = 1 : this.stateNum
                psy = psys(1,:);
                exprNum = this.nextExprNumIndex();
                
                % for an empty constrain
                if isempty(psy)
                    expr = Constraint.createEmptyConstraint();
                    expr.num = exprNum;
                    this.exprs(exprNum) = expr;
                    continue;
                end
                
                expr = Constraint();
                expr.num = exprNum;
                expr.name = 'psy';
                expr.type = 'eq';
                expr.A = [];
                expr.b = [];
                
                currentPhy = this.phys(i);
                currentF = this.fs(:, i);
                phy_d = 0;
                for k = 1 : 1 : length(this.indvars)
                    phy_d = phy_d + diff(currentPhy, this.indvars(k)) * currentF(k);
                end
                
                constraint2 = -phy_d + this.wLambdaExpressions(i) + this.eps(1);
                import lp4.Lp4Config
                de = computeDegree(constraint2, this.indvars) + Lp4Config.C_DEGREE_INC;
                
                name = strcat('c_gama_delta', num2str(i), '_');
                c_gama_delta = sym(name, [1, 1000 * 1000]);
                [constraintDecvars, expression] = constraintExpression(de, psy, c_gama_delta);
                [this, this.decvarsIndexes.cPsyStarts(i), this.decvarsIndexes.cPsyEnds(i)] = this.addDecisionVars(constraintDecvars);
                
                constraint2 = constraint2 + expression;
                
                constraint2 = expand(constraint2);
                expr.polyexpr = constraint2;
                
                this.exprs(exprNum) = expr;
            end  
        end
        
        function this = setGuardConstraint(this, guards)
            this.guards = guards;
            this.guardNum = length(guards);
            
            for i = 1 : this.guardNum
                guard = this.guards(i);
                guardFrom = guard.fromStateIndex;
                guardDest = guard.destStateIndex;
                
                exprNum = this.nextExprNumIndex();
                
                % for an empty constrain
                if isempty(guard.exprs)
                    expr = Constraint.createEmptyConstraint();
                    expr.num = exprNum;
                    this.exprs(exprNum) = expr;
                    continue;
                end
                
                expr = Constraint();
                expr.num = exprNum;
                expr.name = 'guard';
                expr.type = 'eq';
                expr.A = [];
                expr.b = [];
                
                guardFromPhy = this.phys(guardFrom);
                guardDestPhy = this.phys(guardDest);
                
                constraintGuard = -guardDestPhy;
                for i1 = 1 : length(guard.resetVars)
                    resetVar = guard.resetVars(i1);
                    resetVarValue = guard.resetVarValues(i1);
                    constraintGuard = subs(constraintGuard, resetVar, resetVarValue);
                end
                
                constraintGuard = constraintGuard + this.wReExpressions(i);
                
                import lp4.Lp4Config
                de = computeDegree(constraintGuard, this.indvars) + Lp4Config.C_DEGREE_INC;
                
                name = strcat('c_guard', num2str(i), '_');
                c_guard = sym(name, [1, 1000 * 1000]);
                [constraintDecvars, expression] = constraintExpression(de, guard.exprs, c_guard);
                [this, this.decvarsIndexes.cGuardStarts(i), this.decvarsIndexes.cGuardEnds(i)] = this.addDecisionVars(constraintDecvars);
                
                constraintGuard = constraintGuard + expression;
                
                constraintGuard = expand(constraintGuard);
                expr.polyexpr = constraintGuard;
                
                this.exprs(exprNum) = expr;
            end
        end
        
        function this = setReConstraint(this)
            
            for i = 1 : this.guardNum
                guard = this.guards(i);
                
                exprNum = this.nextExprNumIndex();
                
                expr = Constraint();
                expr.num = exprNum;
                expr.name = 're';
                expr.type = 'eq';
                expr.A = [];
                expr.b = [];
                
                constraintRe = - (this.pRePolynomials(i).expression);
                
                import lp4.Lp4Config
                de = computeDegree(constraintRe, this.indvars) + Lp4Config.C_DEGREE_INC;
                
                name = strcat('c_re', num2str(i), '_');
                c_re = sym(name, [1, 1000 * 1000]);
                [constraintDecvars, expression] = constraintExpression(de, guard.exprs, c_re);
                [this, this.decvarsIndexes.cReStarts(i), this.decvarsIndexes.cReEnds(i)] = this.addDecisionVars(constraintDecvars);
                
                constraintRe = constraintRe + expression;
                
                constraintRe = expand(constraintRe);
                expr.polyexpr = constraintRe;
                
                this.exprs(exprNum) = expr;
            end
        end
        
        function this = setZetaConstraint(this, zeta)
            this.zeta = zeta;
            
            exprNum = this.nextExprNumIndex();
            
            % for an empty constrain
            if isempty(zeta)
                expr = Constraint.createEmptyConstraint();
                expr.num = exprNum;
                this.exprs(exprNum) = expr;
                return;
            end
            
            expr = Constraint();
            expr.num = exprNum;
            expr.name = 'zeta';
            expr.type = 'eq';
            expr.A = [];
            expr.b = [];
            
            zetaStatePhy = this.phys(this.zetaStateIndex);
            constraint3 = zetaStatePhy + this.eps(2); % different from lp2
            import lp4.Lp4Config
            de = computeDegree(constraint3, this.indvars) + Lp4Config.C_DEGREE_INC;
            
            c_u_v = sym('c_u_v',[1,10000]);
            [constraintDecvars, expression] = constraintExpression(de, zeta, c_u_v);
            [this, this.decvarsIndexes.cZetaStart, this.decvarsIndexes.cZetaEnd] = this.addDecisionVars(constraintDecvars);
            
            constraint3 = constraint3 + expression;
            
            constraint3 = expand(constraint3);
            expr.polyexpr = constraint3;
            
            this.exprs(exprNum) = expr;
        end
        
        function this = generateEqsForSafetyConstraints(this)
            for k = 1 : length(this.exprs)
                % skip empty constraint
                if this.exprs(k).isEmptyConstraint() || this.exprs(k).isEqGenerated()
                    continue;
                end
                
                [ this.exprs(k).A, this.exprs(k).b ] = eqgenerate( this.indvars, this.decvars, this.exprs(k).polyexpr);
                disp(['constraint ',this.exprs(k).name,' is processed: ',datestr(now,'yyyy-mm-dd HH:MM:SS')]);
            end
        end
        
        function this = setDecVarsConstraint(this)
            exprNum = this.nextExprNumIndex();
            
            expr = Constraint();
            expr.num = exprNum;
            expr.name = 'decvarConstraints';
            expr.type = 'ie';
            expr.polyexpr = [];
            
            cStart = this.decvarsIndexes.cThetaStart;
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
            
            this.exprs(exprNum) = expr;
        end % function setDecVarsConstraint
        
        function this = setWConstraint(this)
            lastExpr = this.exprs(length(this.exprs));
            if (strcmp(lastExpr.name, 'wConstraints'))
                % if the last expr is for w, renew this expr
                exprNum = length(this.exprs);
            else 
                exprNum = this.nextExprNumIndex();
            end
            
            expr = Constraint();
            expr.num = exprNum;
            expr.name = 'wConstraints';
            expr.type = 'ie';
            expr.A = [];
            expr.b = [];
            
            rouIndex = this.getRouIndex();
            for j = 1 : this.stateNum
                phyPolynomial = this.phyPolynomials(j);
                
                for i = 1 : length(phyPolynomial.coefficientVars)
                    p = phyPolynomial.coefficientVars(i);
                    pPartition = this.pPartitions(j, i);
                    
                    expr = pPartition.createConstraintsWithRouAndAddToExpr(p, this.decvars, rouIndex, expr);
                end
            end
            
            for j = 1 : this.stateNum
                pLambdaPolynomial = this.pLambdaPolynomials(j);
                
                for i = 1 : length(pLambdaPolynomial.coefficientVars)
                    p = pLambdaPolynomial.coefficientVars(i);
                    pPartition = this.pLambdaPartitions(j, i);
                    
                    expr = pPartition.createConstraintsWithRouAndAddToExpr(p, this.decvars, rouIndex, expr);
                end
            end
            
            for j = 1 : this.stateNum
                phyPolynomial = this.phyPolynomials(j);
                pLambdaPolynomial = this.pLambdaPolynomials(j);
                wSymbolicVars = this.wLambdas(j).symbolicVars;
                
                for i1 = 1 : 1 : length(phyPolynomial.coefficientVars)
                    for i2 = 1 : 1 : length(pLambdaPolynomial.coefficientVars)
                        p = phyPolynomial.coefficientVars(i1);
                        pPartition = this.pPartitions(j, i1);
                        pLambda = pLambdaPolynomial.coefficientVars(i2);
                        pLambdaPartition = this.pLambdaPartitions(j, i2);
                        w = wSymbolicVars(i1, i2);
                        
                        import lp4.createOverApproximationAsConstraintsWithRou
                        constraints = createOverApproximationAsConstraintsWithRou(p, pPartition, pLambda, pLambdaPartition, w, this.rouVar);
                        
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
            end
            
            for j = 1 : this.guardNum
                pRePolynomial = this.pRePolynomials(j);
                
                for i = 1 : length(pRePolynomial.coefficientVars)
                    p = pRePolynomial.coefficientVars(i);
                    pPartition = this.pRePartitions(j, i);
                    
                    expr = pPartition.createConstraintsWithRouAndAddToExpr(p, this.decvars, rouIndex, expr);
                end
            end
            
            for j = 1 : this.guardNum
                guard = this.guards(j);
                guardFrom = guard.fromStateIndex;
                phyFromPolynomial = this.phyPolynomials(guardFrom);
                pRePolynomial = this.pRePolynomials(j);
                wReSymbolicVars = this.wRes(j).symbolicVars;
                
                for i1 = 1 : 1 : length(phyFromPolynomial.coefficientVars)
                    for i2 = 1 : 1 : length(pRePolynomial.coefficientVars)
                        p = phyFromPolynomial.coefficientVars(i1);
                        pPartition = this.pPartitions(guardFrom, i1);
                        pRe = pRePolynomial.coefficientVars(i2);
                        pRePartition = this.pRePartitions(j, i2);
                        w = wReSymbolicVars(i1, i2);
                        
                        import lp4.createOverApproximationAsConstraintsWithRou
                        constraints = createOverApproximationAsConstraintsWithRou(p, pPartition, pRe, pRePartition, w, this.rouVar);
                        
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
            end
            
            this.exprs(exprNum) = expr;
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
            
            import lp4.HybridLinearProgramSolveRes
            solveRes = HybridLinearProgramSolveRes(this, x, fval, flag, time);
            
            import lp4.Lp4Config
            if Lp4Config.isDebug()
                solveRes.printSolution();
            end
            
        end % function solve
        
        % FIXME `verifyWithLambda`
        
        % `verify` is forward to `verifyWithLambda`
        function [lpVer, solveResVer, resNorms] = verify(lp, solveRes, phyRangeInVerify)
            [lpVer, solveResVer, resNorms] = verifyWithLambda(lp, solveRes, phyRangeInVerify);
        end
        
        function [lpVer, solveResVer, resNorms] = verifyWithLambda(lp, solveRes, phyRangeInVerify)
            
            if ~(solveRes.hasSolution())
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
            
            lpVer = lpVer.setDecVarsConstraint();
            
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
                disp('Not find feasible solution to verify.');
                return;
            end
            
            % verify
            
            import lp4.HybridLinearProgramVerificationWithGivenPhy
            lpVer = HybridLinearProgramVerificationWithGivenPhy.create(lp.indvars, lp.stateNum,...
                lp.fs, lp.eps, lp.thetaStateIndex, lp.zetaStateIndex, lp.theta, lp.psys, lp.zeta, lp.guards, lp.pLambdaDegree, lp.pReDegree,...
                solveRes.getPhyExpressions());
            
            % FIXME `resNorms` is empty
            [lpVer, solveResVer, resNorms] = lpVer.solve();
            
            if ~(solveResVer.hasSolution())
                disp('Verify feasible solution failed.');
            else
                disp('Verify feasible solution succeed, norms ;');
                % FIXME
                % disp(resNorms);
            end
            
            import lp4.Lp4Config
            Lp4Config.displayDelimiterLine();
        end
        
        function res = getPhyCoefficientStart(this, i)
            res = this.decvarsIndexes.phyStarts(i);
        end
        
        function res = getPhyCoefficientEnd(this, i)
            res = this.decvarsIndexes.phyEnds(i);
        end
        
        function res = getRouIndex(this)
            res = this.decvarsIndexes.rouIndex;
        end
        
        function res = getPLambdaCoefficientStart(this, i)
            res = this.decvarsIndexes.lambdaStarts(i);
        end
        
        function res = getPLambdaCoefficientEnd(this, i)
            res = this.decvarsIndexes.lambdaEnds(i);
        end
        
        function res = getPReCoefficientStart(this, i)
            res = this.decvarsIndexes.reStarts(i);
        end
        
        function res = getPReCoefficientEnd(this, i)
            res = this.decvarsIndexes.reEnds(i);
        end
        
        function res = getWLambdaStart(this, i)
            res = this.decvarsIndexes.wLambdaStarts(i);
        end
        
        function res = getWLambdaEnd(this, i)
            res = this.decvarsIndexes.wLambdaEnds(i);
        end
        
        function res = getWReStart(this, i)
            res = this.decvarsIndexes.wReStarts(i);
        end
        
        function res = getWReEnd(this, i)
            res = this.decvarsIndexes.wReEnds(i);
        end
        
        function res = nextExprNumIndex(this)
            res = length(this.exprs) + 1;
        end
        
        function res = wLambdaExpressions(this, i)
            res = this.wLambdas(i).expression;
        end
        
        function res = wReExpressions(this, i)
            res = this.wRes(i).expression;
        end
        
    end % methods
    
    methods (Static)
        
        function hlp = create(indvarsArg, stateNum,...
                fs, eps, thetaStateIndex, zetaStateIndex, theta, psys, zeta,guards, degree, pLambdaDegree, pReDegree,...
                phyRange, pLambdaRange, pReRange)
            import lp4.HybridLinearProgram
            hlp = HybridLinearProgram(indvarsArg, stateNum, size(guards, 2));
            hlp = hlp.init(fs, eps, thetaStateIndex, zetaStateIndex, theta, psys, zeta,guards, degree, pLambdaDegree, pReDegree, phyRange, pLambdaRange, pReRange);
        end
        
        function hlp = createWithoutRanges(indvarsArg, stateNum,...
                fs, eps, thetaStateIndex, zetaStateIndex, theta, psys, zeta, guards, degree, pLambdaDegree, pReDegree)
            import lp4.HybridLinearProgram
            hlp = HybridLinearProgram(indvarsArg, stateNum, size(guards, 2));
            hlp = hlp.initWithoutRanges(fs, eps, thetaStateIndex, zetaStateIndex, theta, psys, zeta, guards, degree, pLambdaDegree, pReDegree);
        end
        
        function [wSymbolicVars, wExpression] = createWExpression(name, pPolynomial, pLambdaPolynomial)
            wSymbolicVars = sym(name, [length(pPolynomial.coefficientVars), length(pLambdaPolynomial.coefficientVars)]);
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
    
end
