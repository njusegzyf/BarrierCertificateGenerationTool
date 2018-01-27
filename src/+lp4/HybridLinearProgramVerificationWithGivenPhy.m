classdef HybridLinearProgramVerificationWithGivenPhy
    
    properties
        indvars % 问题1中的独立变量，例如x = [x1, x2, ..., xn]，类型为符号化变量构成的行向量
        
        stateNum
        thetaStateIndex
        zetaStateIndex
        
        guardNum
        
        phys % 不同状态下的 phy
        
        eps
        
        decvars % 问题1中的决策变量，包括 lambda, re, Cαβ, Cγδ, Cu,v
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
        
        pPartitions
        pLambdaPartitions
        pRePartitions
        
        linprogF
    end % properties
    
    methods
        
        function this = HybridLinearProgramVerificationWithGivenPhy(indvarsArg, stateNum, guardNum)
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
            
            this.pLambdaDegree = 0;
            this.pReDegree = 0;
            
            this.eps = [];
            this.fs = [];
            
            this.exprs = [];
            
            import lp4.HybridLinearProgramDecVarIndexes
            this.decvarsIndexes = HybridLinearProgramDecVarIndexes();
        end
        
        function this = initWithoutRanges(this, fs, eps, thetaStateIndex, zetaStateIndex, theta, psys, zeta, guards, pLambdaDegree, pReDegree, phys)
            
            this.fs = fs;
            this.eps = eps;
            
            this.thetaStateIndex = thetaStateIndex;
            this.zetaStateIndex = zetaStateIndex;
            
            this.phys = phys;
            
            % Note: we needs guards when genereting dec vars
            this.guards = guards;
            this = this.setDegreeAndInit(pLambdaDegree, pReDegree);
            
            this = this.setThetaConstraint(theta);
            this = this.setPsyConstraint(psys);
            this = this.setGuardConstraint(guards);
            this = this.setReConstraint();
            this = this.setZetaConstraint(zeta);
            
            this = this.generateEqsForSafetyConstraints();
            
            this = this.setDecVarsConstraint();
            
            this = this.setLinprogF();
        end
        
        function this = init(this, fs, eps, thetaStateIndex, zetaStateIndex, theta, psys, zeta, guards, pLambdaDegree, pReDegree,...
                phyRange, pLambdaRange, pReRange)
            
            this = this.initWithoutRanges(fs, eps, thetaStateIndex, zetaStateIndex, theta, psys, zeta, guards, pLambdaDegree, pReDegree);
            
            import lp4util.Partition
            this.pPartitions = repmat(phyRange, this.stateNum, 1024);
            this.pLambdaPartitions = repmat(pLambdaRange, this.stateNum, 1024);
            this.pRePartitions = repmat(pReRange, this.guardNum, 1024);
            this = this.setWConstraint();
        end
        
        function this = set.phys(this, phys)
            this.phys = phys;
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
        
        function this = setDegreeAndInit(this, pLambdaDegree, pReDegree)
            this.pLambdaDegree = pLambdaDegree;
            this.pReDegree = pReDegree;
            
            import lp4util.SymbolicPolynomial
            
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
                
                constraint2 = -phy_d + this.pLambdaPolynomials(i).expression * this.phys(i) + this.eps(1);
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
                    constraintGuard = subs(constraintGuard, guard.resetVars(i1), guard.resetVarValues(i1));
                end
                
                constraintGuard = constraintGuard + this.pRePolynomials(i).expression * guardFromPhy;
                
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
            expr.name = 'decvarconstraints';
            expr.type = 'ie';
            expr.polyexpr = [];
            
            cStart = this.decvarsIndexes.cThetaStart;
            cLength = length(this.decvars) - cStart;
            expr.A = zeros(cLength, length(this.decvars));

            for k = 1 : 1 : cLength
                expr.A(k, cStart + k) = -1;
            end
            bc = zeros(cLength, 1);
            expr.b = bc;
            
            this.exprs(exprNum) = expr;
        end % function setDecVarsConstraint
        
        function this = setWConstraint(this)
            exprNum = this.nextExprNumIndex();
            
            expr = Constraint();
            expr.num = exprNum;
            expr.name = 'w';
            expr.type = 'ie';
            expr.A = [];
            expr.b = [];
            
            for j = 1 : this.stateNum
                pLambdaPolynomial = this.pLambdaPolynomials(j);
                
                for i = 1 : length(pLambdaPolynomial.coefficientVars)
                    p = pLambdaPolynomial.coefficientVars(i);
                    pPartition = this.pLambdaPartitions(j, i);
                    
                    expr = pPartition.createConstraintsAndAddToExpr(p, this.decvars, expr);
                end
            end
            
            for j = 1 : this.guardNum
                pRePolynomial = this.pRePolynomials(j);
                
                for i = 1 : length(pRePolynomial.coefficientVars)
                    p = pRePolynomial.coefficientVars(i);
                    pPartition = this.pRePartitions(j, i);
                    
                    expr = pPartition.createConstraintsAndAddToExpr(p, this.decvars, expr);
                end
            end
            
            this.exprs(exprNum) = expr;
        end
        
        function this = setLinprogF(this)
            % init f with all zeros
            this.linprogF = zeros(1, length(this.decvars));
        end
        
        function [this, solveRes, resNorms] = solve(this)
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
            
            import lp4.HybridLinearProgramVerificationWithGivenPhySolveRes
            solveRes = HybridLinearProgramVerificationWithGivenPhySolveRes(this, x, fval, flag, time);
            
            import lp4.Lp4Config
            if Lp4Config.isDebug()
                solveRes.printSolution();
            end
            
            % FIXME
            resNorms = [];
            
        end % function solve
        
        function res = getPhyCoefficientStart(this, i)
            res = this.decvarsIndexes.phyStarts(i);
        end
        
        function res = getPhyCoefficientEnd(this, i)
            res = this.decvarsIndexes.phyEnds(i);
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
                fs, eps, thetaStateIndex, endStateIndex, theta, psys, zeta, guards, pLambdaDegree, pReDegree,...
                phys)
            import lp4.HybridLinearProgramVerificationWithGivenPhy
            hlp = HybridLinearProgramVerificationWithGivenPhy(indvarsArg, stateNum, size(guards, 2));
            hlp = hlp.initWithoutRanges(fs, eps, thetaStateIndex, endStateIndex, theta, psys, zeta,guards, pLambdaDegree, pReDegree, phys);
            
        end
        
    end % methods (Static)
    
end
