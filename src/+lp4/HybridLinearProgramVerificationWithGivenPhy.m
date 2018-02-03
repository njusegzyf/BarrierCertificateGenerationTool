classdef HybridLinearProgramVerificationWithGivenPhy < lp4.HybridLinearProgramVerificationBase
    
    properties
        
        pLambdaDegree
        pLambdaPolynomials % 修改后的问题1中的λ
        
        pReDegree
        pRePolynomials
        
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
        
        function this = init(this, fs, eps, thetaStateIndex, theta, psys, zetas, guards, pLambdaDegree, pReDegree, phys)
            
            this.fs = fs;
            this.eps = eps;
            
            this.thetaStateIndex = thetaStateIndex;
            
            this.phys = phys;
            
            % Note: we needs guards when genereting dec vars
            this.guards = guards;
            this = this.setDegreeAndInit(pLambdaDegree, pReDegree);
            
            this = this.setThetaConstraint(theta);
            this = this.setPsyConstraint(psys);
            this = this.setGuardConstraint(guards);
            this = this.setReConstraint();
            this = this.setZetaConstraint(zetas);
            
            this = this.generateEqsForSafetyConstraints();
            
            this = this.setDecVarsConstraint();
            
            this = this.setLinprogF();
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
            
            if this.isAttachRou
                this.rouVar = sym('rou');
                [this, this.decvarsIndexes.rouIndex, ~] = this.addDecisionVars(this.rouVar);
            end
        end
        
        function this = setThetaConstraint(this, theta)
            this.theta = theta;
            
            % for an empty constrain
            if isempty(theta)
                expr = Constraint.createEmptyConstraint();
                expr.num = 1;
                this.exprs = [this.exprs, expr];
                p    return;
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
            
            c_alpha_beta = sym('c_alpha_beta', [1, lp4.Lp4Config.DEFAULT_DEC_VAR_SIZE]); % pre-defined varibales, only a few of them are the actual variables
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
                psy = psys(i, :);
                exprNum = this.nextExprNumIndex();
                
                % for an empty constrain
                if isempty(psy)
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
                c_gama_delta = sym(name, [1, lp4.Lp4Config.DEFAULT_DEC_VAR_SIZE]);
                [constraintDecvars, expression] = constraintExpression(de, psy, c_gama_delta);
                [this, this.decvarsIndexes.cPsyStarts(i), this.decvarsIndexes.cPsyEnds(i)] = this.addDecisionVars(constraintDecvars);
                
                constraint2 = constraint2 + expression;
                
                constraint2 = expand(constraint2);
                expr.polyexpr = constraint2;
                
                this.exprs(exprNum) = expr;
            end
            
        end
        
        function this = setGuardConstraint(this, guards)
            if length(this.guards) ~= this.guardNum
                error("Wrong guard number.")
            end
            
            this.guards = guards;
            
            for i = 1 : this.guardNum
                guard = this.guards(i);
                guardFrom = guard.fromStateIndex;
                guardDest = guard.destStateIndex;
                
                exprNum = this.nextExprNumIndex();
                
                % for an empty constrain
                if isempty(guard.exprs)
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
                
                constraintGuard = constraintGuard + this.pRePolynomials(i).expression * guardFromPhy;
                
                import lp4.Lp4Config
                de = computeDegree(constraintGuard, this.indvars) + Lp4Config.C_DEGREE_INC;
                
                name = strcat('c_guard', num2str(i), '_');
                c_guard = sym(name, [1, lp4.Lp4Config.DEFAULT_DEC_VAR_SIZE]);
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
                
                % for an empty constrain
                if isempty(guard.exprs)
                    continue;
                end
                
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
                c_re = sym(name, [1, lp4.Lp4Config.DEFAULT_DEC_VAR_SIZE]);
                [constraintDecvars, expression] = constraintExpression(de, guard.exprs, c_re);
                [this, this.decvarsIndexes.cReStarts(i), this.decvarsIndexes.cReEnds(i)] = this.addDecisionVars(constraintDecvars);
                
                constraintRe = constraintRe + expression;
                
                constraintRe = expand(constraintRe);
                expr.polyexpr = constraintRe;
                
                this.exprs(exprNum) = expr;
            end
        end
        
        function this = setZetaConstraint(this, zetas)
            this.zetas = zetas;
            
            for zetaIndex = 1 : length(this.zetas)
                zeta = this.zetas(zetaIndex);
                
                % for an empty constrain
                if isempty(zeta.exprs)
                    return;
                end
                
                exprNum = this.nextExprNumIndex();
                
                expr = Constraint();
                expr.num = exprNum;
                expr.name = 'zeta';
                expr.type = 'eq';
                expr.A = [];
                expr.b = [];
                
                zetaStatePhy = this.phys(zeta.stateIndex);
                constraint3 = zetaStatePhy + this.eps(2); % different from lp2
                import lp4.Lp4Config
                de = computeDegree(constraint3, this.indvars) + Lp4Config.C_DEGREE_INC;
                
                name = strcat('c_u_v', num2str(zetaIndex));
                c_u_v = sym(name,[1, lp4.Lp4Config.DEFAULT_DEC_VAR_SIZE]);
                [constraintDecvars, expression] = constraintExpression(de, zeta.exprs, c_u_v);
                [this, this.decvarsIndexes.cZetaStarts(zetaIndex), this.decvarsIndexes.cZetaEnds(zetaIndex)] = this.addDecisionVars(constraintDecvars);
                
                constraint3 = constraint3 + expression;
                
                constraint3 = expand(constraint3);
                expr.polyexpr = constraint3;
                
                this.exprs(exprNum) = expr;
            end
        end
        
        function this = generateEqsForSafetyConstraints(this)
            this = lp4.Lp4AndHlpVerificationBase.generateConstraintEqsParallelly(this);
        end
        
        function this = setDecVarsConstraint(this)
            exprNum = this.nextExprNumIndex();
            
            expr = Constraint();
            expr.num = exprNum;
            expr.name = 'decvarconstraints';
            expr.type = 'ie';
            expr.polyexpr = [];
            
            cStart = this.decvarsIndexes.cThetaStart - 1;  % Note: `cStart + 1` is the actual start index
            cLength = length(this.decvars) - cStart;
            expr.A = zeros(cLength, length(this.decvars));
            
            for k = 1 : 1 : cLength
                expr.A(k, cStart + k) = -1;
            end
            if this.isAttachRou
                rouIndex = this.decvarsIndexes.rouIndex;
                for k = 1 : cLength
                    % - rou
                    expr.A(k, rouIndex) = -1;
                end
            end
            expr.b = zeros(cLength, 1);
            
            this.exprs(exprNum) = expr;
        end % function setDecVarsConstraint
        
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
        
        function res = nextExprNumIndex(this)
            res = length(this.exprs) + 1;
        end
        
        function solveRes = createSolveRes(this, x, fval, flag, time)
            solveRes = lp4.HybridLinearProgramVerificationWithGivenPhySolveRes(this, x, fval, flag, time);
        end
        
    end % methods
    
    methods (Static)
        
        function hlp = create(indvarsArg, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards, pLambdaDegree, pReDegree,...
                phys)
            
            import lp4.HybridLinearProgramVerificationWithGivenPhy
            hlp = HybridLinearProgramVerificationWithGivenPhy(indvarsArg, stateNum, size(guards, 2));
            hlp = hlp.init(fs, eps, thetaStateIndex, theta, psys, zetas, guards, pLambdaDegree, pReDegree, phys);
        end
        
        function hlp = createWithRou(indvarsArg, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards, pLambdaDegree, pReDegree,...
                phys)
            
            import lp4.HybridLinearProgramVerificationWithGivenPhy
            hlp = HybridLinearProgramVerificationWithGivenPhy(indvarsArg, stateNum, size(guards, 2));
            hlp.isAttachRou = true;
            hlp = hlp.init(fs, eps, thetaStateIndex, theta, psys, zetas, guards, pLambdaDegree, pReDegree, phys);
        end
        
    end % methods (Static)
    
end
