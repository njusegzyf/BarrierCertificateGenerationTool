classdef HybridLinearProgramVerificationWithGivenLambdaAndRe
    
    properties
        indvars % 问题1中的独立变量，例如x = [x1, x2, ..., xn]，类型为符号化变量构成的行向量
        degree % 问题1中的带求多项式函数φ的次数，类型为正整数
        
        stateNum
        thetaStateIndex
        
        guardNum
        
        phys % 不同状态下的 phy
        phyPolynomials
        
        eps
        
        lambdas
        res
        
        fs % 在不同状态下的 f
        
        theta
        psys
        zetas
        guards
        
        decvars % 问题1中的决策变量，包括 不同状态下Phy的决策变量, Cαβ, Cγδ, Cu,v
        decvarsIndexes
        
        exprs
        
        linprogF

    end % properties
    
    methods
        
        function this = HybridLinearProgramVerificationWithGivenLambdaAndRe(indvarsArg, stateNum, guardNum)
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
            
            this.eps = [];
            this.fs = [];
            
            this.exprs = [];
            
            import lp4.HybridLinearProgramDecVarIndexes
            this.decvarsIndexes = HybridLinearProgramDecVarIndexes();
        end
        
        function this = init(this, fs, eps, thetaStateIndex, theta, psys, zetas, guards, degree,...
                             lambdas, res)
            
            this.fs = fs;
            this.eps = eps;
            
            this.thetaStateIndex = thetaStateIndex;
           
            this.lambdas = lambdas;
            this.res = res;
            
            % Note: we needs guards when genereting dec vars
            this.guards = guards;
            this = this.setDegreeAndInit(degree);
            
            this = this.setThetaConstraint(theta);
            this = this.setPsyConstraint(psys);
            this = this.setGuardConstraint(guards);
            this = this.setReConstraint();
            this = this.setZetaConstraint(zetas);
            
            this = this.generateEqsForSafetyConstraints();
            
            this = this.setDecVarsConstraint();
            
            this = this.setLinprogF();
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
                error('Wrong decvars.');
            end
        end
        
        function this = setDegreeAndInit(this, degree)
            this.degree = degree;
            
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
            
            import lp4.Lp4Config
            c_alpha_beta = sym('c_alpha_beta', [1, Lp4Config.DEFAULT_DEC_VAR_SIZE]); % pre-defined varibales, only a few of them are the actual variables
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
                error('Psys are of wrong number.')
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
                
                constraint2 = -phy_d + this.lambdas(i) * currentPhy + this.eps(1);
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
                error('Wrong guard number.')
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
                
                constraintGuard = constraintGuard + this.res(i) * guardFromPhy;
                
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
                
                constraintRe = - this.res(i);
                
                import lp4.Lp4Config
                de = computeDegree(constraintRe, this.indvars) + Lp4Config.C_DEGREE_INC;
                
                % if re is a constant, ignore it
                if de == 0
                    continue;
                end
                
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
            
            cStart = this.decvarsIndexes.cThetaStart - 1;  % Note: `cStart + 1` is the actual start index
            cLength = length(this.decvars) - cStart;
            expr.A = zeros(cLength, length(this.decvars));

            for k = 1 : 1 : cLength
                expr.A(k, cStart + k) = -1;
            end
            expr.b = zeros(cLength, 1);
            
            this.exprs(exprNum) = expr;
        end % function setDecVarsConstraint

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
            
            import lp4.HybridLinearProgramVerificationWithGivenLambdaAndReSolveRes
            solveRes = HybridLinearProgramVerificationWithGivenLambdaAndReSolveRes(this, x, fval, flag, time);
            
            import lp4.Lp4Config
            if Lp4Config.isDebug()
                solveRes.printSolution();
            end
            
            if solveRes.hasSolution()
                resNorms = solveRes.computeAllExprsNorms();
            else
                resNorms = [];
            end
            
        end % function solve
        
        function res = getPhyCoefficientStart(this, i)
            res = this.decvarsIndexes.phyStarts(i);
        end
        
        function res = getPhyCoefficientEnd(this, i)
            res = this.decvarsIndexes.phyEnds(i);
        end
        
        function res = nextExprNumIndex(this)
            res = length(this.exprs) + 1;
        end
      
    end % methods
    
    methods (Static)
        
        function hlp = create(indvarsArg, stateNum,...
                fs, eps, thetaStateIndex, theta, psys, zetas, guards, degree,...
                lambdas, res)
            
            import lp4.HybridLinearProgramVerificationWithGivenLambdaAndRe
            hlp = HybridLinearProgramVerificationWithGivenLambdaAndRe(indvarsArg, stateNum, size(guards, 2));
            hlp = hlp.init(fs, eps, thetaStateIndex, theta, psys, zetas, guards, degree, lambdas, res);
        end
        
    end % methods (Static)
    
end
