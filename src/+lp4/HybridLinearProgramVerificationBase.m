classdef HybridLinearProgramVerificationBase < lp4.Lp4AndHlpVerificationBase
    
    properties
        stateNum
        thetaStateIndex
        guardNum
        
        phys % 在不同状态下的 phy
        fs % 在不同状态下的 f
        
        theta
        psys
        zetas
        guards
        
        decvarsIndexes
    end % properties
    
    methods

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

            leftDegree = computeDegree(constraint1, this.indvars);
            de = lp4.Lp4Config.getVerificationCDegree(leftDegree);
            
            if (lp4.Lp4Config.IS_USE_NEW_CONSTRAINT_GENERATION_FUNC) 
                [constraintDecvars, expression] = lp4util.generateConstraintExpression(de, theta, 'c_alpha_beta');
            else 
                % c_alpha_beta = sym('c_alpha_beta', [1, lp4.Lp4Config.DEFAULT_DEC_VAR_SIZE]); % pre-defined varibales, only a few of them are the actual variables
                c_alpha_beta = sym('c_alpha_beta', [1, lp4.Lp4Config.getDecVarArraySize(length(theta) * 2, de)]);
                [constraintDecvars, expression] = constraintExpression(de, theta, c_alpha_beta);
            end

            [this, this.decvarsIndexes.cThetaStart, this.decvarsIndexes.cThetaEnd] = this.addDecisionVars(constraintDecvars);
            
            constraint1 = constraint1 + expression;
            
            constraint1 = expand(constraint1);
            expr.polyexpr = constraint1;
            
            % Note: This is the first expression, and use `this.exprs(expr.num) = expr;` here
            % will cause an error in Matlab 2014b but is Ok in Matlab 2017b.
            this.exprs = [this.exprs, expr];
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
                
                leftDegree = computeDegree(constraint3, this.indvars);
                de = lp4.Lp4Config.getVerificationCDegree(leftDegree);
                
                name = strcat('c_u_v', num2str(zetaIndex));
                
                if (lp4.Lp4Config.IS_USE_NEW_CONSTRAINT_GENERATION_FUNC)
                    [constraintDecvars, expression] = lp4util.generateConstraintExpression(de, zeta.exprs, name);
                else
                    c_u_v = sym(name,[1, lp4.Lp4Config.getDecVarArraySize(length(zeta.exprs) * 2, de)]);
                    [constraintDecvars, expression] = constraintExpression(de, zeta.exprs, c_u_v);
                end
                
                [this, this.decvarsIndexes.cZetaStarts(zetaIndex), this.decvarsIndexes.cZetaEnds(zetaIndex)] = this.addDecisionVars(constraintDecvars);
                
                constraint3 = constraint3 + expression;
                
                constraint3 = expand(constraint3);
                expr.polyexpr = constraint3;
                
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
                
                leftDegree = computeDegree(constraintRe, this.indvars);
                de = lp4.Lp4Config.getVerificationCDegree(leftDegree);
                
                % if re is a constant, ignore it
                %                 if de == 0
                %                     continue;
                %                 end
                
                name = strcat('c_re', num2str(i), '_');
                
                if (lp4.Lp4Config.IS_USE_NEW_CONSTRAINT_GENERATION_FUNC)
                    [constraintDecvars, expression] = lp4util.generateConstraintExpression(de, guard.exprs, name);
                else
                    c_re = sym(name, [1, lp4.Lp4Config.getDecVarArraySize(length(guard.exprs) * 2, de)]);
                    [constraintDecvars, expression] = constraintExpression(de, guard.exprs, c_re);
                end

                [this, this.decvarsIndexes.cReStarts(i), this.decvarsIndexes.cReEnds(i)] = this.addDecisionVars(constraintDecvars);
                
                constraintRe = constraintRe + expression;
                
                constraintRe = expand(constraintRe);
                expr.polyexpr = constraintRe;
                
                this.exprs(exprNum) = expr;
            end
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
            
            for k = 1 : cLength
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
        
        % override
        function  res = getRouIndex(this)
            res = this.decvarsIndexes.rouIndex;
        end
        
        % override
        function res = getCStart(this)
            res = this.decvarsIndexes.cThetaStart(1);
        end
        
        function res = nextExprNumIndex(this)
            res = length(this.exprs) + 1;
        end
        
        function cvxSolveRes = createCvxSolveRes(linearProgram, x, cvxOptval, cvxStatus, cvxCpuTime)
            cvxSolveRes = lp4.HybridLinearProgramCvxVerificationSolveRes(linearProgram, x, cvxOptval, cvxStatus, cvxCpuTime);
        end
        
    end % methods
    
    methods (Static)
        
    end % methods (Static)
    
end
