classdef HybridLinearProgramVerificationBase < lp4.Lp4AndHlpVerificationBase
    
    properties
        stateNum
        thetaStateIndex
        guardNum
        
        phys % �ڲ�ͬ״̬�µ� phy
        fs % �ڲ�ͬ״̬�µ� f
        
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
            import lp4.Lp4Config
            de = computeDegree(constraint1, this.indvars) + Lp4Config.C_DEGREE_INC;
            
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
                import lp4.Lp4Config
                de = computeDegree(constraint3, this.indvars) + Lp4Config.C_DEGREE_INC;
                
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
