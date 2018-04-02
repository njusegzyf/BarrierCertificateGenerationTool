classdef LinearProgram4VerificationBase  < lp4.Lp4AndHlpVerificationBase
    
    properties
        f % 问题1中的f
        
        phy % 问题1中的多项式φ
        lambda % 问题1中的多项式λ
        
        c1Length % 记录约束 1-3 的决策变量 C 的长度
        c2Length
        c3Length
        
        rouIndex = -1
    end % properties
    
    methods
        
        function this = setThetaConstraint(this, theta)
            % for an empty constrain
            if isempty(theta)
                expr = Constraint.createEmptyConstraint();
                expr.num = 1;
                this.exprs = [expr];
                return;
            end
            
            expr = Constraint();
            expr.num = length(this.exprs) + 1;
            expr.name = 'theta';
            expr.type = 'eq';
            expr.A = [];
            expr.b = [];
            
            constraint1 = -this.phy;
            if lp4.Lp4Config.IS_ADD_EPS_IN_THETA_EXP 
                constraint1 = constraint1 + this.eps(2);
            end
                
            de = computeDegree(constraint1, this.indvars) + lp4.Lp4Config.VERIFICATION_C_DEGREE_INC;
            
            c_alpha_beta = sym('c_alpha_beta', [1, lp4.Lp4Config.DEFAULT_DEC_VAR_SIZE]); % pre-defined varibales, only a few of them are the actual variables
            [constraintDecvars, expression] = constraintExpression(de, theta, c_alpha_beta);
            this = this.addDecisionVars(constraintDecvars);
            this.c1Length = length(constraintDecvars);
            
            constraint1 = constraint1 + expression;
            expr.polyexpr = constraint1;
            
            this.exprs = [this.exprs expr];
        end
        
        % This is the main difference from `LinearProgram4`.
        function this = setPsyConstraint(this, psy)
            % for an empty constrain
            if isempty(psy)
                expr = Constraint.createEmptyConstraint();
                expr.num = 2;
                this.exprs(2) = expr;
                return;
            end
            
            expr = Constraint();
            expr.num = length(this.exprs) + 1;
            expr.name = 'psy';
            expr.type = 'eq';
            expr.A = [];
            expr.b = [];
            
            phy_d = 0;
            for k = 1 : 1 : length(this.indvars)
                phy_d = phy_d + diff(this.phy, this.indvars(k)) * this.f(k);
            end
            
            constraint2 = -phy_d + this.phy * this.lambda + this.eps(1);
            de = computeDegree(constraint2, this.indvars) + lp4.Lp4Config.VERIFICATION_C_DEGREE_INC;
            
            c_gama_delta = sym('c_gama_delta', [1, lp4.Lp4Config.DEFAULT_DEC_VAR_SIZE]);
            [constraintDecvars, expression] = constraintExpression(de, psy, c_gama_delta);
            this = this.addDecisionVars(constraintDecvars);
            this.c2Length = length(constraintDecvars);
            
            constraint2 = constraint2 + expression;
            
            expr.polyexpr = constraint2;
            
            this.exprs = [this.exprs expr];
        end
        
        function this = setZetaConstraint(this, zeta)
            % for an empty constrain
            if isempty(zeta)
                expr = Constraint.createEmptyConstraint();
                expr.num = 3;
                this.exprs(3) = expr;
                return;
            end
            
            expr = Constraint();
            expr.num = length(this.exprs) + 1;
            expr.name = 'zeta';
            expr.type = 'eq';
            expr.A = [];
            expr.b = [];
            
            constraint3 = this.phy + this.eps(2); % different from lp2
            de = computeDegree(constraint3, this.indvars) + lp4.Lp4Config.VERIFICATION_C_DEGREE_INC;
            
            c_u_v = sym('c_u_v', [1, lp4.Lp4Config.DEFAULT_DEC_VAR_SIZE]);
            [constraintDecvars, expression] = constraintExpression(de,zeta,c_u_v);
            this = this.addDecisionVars(constraintDecvars);
            this.c3Length = length(constraintDecvars);
            
            constraint3 = constraint3 + expression;
            
            expr.polyexpr = constraint3;
            
            this.exprs = [this.exprs expr];
        end
        
        function this = generateEqsForConstraint1To3(this)
            this = lp4.Lp4AndHlpVerificationBase.generateConstraintEqsParallelly(this);
        end
        
        function res = getRouIndex(this)
            res = this.rouIndex;
        end
        
        function [constraintDecvars, expression, indexVectors] = getThetaConstraintRightExpr(this, theta)            
            constraint1 = -this.phy;
            de = computeDegree(constraint1, this.indvars);
            
            c_alpha_beta = sym('c_alpha_beta', [1, lp4.Lp4Config.DEFAULT_DEC_VAR_SIZE]);
            [constraintDecvars, expression, indexVectors] = getConstraintExpressionAsVectors(de, theta, c_alpha_beta);
        end

        function [constraintDecvars, expression, indexVectors] = getPsyConstraintRightExpr(this, psy)
            phy_d = 0;
            for k = 1 : 1 : length(this.indvars)
                phy_d = phy_d + diff(this.phy, this.indvars(k)) * this.f(k);
            end
            
            constraint2 = -phy_d + this.phy * this.lambda + this.eps(1);
            de = computeDegree(constraint2, this.indvars) + lp4.Lp4Config.VERIFICATION_C_DEGREE_INC;
            
            c_gama_delta = sym('c_gama_delta', [1, lp4.Lp4Config.DEFAULT_DEC_VAR_SIZE]);
            [constraintDecvars, expression, indexVectors] = getConstraintExpressionAsVectors(de, psy, c_gama_delta);
        end
        
        function  [constraintDecvars, expression, indexVectors] = getZetaConstraintRightExpr(this, zeta)
            constraint3 = this.phy + this.eps(2);
            de = computeDegree(constraint3, this.indvars) + lp4.Lp4Config.VERIFICATION_C_DEGREE_INC;
            
            c_u_v = sym('c_u_v', [1, lp4.Lp4Config.DEFAULT_DEC_VAR_SIZE]);
            [constraintDecvars, expression, indexVectors] = getConstraintExpressionAsVectors(de, zeta, c_u_v);
        end

    end % methods
    
    methods (Static)
        
    end % methods (Static)
    
end
