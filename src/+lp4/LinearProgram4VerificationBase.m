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
        
        function this = LinearProgram4VerificationBase(indvarsArg)
                        
            % vars can only be a vector of a matrix of symbolic variables            
            if ~isa(indvarsArg, 'sym')
                error('vars is not a vector of a matrix of symbolic variables');
            end
            
            this.indvars = lp4util.reshapeToVector(indvarsArg);
            
            this.eps = [];
            this.f = [];
            this.decvars = [];
        end
        
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
                
            leftDegree = computeDegree(constraint1, this.indvars);
            de = lp4.Lp4Config.getVerificationCDegree(leftDegree);
            
            if (lp4.Lp4Config.IS_USE_NEW_CONSTRAINT_GENERATION_FUNC) 
                [constraintDecvars, expression] = lp4util.generateConstraintExpression(de, theta, 'c_alpha_beta');
            else 
                % c_alpha_beta = sym('c_alpha_beta', [1, lp4.Lp4Config.DEFAULT_DEC_VAR_SIZE]); % pre-defined varibales, only a few of them are the actual variables
                c_alpha_beta = sym('c_alpha_beta', [1, lp4.Lp4Config.getDecVarArraySize(length(theta) * 2, de)]);
                [constraintDecvars, expression] = constraintExpression(de, theta, c_alpha_beta);
            end
            
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
            
            leftDegree = computeDegree(constraint2, this.indvars);
            de = lp4.Lp4Config.getVerificationCDegree(leftDegree);
            
            if (lp4.Lp4Config.IS_USE_NEW_CONSTRAINT_GENERATION_FUNC) 
                [constraintDecvars, expression] = lp4util.generateConstraintExpression(de, psy, 'c_gama_delta');
            else 
                c_gama_delta = sym('c_gama_delta', [1, lp4.Lp4Config.getDecVarArraySize(length(psy) * 2, de)]);
                [constraintDecvars, expression] = constraintExpression(de, psy, c_gama_delta);
            end
            
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
            
            leftDegree = computeDegree(constraint3, this.indvars);
            de = lp4.Lp4Config.getVerificationCDegree(leftDegree);
            
            if (lp4.Lp4Config.IS_USE_NEW_CONSTRAINT_GENERATION_FUNC) 
                [constraintDecvars, expression] = lp4util.generateConstraintExpression(de, zeta, 'c_u_v');
            else 
                c_u_v = sym('c_u_v', [1, lp4.Lp4Config.getDecVarArraySize(length(zeta) * 2, de)]);
                [constraintDecvars, expression] = constraintExpression(de,zeta,c_u_v);
            end

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
        
        function res = getCStart(this)
            res = this.getC1Start();
        end
        
%         function [constraintDecvars, expression, indexVectors] = getThetaConstraintExpr(this, theta)            
%             constraint1 = -this.phy;
%             if lp4.Lp4Config.IS_ADD_EPS_IN_THETA_EXP 
%                 constraint1 = constraint1 + this.eps(2);
%             end
%             de = computeDegree(constraint1, this.indvars);
%             
%             c_alpha_beta = sym('c_alpha_beta', [1, lp4.Lp4Config.DEFAULT_DEC_VAR_SIZE]);
%             [constraintDecvars, expression, indexVectors] = getConstraintExpressionAsVectors(de, theta, c_alpha_beta);
%         end
% 
%         function [constraintDecvars, expression, indexVectors] = getPsyConstraintExpr(this, psy)
%             phy_d = 0;
%             for k = 1 : 1 : length(this.indvars)
%                 phy_d = phy_d + diff(this.phy, this.indvars(k)) * this.f(k);
%             end
%             
%             constraint2 = -phy_d + this.phy * this.lambda + this.eps(1);
%             de = computeDegree(constraint2, this.indvars) + lp4.Lp4Config.VERIFICATION_C_DEGREE_INC;
%             
%             c_gama_delta = sym('c_gama_delta', [1, lp4.Lp4Config.DEFAULT_DEC_VAR_SIZE]);
%             [constraintDecvars, expression, indexVectors] = getConstraintExpressionAsVectors(de, psy, c_gama_delta);
%         end
%         
%         function  [constraintDecvars, expression, indexVectors] = getZetaConstraintExpr(this, zeta)
%             constraint3 = this.phy + this.eps(2);
%             de = computeDegree(constraint3, this.indvars) + lp4.Lp4Config.VERIFICATION_C_DEGREE_INC;
%             
%             c_u_v = sym('c_u_v', [1, lp4.Lp4Config.DEFAULT_DEC_VAR_SIZE]);
%             [constraintDecvars, expression, indexVectors] = getConstraintExpressionAsVectors(de, zeta, c_u_v);
%         end
        
        % @override
        function cvxSolveRes = createCvxSolveRes(linearProgram, x, cvxOptval, cvxStatus, cvxCpuTime)
            exitflag = lp4util.CvxSolveRes.convertCvxStatusToExitflag(cvxStatus);
            cvxSolveRes = linearProgram.createSolveRes(x, cvxOptval, exitflag, cvxCpuTime);
            cvxSolveRes.output = cvxStatus;
        end

    end % methods
    
    methods (Static)

    end % methods (Static)
    
end
