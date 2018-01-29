classdef LinearProgram4Verification2
    %LinearProgram4Verification2 A linear program used to verify the solution of LinearProgram4
    % with given lambda.
    
    properties
        indvars % 问题1中的独立变量，例如x = [x1, x2, ..., xn]，类型为符号化变量构成的行向量
        degree % 问题1中的带求多项式函数φ的次数，类型为正整数
        
        phy % 问题1中的多项式φ
        phySize
        lambda % 问题1中的多项式λ
        
        eps % 问题1中的?1和?1构成的向量
        f % 问题1中的f
        decvars % 待验证解中的决策变量，包括 phy, Cα,β, Cγ,δ, Cu,v
        exprs
        
        phyPolynomial
        
        c1Length
        c2Length
        c3Length
        
        pPartitions = []
    end % properties
    
    methods
        
        function this = LinearProgram4Verification2(indvarsArg)
            
            % vars can only be a vector of a matrix of symbolic variables
            if ~isa(indvarsArg, 'sym')
                error('vars is not a vector of a matrix of symbolic variables');
            end
            
            this.indvars = reshape(indvarsArg, 1, size(indvarsArg, 1) * size(indvarsArg, 2));
            
            this.degree = 0;
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
        
        function lp = set.lambda(lp, lambda)
            lp.lambda = lambda;
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
        
        function this = setDegreeAndInit(this, degree)
            this.degree = degree;
            
            % set phy
            this.phySize = monomialNumber(length(this.indvars), degree);
            p = sym('p', [1, this.phySize]);
            this = this.addDecisionVars(p);
            this.phy = p * monomials(this.indvars, 0 : degree);
            
            import lp4util.SymbolicPolynomial
            this.phyPolynomial = SymbolicPolynomial(this.indvars, degree, p, this.phy);
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
            import lp4.Lp4Config
            de = computeDegree(constraint1, this.indvars) + Lp4Config.VERIFICATION_C_DEGREE_INC;
            
            c_alpha_beta = sym('c_alpha_beta', [1,10000]); % pre-defined varibales, only a few of them are the actual variables
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
            import lp4.Lp4Config
            de = computeDegree(constraint2, this.indvars) + Lp4Config.VERIFICATION_C_DEGREE_INC;
            
            c_gama_delta = sym('c_gama_delta',[1,10000]);
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
            import lp4.Lp4Config;
            de = computeDegree(constraint3, this.indvars) + Lp4Config.VERIFICATION_C_DEGREE_INC;
            
            c_u_v = sym('c_u_v',[1,10000]);
            [constraintDecvars, expression] = constraintExpression(de,zeta,c_u_v);
            this = this.addDecisionVars(constraintDecvars);
            this.c3Length = length(constraintDecvars);
            
            constraint3 = constraint3 + expression;
            
            expr.polyexpr = constraint3;
            
            this.exprs = [this.exprs expr];
        end
        
        function this = generateEqsForConstraint1To3(this)
            for k = 1 : 1 : 3
                if this.exprs(k).isEmptyConstraint()
                    continue;
                end
                
                [ this.exprs(k).A, this.exprs(k).b ] = eqgenerate( this.indvars, this.decvars, this.exprs(k).polyexpr);
                disp(['constraint ',this.exprs(k).name,' is processed: ',datestr(now,'yyyy-mm-dd HH:MM:SS')]);
            end
        end
        
        function this = setPhyConstraint(this)
            if isempty(this.pPartitions)
                return;
            end
            
            expr = Constraint();
            expr.num = length(this.exprs) + 1;
            expr.name = 'phy';
            expr.type = 'ie';
            expr.A = [];
            expr.b = [];
            
            for i = 1 : 1 : length(this.phyPolynomial.coefficientVars)
                p = this.phyPolynomial.coefficientVars(i);
                pPartition = this.pPartitions(i);
                
                expr = pPartition.createConstraintsAndAddToExpr(p, this.decvars, expr);
            end
            
            this.exprs = [this.exprs expr];
        end
        
        function this = setDevVarsConstraint(this)
            decexpr = Constraint();
            decexpr.num = length(this.exprs) + 1;
            decexpr.name = 'decvarconstraints';
            decexpr.type = 'ie';
            decexpr.polyexpr = [];
            
            cStart = this.phySize;
            cLength = length(this.decvars) - cStart;
            decexpr.A = zeros(cLength, length(this.decvars));
            for k = 1 : 1 : cLength
                decexpr.A(k, cStart + k) = -1;
            end
            bc = zeros(cLength, 1);
            decexpr.b = bc;
            
            this.exprs = [this.exprs decexpr];
        end % function setDevVarsConstraint
        
        function [this, solveRes, resNorms] = solve(this)
            Aeq = [];
            beq = [];
            Aie = [];
            bie = [];
            
            for k = 1 : 1 : length(this.exprs)
                if this.exprs(k).isEmptyConstraint()
                    continue;
                end
                
                if this.exprs(k).type == 'eq'
                    Aeq = [Aeq; this.exprs(k).A];
                    beq = [beq; this.exprs(k).b];
                else
                    Aie = [Aie; this.exprs(k).A];
                    bie = [bie; this.exprs(k).b];
                end
            end
            
            linprogF = zeros(1, length(this.decvars));
            tic;
            [x, fval, flag, ~] = linprog(linprogF, Aie, bie, Aeq, beq);
            time = toc;
            
            import lp4.LinearProgram4Verification2SolveResult
            solveRes = LinearProgram4Verification2SolveResult(this, x, fval, flag, time);
            
            import lp4.Lp4Config
            if Lp4Config.isDebug()
                solveRes.printSolution();
            end
            
            if solveRes.hasSolution()
                resNorms = solveRes.verifyNorms();
            else
                resNorms = [];
            end
        end % function solve
        
        function res = getPhyCoefficientStart(this)
            res = 1;
        end
        
        function res = getPhyCoefficientLength(this)
            res = length(this.phyPolynomial.coefficientVars);
        end
        
        function res = getC1Start(this)
            res = this.getPhyCoefficientStart() + this.getPhyCoefficientLength();
        end
        
        function res = getC2Start(this)
            res = this.getC1Start() + this.c1Length;
        end
        
        function res = getC3Start(this)
            res = this.getC2Start() + this.c2Length;
        end
        
    end % methods
    
end % classdef
