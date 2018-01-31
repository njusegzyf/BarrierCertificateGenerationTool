classdef LinearProgram4Verification2 < lp4.LinearProgram4VerificationBase
    %LinearProgram4Verification2 A linear program used to verify the solution of LinearProgram4 with given lambda.
    
    properties
        degree % 问题1中的带求多项式函数φ的次数，类型为正整数
        phySize
        phyPolynomial
        pPartitions = []
    end % properties
    
    methods
        
        function this = LinearProgram4Verification2(indvarsArg)
            
            % vars can only be a vector of a matrix of symbolic variables
            if ~isa(indvarsArg, 'sym')
                error('vars is not a vector of a matrix of symbolic variables');
            end
            
            this.indvars = lp4util.reshapeToVector(indvarsArg);
            
            this.degree = 0;
            this.eps = [];
            this.f = [];
            this.decvars = [];
        end
        
        function this = init(this, f, eps, theta, psy, zeta, degree, lambda, phyRangeInVerify)
            this.f = f;
            this.eps = eps;
            
            % set the degree of phy
            this = this.setDegreeAndInit(degree + lp4.Lp4Config.VERIFICATION_PHY_DEGREE_INC);
            
            % set lambda expression
            this.lambda = lambda;
            
            this = this.setThetaConstraint(theta);
            this = this.setPsyConstraint(psy);
            this = this.setZetaConstraint(zeta);
            this = this.generateEqsForConstraint1To3();
            
            if isa(phyRangeInVerify, 'lp4util.Partition')
                this.pPartitions = repmat(phyRangeInVerify, 1024, 1);
                this = this.setPhyConstraint();
            end
            
            this = this.setDevVarsConstraint();
            
            this = this.setLinprogF();
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
            
            if this.isAttachRou
                this.rouVar = sym('rou');
                [this, this.rouIndex, ~] = this.addDecisionVars(this.rouVar);
            end
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
            import lp4.Lp4Config
            de = computeDegree(constraint2, this.indvars) + Lp4Config.VERIFICATION_C_DEGREE_INC;
            
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
            import lp4.Lp4Config;
            de = computeDegree(constraint3, this.indvars) + Lp4Config.VERIFICATION_C_DEGREE_INC;
            
            c_u_v = sym('c_u_v', [1, lp4.Lp4Config.DEFAULT_DEC_VAR_SIZE]);
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
            if this.isAttachRou
                rouIndex = this.getRouIndex();
                for k = 1 : cLength
                    % - rou
                    decexpr.A(k, rouIndex) = -1;
                end
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
            
            tic;
            [x, fval, flag, ~] = linprog(this.linprogF, Aie, bie, Aeq, beq);
            time = toc;
            
            solveRes = this.createSolveRes(x, fval, flag, time);
            
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
        
        function res = getPhyCoefficientStart(this)
            res = 1;
        end
        
        function res = getPhyCoefficientLength(this)
            res = length(this.phyPolynomial.coefficientVars);
        end
        
        function res = getC1Start(this)
            res = this.getPhyCoefficientStart() + this.getPhyCoefficientLength();
            if this.isAttachRou
                res = res + 1;
            end
        end
        
        function res = getC2Start(this)
            res = this.getC1Start() + this.c1Length;
        end
        
        function res = getC3Start(this)
            res = this.getC2Start() + this.c2Length;
        end
        
        function solveRes = createSolveRes(this, x, fval, flag, time)
            solveRes = lp4.LinearProgram4Verification2SolveResult(this, x, fval, flag, time);
        end
        
    end % methods
    
    methods (Static)
        
        function lp = create(indvars, f, eps, theta, psy, zeta, degree, lambda, phyRangeInVerify)
            lp = lp4.LinearProgram4Verification2(indvars);
            lp = lp.init(f, eps, theta, psy, zeta, degree, lambda, phyRangeInVerify);
        end
        
        function lp = createWithRou(indvars, f, eps, theta, psy, zeta, degree, lambda, phyRangeInVerify)
            lp = lp4.LinearProgram4Verification2(indvars);
            lp.isAttachRou = true;
            lp = lp.init(f, eps, theta, psy, zeta, degree, lambda, phyRangeInVerify);
        end
        
    end % methods (Static)
    
end % classdef
