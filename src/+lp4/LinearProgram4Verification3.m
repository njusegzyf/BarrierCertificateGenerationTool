classdef LinearProgram4Verification3 < lp4.LinearProgram4VerificationBase
    %LinearProgram4Verification3 A linear program used to verify the solution of LinearProgram4 with given phy.
    
    properties
        lambdaDegree % 问题1中的带求多项式函数λ的次数，类型为正整数
        lambdaSize
        lambdaPolynomial
    end % properties
    
    methods
        
        function this = LinearProgram4Verification3(indvarsArg)
            
            % vars can only be a vector of a matrix of symbolic variables
            if ~isa(indvarsArg, 'sym')
                error('vars is not a vector of a matrix of symbolic variables');
            end
            
            this.indvars = lp4util.reshapeToVector(indvarsArg);
            
            this.lambdaDegree = 0;
            this.eps = [];
            this.f = [];
            this.decvars = [];
        end
        
        function this = init(this, f, eps, theta, psy, zeta, pLambdaDegree, phy)
            this.f = f;
            this.eps = eps;
            
            % set the degree of lambda
            this = this.setDegreeAndInit(pLambdaDegree + lp4.Lp4Config.VERIFICATION_LAMBDA_DEGREE_INC);
            
            % set phy expression
            this.phy = phy;
            
            this = this.setThetaConstraint(theta);
            this = this.setPsyConstraint(psy);
            this = this.setZetaConstraint(zeta);
            this = this.generateEqsForConstraint1To3();
            
            this = this.setDevVarsConstraint();
            
            this = this.setLinprogF();
        end
        
        function this = setDegreeAndInit(this, lambdaDegree)
            this.lambdaDegree = lambdaDegree;
            
            % set lambda
            this.lambdaSize = monomialNumber(length(this.indvars), lambdaDegree);
            pLambda = sym('pLambda', [1, this.lambdaSize]);
            this = this.addDecisionVars(pLambda);
            this.lambda = pLambda * monomials(this.indvars, 0 : lambdaDegree);
            
            this.lambdaPolynomial = lp4util.SymbolicPolynomial(this.indvars, lambdaDegree, pLambda, this.lambda);
            
            if this.isAttachRou
                this.rouVar = sym('rou');
                [this, this.rouIndex, ~] = this.addDecisionVars(this.rouVar);
            end
        end
        
        function this = setThetaConstraint(this, thetaVars)
            expr = Constraint();
            expr.num = length(this.exprs) + 1;
            expr.name = 'theta';
            expr.type = 'eq';
            expr.A = [];
            expr.b = [];
            
            constraint1 = -this.phy;
            de = computeDegree(constraint1, this.indvars);
            
            c_alpha_beta = sym('c_alpha_beta', [1, lp4.Lp4Config.DEFAULT_DEC_VAR_SIZE]); % pre-defined varibales, only a few of them are the actual variables
            [constraintDecvars, expression] = constraintExpression(de, thetaVars, c_alpha_beta);
            this = this.addDecisionVars(constraintDecvars);
            this.c1Length = length(constraintDecvars);
            
            constraint1 = constraint1 + expression;
            expr.polyexpr = constraint1;
            
            this.exprs = [this.exprs expr];
        end
        
        % This is the main difference from `LinearProgram4`.
        function this = setPsyConstraint(this, psyVars)
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
            de = computeDegree(constraint2, this.indvars);
            
            c_gama_delta = sym('c_gama_delta', [1, lp4.Lp4Config.DEFAULT_DEC_VAR_SIZE]);
            [constraintDecvars, expression] = constraintExpression(de, psyVars, c_gama_delta);
            this = this.addDecisionVars(constraintDecvars);
            this.c2Length = length(constraintDecvars);
            
            constraint2 = constraint2 + expression;
            
            expr.polyexpr = constraint2;
            
            this.exprs = [this.exprs expr];
        end
        
        function this = setZetaConstraint(this, zetaVars)
            expr = Constraint();
            expr.num = length(this.exprs) + 1;
            expr.name = 'zeta';
            expr.type = 'eq';
            expr.A = [];
            expr.b = [];
            
            constraint3 = this.phy + this.eps(2); % different from lp2
            de = computeDegree(constraint3, this.indvars);
            
            c_u_v = sym('c_u_v', [1, lp4.Lp4Config.DEFAULT_DEC_VAR_SIZE]);
            [constraintDecvars, expression] = constraintExpression(de,zetaVars,c_u_v);
            this = this.addDecisionVars(constraintDecvars);
            this.c3Length = length(constraintDecvars);
            
            constraint3 = constraint3 + expression;
            
            expr.polyexpr = constraint3;
            
            this.exprs = [this.exprs expr];
        end
        
        function this = generateEqsForConstraint1To3(this)
            for k = 1 : 1 : 3
                [ this.exprs(k).A, this.exprs(k).b ] = eqgenerate( this.indvars, this.decvars, this.exprs(k).polyexpr);
                disp(['constraint ',this.exprs(k).name,' is processed: ',datestr(now,'yyyy-mm-dd HH:MM:SS')]);
            end
        end
        
        function this = setDevVarsConstraint(this)
            decexpr = Constraint();
            decexpr.num = length(this.exprs) + 1;
            decexpr.name = 'decvarconstraints';
            decexpr.type = 'ie';
            decexpr.polyexpr = [];
            
            cStart = this.lambdaSize;
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
            
            solveRes = lp4.LinearProgram4Verification3SolveResult(this, x, fval, flag, time);
            
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
        
        function res = getLambdaCoefficientStart(this)
            res = 1;
        end
        
        function res = getLambdaCoefficientLength(this)
            res = length(this.lambdaPolynomial.coefficientVars);
        end
        
        function res = getC1Start(this)
            res = this.getLambdaCoefficientStart() + this.getLambdaCoefficientLength();
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
            solveRes = lp4.LinearProgram4Verification3SolveResult(this, x, fval, flag, time);
        end
        
    end % methods
    
    methods (Static)
        
        function lp = create(indvars, f, eps, theta, psy, zeta, pLambdaDegree, phy)
            lp = lp4.LinearProgram4Verification3(indvars);
            lp = lp.init(f, eps, theta, psy, zeta, pLambdaDegree, phy);
        end
        
        function lp = createWithRou(indvars, f, eps, theta, psy, zeta, pLambdaDegree, phy)
            lp = lp4.LinearProgram4Verification3(indvars);
            lp.isAttachRou = true;
            lp = lp.init(f, eps, theta, psy, zeta, pLambdaDegree, phy);
        end
        
    end % methods (Static)
    
end % classdef
