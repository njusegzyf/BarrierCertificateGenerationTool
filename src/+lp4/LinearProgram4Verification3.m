classdef LinearProgram4Verification3
    %LinearProgram4Verification3 A linear program used to verify the solution of LinearProgram4
    % with given phy.
    
    properties
        indvars % 问题1中的独立变量，例如x = [x1, x2, ..., xn]，类型为符号化变量构成的行向量
        lambdaDegree % 问题1中的带求多项式函数λ的次数，类型为正整数
        
        phy % 问题1中的多项式φ
        lambda % 问题1中的多项式λ
        lambdaSize
        
        eps % 问题1中的?1和?1构成的向量
        f % 问题1中的f
        decvars % 待验证解中的决策变量，包括 λ, Cα,β, Cγ,δ, Cu,v
        exprs
        
        lambdaPolynomial
        
        solution
        
        c1Length
        c2Length
        c3Length
    end % properties
    
    methods
        
        function this = LinearProgram4Verification3(indvarsArg)
            
            % vars can only be a vector of a matrix of symbolic variables
            if ~isa(indvarsArg, 'sym')
                error('vars is not a vector of a matrix of symbolic variables');
            end
            
            this.indvars = reshape(indvarsArg, 1, size(indvarsArg, 1) * size(indvarsArg, 2));
            
            this.lambdaDegree = 0;
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
        
        function lp = set.phy(lp, phy)
            lp.phy = phy;
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
        
        function this = setDegreeAndInit(this, lambdaDegree)
            this.lambdaDegree = lambdaDegree;
            
            % set lambda
            this.lambdaSize = monomialNumber(length(this.indvars), lambdaDegree);
            pLambda = sym('pLambda', [1, this.lambdaSize]);
            this = this.addDecisionVars(pLambda);
            this.lambda = pLambda * monomials(this.indvars, 0 : lambdaDegree);
            
            import lp4util.SymbolicPolynomial
            this.lambdaPolynomial = SymbolicPolynomial(this.indvars, lambdaDegree, pLambda, this.lambda);
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
            
            c_alpha_beta = sym('c_alpha_beta', [1,10000]); % pre-defined varibales, only a few of them are the actual variables
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
            
            c_gama_delta = sym('c_gama_delta',[1,10000]);
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
            
            c_u_v = sym('c_u_v',[1,10000]);
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
            
            linprogF = zeros(1, length(this.decvars));
            tic;
            [x, fval, flag, ~] = linprog(linprogF, Aie, bie, Aeq, beq);
            time = toc;
            
            import lp4.LinearProgram4Verification3SolveResult
            solveRes = LinearProgram4Verification3SolveResult(this, x, fval, flag, time);
            
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
        
        function res = getLambdaCoefficientStart(this)
            res = 1;
        end
        
        function res = getLambdaCoefficientLength(this)
            res = length(this.lambdaPolynomial.coefficientVars);
        end
        
        function res = getC1Start(this)
            res = this.getLambdaCoefficientStart() + this.getLambdaCoefficientLength();
        end
        
        function res = getC2Start(this)
            res = this.getC1Start() + this.c1Length;
        end
        
        function res = getC3Start(this)
            res = this.getC2Start() + this.c2Length;
        end
        
    end % methods
    
end % classdef
