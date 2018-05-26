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
