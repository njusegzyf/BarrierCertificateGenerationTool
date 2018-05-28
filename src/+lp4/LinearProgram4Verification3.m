classdef LinearProgram4Verification3 < lp4.LinearProgram4VerificationBase
    %LinearProgram4Verification3 A linear program used to verify the solution of LinearProgram4 with given phy.
    
    properties
        lambdaDegree % ����1�еĴ������ʽ�����˵Ĵ���������Ϊ������
        lambdaSize
        lambdaPolynomial
    end % properties
    
    methods
        
        function this = LinearProgram4Verification3(indvarsArg)
            this@lp4.LinearProgram4VerificationBase(indvarsArg);
            
            this.lambdaDegree = 0;
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
