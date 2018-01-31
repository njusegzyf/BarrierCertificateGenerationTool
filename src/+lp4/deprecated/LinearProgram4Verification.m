classdef LinearProgram4Verification
    %LINEARPROGRAM A linear program used to verify the solution of LinearProgram4 with given phy and lambda.
    
    properties
        indvars % ����1�еĶ�������������x = [x1, x2, ..., xn]������Ϊ���Ż��������ɵ�������
        degree % ����1�еĴ������ʽ�����յĴ���������Ϊ������
        
        phy % ����1�еĶ���ʽ��
        lambda % ����1�еĶ���ʽ��
        
        eps % ����1�е�?1��?1���ɵ�����
        f % ����1�е�f
        decvars % ����֤���еľ��߱��������� C��,��, C��,��, Cu,v
        exprs
        
        solution
    end % properties
    
    methods
        
        function this = LinearProgram4Verification(indvarsArg)
            
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
        
        function lp = set.degree(lp, degree)
            lp.degree = degree;
        end
            
        function lp = set.phy(lp, phy)
            lp.phy = phy;
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
                % ���µľ��߱�����Ϊ����������ʽ��ӵ�lp��lp.decvars������
                % decision vars can only be a vector of a matrix of symbolic variables
                % reshape `decvars` to of dim `[1, size(decvars, 1) * size(decvars, 2)]`
                this.decvars = [this.decvars reshape(decvars, 1, size(decvars, 1) * size(decvars, 2))];
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
            
            c_alpha_beta = sym('c_alpha_beta', [1,10000]); % pre-defined varibales, only a few of them are the actual variables
            [constraintDecvars, expression] = constraintExpression(de, thetaVars, c_alpha_beta);
            this = this.addDecisionVars(constraintDecvars);
            
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

			cStart = 0;
            cLength = length(this.decvars) - cStart;
            decexpr.A = zeros(cLength, length(this.decvars));
            for k = 1 : 1 : cLength
                decexpr.A(k, cStart + k) = -1;
            end
            bc = zeros(cLength, 1);
            decexpr.b = bc;
			
			this.exprs = [this.exprs decexpr];
        end % function setDevVarsConstraint
        
        function [this, solveRes] = solve(this)
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
            
            % for debug, print items in `Aie` and `bie` that is not of type double
            %             Aieflat = Aie(:);
            %             for i = 1 : 1 : length(Aieflat)
            %                 if ~isa(Aieflat(i), 'double')
            %                     disp(Aieflat(i));
            %                 end
            %             end
            %
            %             bieflat = bie(:);
            %             for i = 1 : 1 : length(bieflat)
            %                 if ~isa(bieflat(i), 'double')
            %                     disp(bieflat(i));
            %                 end
            %             end
            
            linprogF = zeros(1, length(this.decvars));
            tic;
            [x, fval, flag, ~] = linprog(linprogF, Aie, bie, Aeq, beq);
            time = toc;
            
            import lp4.LinearProgram4SolveResult
            solveRes = LinearProgram4SolveResult(this, x, fval, flag, time);
            
            import lp4.Lp4Config
            if Lp4Config.isDebug()
                solveRes.printSolution();
            end
            
        end % function solve
        
    end % methods
    
end % classdef
