classdef (Abstract) Lp4AndHlpVerificationBase
    %LP4ANDHLLPVERIFICATIONBASE
    
    properties
        indvars % ����1�еĶ�������������x = [x1, x2, ..., xn]������Ϊ���Ż��������ɵ�������
        eps
        theta
        decvars % �����еľ��߱���
        exprs % ���еĵ�ʽ/����ʽԼ��
        linprogF % ���Թ滮��Ŀ�꺯��
        isAttachRou = false; % �Ƿ��ڲ���ʽԼ���м��� rou������ rou ��ΪĿ�꺯��
        rouVar
    end
    
    methods
        
        function [this, startIndex, endIndex] = addDecisionVars(this, decvars)
            % addDecisionVars add decision variables.
            %
            % lp is a linear program.
            % decvars are decsion variables and should be symbolic.
            
            if isa(decvars, 'sym')
                % ���µľ��߱�����Ϊ����������ʽ��ӵ�lp��lp.decvars������
                % decision vars can only be a vector of a matrix of symbolic variables
                % reshape `decvars` to of dim `[1, size(decvars, 1) * size(decvars, 2)]`
                startIndex = size(this.decvars, 2) + 1;
                this.decvars = [this.decvars reshape(decvars, 1, size(decvars, 1) * size(decvars, 2))];
                endIndex = size(this.decvars, 2);
            else
                error('Wrong decvars.');
            end
        end
        
        function this = setLinprogF(this)
            % init f with all zeros
            this.linprogF = zeros(1, length(this.decvars));
            
            if this.isAttachRou
                this.linprogF(this.getRouIndex()) = 1;
            end
        end
        
        function [Aeq, beq, Aie, bie] = collectEqsAndIes(this)
            Aeq = [];
            beq = [];
            Aie = [];
            bie = [];
            
            for k = 1 : 1 : length(this.exprs)
                if this.exprs(k).isEmptyConstraint()
                    continue;
                end
                
                if strcmp(this.exprs(k).type, 'eq')
                    Aeq = [Aeq; this.exprs(k).A];
                    beq = [beq; this.exprs(k).b];
                else
                    Aie = [Aie; this.exprs(k).A];
                    bie = [bie; this.exprs(k).b];
                end
            end
        end
        
        function this = generateEqsForSafetyConstraints(this)
            this = lp4.Lp4AndHlpVerificationBase.generateConstraintEqsParallelly(this);
        end
        
        function [this, solveRes, resNorms] = solve(this)
            if lp4.Lp4Config.IS_USE_CVX
                [this, solveRes, resNorms] = this.solveWithCvx();
            else
                [this, solveRes, resNorms] = this.solveWithLinprog();
            end
        end
        
        function [this, solveRes, resNorms] = solveWithLinprog(this)
            
            [Aeq, beq, Aie, bie] = this.collectEqsAndIes();
            
            if ~lp4.Lp4Config.IS_SET_LINPROG_LOWERBOUND
                tic;
                [x, fval, flag, ~] = linprog(this.linprogF, Aie, bie, Aeq, beq);
                time = toc;
            else
                lb = repmat(lp4.Lp4Config.LINPROG_LOWERBOUND, 1, length(this.decvars));
                tic;
                [x, fval, flag, ~] = linprog(this.linprogF, Aie, bie, Aeq, beq, lb);
                time = toc;
            end
            
            if lp4.Lp4Config.IS_DROP_NEGATIVE_C
                x = this.dropNegativeC(x);
            end
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
        
        function [this, solveRes, resNorms] = solveWithCvx(this)
            
            [Aeq, beq, Aie, bie] = this.collectEqsAndIes();
            
            if ~this.isAttachRou % if we do not need to minimize rou (not use linprogF)
                
                cvx_begin
                
                variable x(length(this.decvars));
                subject to
                Aie * x <= bie
                Aeq * x == beq
                
                cvx_end
                
            else % if we need to minimize rou (use linprogF)
                
                rouIndex = this.getRouIndex();
                
                cvx_begin
                
                variable x(length(this.decvars));
                minimize( x(rouIndex) )
                subject to
                Aie * x <= bie
                Aeq * x == beq
                
                cvx_end
                
            end
            
            if lp4.Lp4Config.IS_DROP_NEGATIVE_C
                x = this.dropNegativeC(x);
            end
            solveRes = this.createCvxSolveRes(x, cvx_optval, cvx_status, cvx_cputime);
            
            import lp4.Lp4Config
            if Lp4Config.isDebug()
                solveRes.printSolution();
            end
            
            if solveRes.hasSolution()
                resNorms = solveRes.computeAllExprsNorms();
            else
                resNorms = [];
            end
            
        end % function solveWithCvx
        
        function x = dropNegativeC(this, x)
            xLen = length(x);
            for i = this.getCStart() : xLen 
                if x(i) < 0 
                    x(i) = 0;
                end
            end
        end
        
    end
    
    methods (Abstract)
        
        res = getRouIndex(this)
        
        solveRes = createSolveRes(this, x, fval, flag, time)
        
        res = getCStart(this)
        
    end
    
    methods (Static)
        
		function this = generateConstraintEqs(this)
            
            indVars = this.indvars;
            decVars = this.decvars;
            
            for i = 1 : length(this.exprs) 
			    expr = this.exprs;

                % skip non `eq` / empty constraint
                if ~strcmp(expr.type, 'eq') || expr.isEmptyConstraint() || expr.isEqGenerated()
                    continue;
                end
                
                [ expr.A, expr.b ] = eqgenerate(indVars, decVars, expr.polyexpr);
                disp(['constraint ', expr.name,' is processed: ',datestr(now,'yyyy-mm-dd HH:MM:SS')]);
                
                 this.exprs(i) = expr;
            end % for
        end % generateConstraintEqs
		
        function this = generateConstraintEqsParallelly(this)
            % ���л� eqgenerate ������ÿ��Լ��������ϵ��������ɵ�ʽԼ����
            
            % Note: Matlab ������ parfor ѭ���г������� this.expr(1) �ı���
            indVars = this.indvars;
            decVars = this.decvars;
            % Note: ������鲻�ɸ�ֵΪ`[]`�����򱨴� ������������ά�ȣ�
            outExprs = this.exprs;
            
            parfor i = 1 : length(this.exprs)
                expr = outExprs(i);
                % skip non `eq` / empty constraint
                if ~strcmp(expr.type, 'eq') || expr.isEmptyConstraint() || expr.isEqGenerated()
                    continue;
                end
                
                [ expr.A, expr.b ] = eqgenerate(indVars, decVars, expr.polyexpr);
                disp(['constraint ', expr.name,' is processed: ',datestr(now,'yyyy-mm-dd HH:MM:SS')]);
                
                outExprs(i) = expr;
            end % parfor
            
            this.exprs = outExprs;
        end % generateConstraintEqsParallelly
        
    end % methods (Static)
    
end
