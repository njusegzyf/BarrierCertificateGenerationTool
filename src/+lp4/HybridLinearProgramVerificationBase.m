classdef HybridLinearProgramVerificationBase
    
    % TODO
    
    properties
        indvars % 问题1中的独立变量，例如x = [x1, x2, ..., xn]，类型为符号化变量构成的行向量
        
        stateNum
        thetaStateIndex
        
        guardNum
        
        eps
        
        phys
        
        fs % 在不同状态下的 f
        
        theta
        psys
        zetas
        guards
        
        decvars % 问题1中的决策变量，包括 不同状态下Phy的决策变量, Cαβ, Cγδ, Cu,v
        decvarsIndexes
        
        exprs
        
        linprogF
        
        isAttachRou
        rouVar
    end % properties
    
    methods
        
        function this = setLinprogF(this)
            % init f with all zeros
            this.linprogF = zeros(1, length(this.decvars));
            
            if this.isAttachRou
                this.linprogF(this.decvarsIndexes.rouIndex) = 1;
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
        
        function [this, solveRes, resNorms] = solve(this)
            if lp4.Lp4Config.IS_USE_CVX
                [this, solveRes, resNorms] = this.solveWithCvx();
            else 
                [this, solveRes, resNorms] = this.solveWithLinprog();
            end
        end
        
        function [this, solveRes, resNorms] = solveWithLinprog(this)
            
            [Aeq, beq, Aie, bie] = collectEqsAndIes(this);
            
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
            
            [Aeq, beq, Aie, bie] = collectEqsAndIes(this);
            
            if ~this.isAttachRou % if we do not need to minimize rou (not use linprogF)
                
                cvx_begin
                
                variable x(length(this.decvars));
                subject to
                Aie * x <= bie
                Aeq * x == beq
                
                cvx_end
                
            else % if we need to minimize rou (use linprogF)
                
                rouIndex = this.decvarsIndexes.rouIndex;
                
                cvx_begin
                
                variable x(length(this.decvars));
                minimize( x(rouIndex) )
                subject to
                Aie * x <= bie
                Aeq * x == beq
                
                cvx_end
                
            end
            
            solveRes = this.createCvxSolveRes(x, cvx_optval, cvx_status, cvx_cputime);
            
            import lp4.Lp4Config
            if Lp4Config.isDebug()
                solveRes.printSolution();
            end
            
            % FIXME
            resNorms = [];
%             if solveRes.hasSolution()
%                 resNorms = solveRes.computeAllExprsNorms();
%             else
%                 resNorms = [];
%             end
             
        end % function solveWithCvx
        
        function solveRes = createSolveRes(this, x, fval, flag, time)
            error("Sub classes should override createSolveRes method.");
            solveRes = 0;
        end
        
        function cvxSolveRes = createCvxSolveRes(linearProgram, x, cvxOptval, cvxStatus, cvxCpuTime)
            cvxSolveRes = lp4.HybridLinearProgramCvxVerificationSolveRes(linearProgram, x, cvxOptval, cvxStatus, cvxCpuTime);
        end
        
    end % methods
    
    methods (Static)
        
    end % methods (Static)
    
end
