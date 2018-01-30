classdef HybridLinearProgramCvxVerificationSolveRes < lp4util.CvxSolveRes & lp4.HybridLinearProgramVerificationSolveResBase
    
    methods
        
        function this = HybridLinearProgramCvxVerificationSolveRes(linearProgram, x, cvxOptval, cvxStatus, cvxCpuTime)
            % `linearProgram` should be of type 'lp4.HybridLinearProgramVerificationBase'
            errorIfWrongType(linearProgram, 'lp4.HybridLinearProgramVerificationBase', 'linearProgram');
            
            % Note: We must initialize the subclass object with calls to each superclass constructor
            % (call both constructors of CvxSolveRes and HybridLinearProgramVerificationSolveResBase),
            % so we define a constructor that supports the no argument case in HybridLinearProgramVerificationBase.
            % @see http://cn.mathworks.com/help/matlab/matlab_oop/subclassing-multiple-classes.html?searchHighlight=multiple%20inheritance&s_tid=doc_srchtitle
            
            this@lp4util.CvxSolveRes(linearProgram, x, cvxOptval, cvxStatus, cvxCpuTime);
            this@lp4.HybridLinearProgramVerificationSolveResBase();
            
        end % end ctor HybridLinearProgramCvxVerificationRes
        
    end % methods
end
