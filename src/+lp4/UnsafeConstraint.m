classdef UnsafeConstraint
    %UNSAFECONSTRAINT 
    
    properties
        stateIndex
        exprs
    end
    
    methods
        function this = UnsafeConstraint(stateIndex, exprs)
            this.stateIndex = stateIndex;
            this.exprs = exprs;
        end
    end
end

