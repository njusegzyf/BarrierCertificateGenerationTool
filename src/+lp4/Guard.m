classdef Guard
    
    properties
        fromStateIndex
        destStateIndex
        exprs
        resetVars
        resetVarValues
    end
    
    methods
        function this = Guard(fromStateIndex, destStateIndex, exprs, resetVars, resetVarValues)
            this.fromStateIndex = fromStateIndex;
            this.destStateIndex = destStateIndex;
            this.exprs = exprs;
            this.resetVars = resetVars;
            this.resetVarValues = resetVarValues;
        end
    end
end

