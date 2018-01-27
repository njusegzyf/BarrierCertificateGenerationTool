classdef DecVarsAndExpr

    properties
        symbolicVars % of matrix
        expression
    end
    
    methods
        function this = DecVarsAndExpr(symbolicVars, expression)
            this.symbolicVars = symbolicVars;
            this.expression = expression;
        end

    end
end

