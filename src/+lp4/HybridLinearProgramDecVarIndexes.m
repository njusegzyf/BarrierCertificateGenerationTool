classdef HybridLinearProgramDecVarIndexes
    %HYBRIDLINEARPROGRAMDECVARINDEXES
    
    properties
        phyStarts
        phyEnds
        
        rouIndex
        
        lambdaStarts
        lambdaEnds
        
        reStarts
        reEnds
        
        wLambdaStarts
        wLambdaEnds
        
        wReStarts
        wReEnds
        
        cThetaStart
        cThetaEnd
        
        cPsyStarts
        cPsyEnds
        
        cGuardStarts
        cGuardEnds
        
        cReStarts
        cReEnds
        
        cZetaStarts
        cZetaEnds
        
    end
    
    methods
        function this = HybridLinearProgramDecVarIndexes()
            this.phyStarts = [];
            this.phyEnds = [];
            
            this.lambdaStarts = [];
            this.lambdaEnds = [];
            
            this.reStarts = [];
            this.reEnds = [];
            
            this.wLambdaStarts = [];
            this.wLambdaEnds = [];
            
            this.wReStarts = [];
            this.wReEnds = [];
            
            this.cPsyStarts = [];
            this.cPsyEnds = [];
            
            this.cGuardStarts = [];
            this.cGuardEnds = [];
            
            this.cReStarts = [];
            this.cReEnds = [];
            
            this.cZetaStarts = [];
            this.cZetaEnds = [];
        end
        
    end
end

