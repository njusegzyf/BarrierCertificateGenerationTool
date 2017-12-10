classdef Partition2
  %PARTITION Represents a partition of P and C.
  
  properties
    pBoundLow
    pBoundHigh
    cBoundLow
    cBoundHigh
  end
  
  methods
    % Constructor
    function this = Partition(pBoundLowArg, pBoundHighArg, cBoundLowArg, cBoundHighArg)
      this.pBoundLow = pBoundLowArg;
      this.pBoundHigh = pBoundHighArg;
      this.cBoundLow = cBoundLowArg;
      this.cBoundHigh = cBoundHighArg;
    end 
    
    function res = splitPartition(this)
      pBoundMiddle = (this.pBoundLow + this.pBoundHigh) / 2;
      cBoundMiddle = (this.cBoundLow + this.cBoundHigh) / 2;
      res = [Partition(this.pBoundLow, pBoundMiddle, this.cBoundLow, cBoundMiddle)...
             Partition(pBoundMiddle, this.pBoundHigh, this.cBoundLow, cBoundMiddle)];
    end
  end
  
  methods (Static = true)  
    function a = test(b, c)  
      a = b + c;  
    end  
  end  
  
end

