classdef Constraint
  %CONSTRAINT 
  
  properties
    num % index
    name
    type % `eq`, `ie`
    polyexpr
    A
    b
  end
  
  methods 
      function res = isEmptyConstraint(this)
          res = strcmp(this.name, 'empty');
      end
  end
  
  methods (Static)
      
      function res = createEmptyConstraint() 
          res = Constraint();
          res.name = 'empty';
      end    
  end
  
end

