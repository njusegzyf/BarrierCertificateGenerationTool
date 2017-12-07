function lp = lpprogram(indvars)
% LPPROGRAM --- Initialize a new linear program.
%
% lp is a new linear program. 
% vars is a vector of independent variables and should be symbolic.
%

if isa(indvars,'sym') 
  % vars can only be a vector of a matrix of symbolic variables
  lp.indvars = reshape(indvars,1,size(indvars,1)*size(indvars,2));
    
  lp.degree = 0;
  lp.phy = 0;
    
  lp.M = 0;
  lp.r = 0;
  lp.eps=[];
    
  lp.f = [];
  lp.decvars = [];
	
  lp.expr.num = 0;
  lp.expr.name = '';
  lp.expr.type = {};
  lp.expr.polyexpr = [];
  lp.expr.A = [];
  lp.expr.b = [];
  
  lp.objective = [];       
	  
end

end

