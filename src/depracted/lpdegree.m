function lp = lpdegree(lp, degree)
% lpdecvars --- add decision variables.
%
% lp is a linear program. 
% vars are decsion variables and should be symbolic.
%

lp.degree = degree;

p = sym('p',[1, monomialNumber(length(lp.indvars), degree)]);

lp = lpdecvars(lp, p);

%% different from lp2
% %r = sym('r');   
% eps1 = sym('eps1');
% eps2 = sym('eps2');
% 
% lp = lpdecvars(lp, [eps1,eps2]);

%%

lp.phy = p * monomials(lp.indvars,[0:degree]);

end

