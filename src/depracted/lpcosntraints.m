function lp = lpcosntraints(lp, constraints, name, type)
% lpcosntraints --- add decision variables.
%
% lp is a linear program. 
% vars are decsion variables and should be symbolic.
%

if (nargin == 3)
    expr.num = lp.expr(end).num+1;
    expr.name = name;
    expr.type = 'eq';
    expr.A = [];
	expr.b = [];
    if (strcmp(name,'theta'))
        constraint1 = -lp.phy;  
        de = computeDegree(constraint1,lp.indvars);        
        
        c_alpha_beta = sym('c_alpha_beta',[1,10000]); % pre-defined varibales, only a few of them are the actual variables
        [decvars, expression] = constraintExpression(de,constraints,c_alpha_beta,lp.M);
        lp = lpdecvars(lp,decvars);
             
        constraint1 = constraint1 + expression;
        expr.polyexpr = constraint1;
    end
    if (strcmp(name,'psy'))
        phy_d = 0;
        for k = 1:1:length(lp.indvars)
            phy_d = phy_d + diff(lp.phy,lp.indvars(k))*lp.f(k);
        end
        constraint2 = -phy_d + lp.r*lp.phy + lp.eps(1);   % different from lp2    
        de = computeDegree(constraint2,lp.indvars);
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         de=mande;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        c_gama_delta = sym('c_gama_delta',[1,10000]);
        [decvars, expression] = constraintExpression(de,constraints,c_gama_delta,lp.M);
        lp = lpdecvars(lp,decvars);
         
        constraint2 = constraint2 + expression;
        expr.polyexpr = constraint2;
    end
    if (strcmp(name,'zeta'))
        constraint3 = lp.phy + lp.eps(2);        % different from lp2     
        de = computeDegree(constraint3,lp.indvars);        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         de=mande;
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        c_u_v = sym('c_u_v',[1,10000]);
        [decvars, expression] = constraintExpression(de,constraints,c_u_v,lp.M);
        lp = lpdecvars(lp,decvars);
        
        constraint3 = constraint3 + expression;
        expr.polyexpr = constraint3;
    end
    disp(['constraint ',expr.name,' is generated: ',datestr(now,'yyyy-mm-dd HH:MM:SS')]);
    lp.expr(expr.num).num = expr.num;
    lp.expr(expr.num).name = expr.name;
	lp.expr(expr.num).type = expr.type;
    lp.expr(expr.num).polyexpr = expr.polyexpr;
	lp.expr(expr.num).A = [];
	lp.expr(expr.num).b = [];
end

if (lp.expr(end).num == 3)
    for k = 1:1:length(lp.expr)
        [ lp.expr(k).A, lp.expr(k).b ] = eqgenerate( lp.indvars, lp.decvars, lp.expr(k).polyexpr);
        disp(['constraint ',lp.expr(k).name,' is processed: ',datestr(now,'yyyy-mm-dd HH:MM:SS')]);
    end
    % each c decision variable should be not less than M, -c <= -M
    M = lp.M;
    decexpr.num = lp.expr(end).num+1;
    decexpr.name = 'decvarconstraints';
    decexpr.type = 'ie';
    decexpr.polyexpr = [];
%     Ap = zeros(monomialNumber(length(lp.indvars), lp.degree)); % the p variables have no bounds
%     Ac = -eye(length(lp.decvars)-monomialNumber(length(lp.indvars), lp.degree));
%     decexpr.A = [Ap; Ac];
    decexpr.A = zeros(length(lp.decvars));
    for k = 1:1:length(lp.decvars)
        if (k >= monomialNumber(length(lp.indvars), lp.degree)+1)
            decexpr.A(k,k) = -1;
        end
    end
    bp = zeros(monomialNumber(length(lp.indvars), lp.degree), 1); % the p variables have no bounds
    bc = -ones(length(lp.decvars)-monomialNumber(length(lp.indvars), lp.degree), 1)*lp.M;
    decexpr.b = [bp; bc];
    lp.expr(decexpr.num).num = decexpr.num;
    lp.expr(decexpr.num).name = decexpr.name;
	lp.expr(decexpr.num).type = decexpr.type;
    lp.expr(decexpr.num).polyexpr = [];
	lp.expr(decexpr.num).A = decexpr.A;
	lp.expr(decexpr.num).b = decexpr.b;
end

if (nargin == 4 && strcmp(type, 'eq'))
end

if (nargin == 4 && strcmp(type, 'ie'))
end
    


end

