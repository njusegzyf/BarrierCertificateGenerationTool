function lp = lpsolve2(lp)
% lpdecvars --- add decision variables.
%
% lp is a linear program. 
% vars are decsion variables and should be symbolic.
%

Aeq = [];
beq = [];
for k = 1:1:3
    Aeq = [Aeq; lp.expr(k).A];
    beq = [beq; lp.expr(k).b];
end

Aie = [];
bie = [];
for k = 4:1:length(lp.expr)
    Aie = [Aie; lp.expr(k).A];
    bie = [bie; lp.expr(k).b];
end

cvx_begin

    variable x(length(lp.decvars));

    subject to

        Aeq*x==beq
        Aie*x<=bie

cvx_end

% disp('--------------------------------------------------------------');
% disp('The parameter setting:');
% disp(['degree: ',num2str(lp.degree),'; r: ',num2str(lp.r),'; M: ', num2str(lp.M),'; eps1: ',num2str(lp.eps(1)),'; eps2: ',num2str(lp.eps(2))]);
% 
% if (flag == 1)
%     pvar = lp.decvars(1:monomialNumber(length(lp.indvars),lp.degree));
%     pvar_normal = reshape(pvar,1,size(pvar,1)*size(pvar,2));
%     pval = decval(1:monomialNumber(length(lp.indvars),lp.degree));
%     pval_normal = reshape(pval,1,size(pval,1)*size(pval,2));
%     disp('--------------------------------------------------------------');
%     disp('The coefficients of function phy is:');
%     disp(pval);
%     disp('--------------------------------------------------------------');
%     disp('The function phy is:');
%     disp(subs(lp.phy,pvar_normal,pval_normal));
%     disp('--------------------------------------------------------------');
%     disp('The computation time is:');
%     disp(time);
%     disp('--------------------------------------------------------------');
% end
% 
% if (flag == 0)
%     disp('--------------------------------------------------------------');
%     disp('Maximum number of iterations reached.');
%     disp('--------------------------------------------------------------');
% end
% 
% if (flag < 0)
%     disp('--------------------------------------------------------------');
%     disp(['The problem with degree ', num2str(lp.degree),' maybe have no solution.']);
%     disp('--------------------------------------------------------------');
% end

end

