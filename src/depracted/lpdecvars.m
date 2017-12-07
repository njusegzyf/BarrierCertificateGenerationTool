function lp = lpdecvars(lp, decvars)
% lpdecvars --- add decision variables.
%
% lp is a linear program. 
% vars are decsion variables and should be symbolic.
%

if isa(decvars,'sym')
    % ���µľ��߱�����Ϊ����������ʽ��ӵ�lp��lp.decvars������
    % decision vars can only be a vector of a matrix of symbolic variables
    % reshape `decvars` to of dim `[1, size(decvars, 1) * size(decvars, 2)]`
    lp.decvars = [lp.decvars reshape(decvars, 1, size(decvars, 1) * size(decvars, 2))];
end

end

