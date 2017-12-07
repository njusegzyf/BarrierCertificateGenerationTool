function degree = computeDegree(expr, indvars)
% ����һ������ʽ�Ĵ���
% expr�������ű����Ķ���ʽ���ʽ
% indvars�������������ɵ�����
% degree����ʾ����ʽ������������

% ����mupad��poly2list����������ʽ�ֽ�Ϊ��������ʽԪ�ء�
% ÿ��Ԫ������һ��1 �� 2����������һ��Ԫ�ر�ʾ�õ���ʽ��ϵ�����ڶ���Ԫ�ر�ʾ�õ���ʽ����ÿ�����ű����Ĵ�����
% �ڵõ����е���ʽ�󣬼����ÿ������ʽ�Ĵ����������ж��������Ĵ���֮�ͣ�����ȡ���Ĵ�����Ϊ��������ʽ�Ĵ�����

% charindvars = converttochar(indvars);
% Note: In Matlab 9.3.0.713579 (R2017b), `converttochar` returns a string
% like `(matrix([[x1, x2]]`, which does not work.

import lp4util.symbolArrayToString
charindvars = symbolArrayToString(indvars);

coefmon = feval(symengine, 'poly2list', expr, charindvars);

degree = 0;

for k = 1:1:length(coefmon)
    dummyvar = reshape(coefmon(k), 2, 1);
    
    % `mon(k,:)` defines the k-th row of the matrix mon.
    % mon(k,:) = double(dummyvar(2));
    % if (sum(mon(k,:)) > degree)
    %     degree = sum(mon(k,:));
    % end
    
    % changed to
    tempDegree = sum(double(dummyvar(2)));
    if (degree < tempDegree)
        degree = tempDegree;
    end
end

end

