function number = monomialNumber(variableNumber, degree)
%   monomialNumber returns the number of monomials which consists of 'variableNumber' variables and the degree is less than 'degree'
%   ����monimialNumber�������Ǹ��ݱ����ĸ�����degree������ɱ������ɵĴ���������degree�����е���ʽ������
%   ������variableNumber������������ʾ����������
%   degree������������ʾ����ʽ������������
%   ����ֵ��number������ʽ������

number = 0;
for i = 0 : degree
    % nchoosek(variableNumber+i-1, variableNumber-1) ���� �ܴ���Ϊ i ʱ�Ķ���ʽ����
    % �����������������е� N��С���M�����ӣ���ͬ���в�ͬ��������� ���⣬�� ��巨 �ɽ�
    % @see http://chensmiles.blog.163.com/blog/static/121463991200962113136292/
    number = number + nchoosek(variableNumber+i-1, variableNumber-1);
end

end

