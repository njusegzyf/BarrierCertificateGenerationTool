function [] = cvxSample()

echo on

% n = 100;
% A = randn(0.5*n,n);
% b = randn(0.5*n,1);
% c = randn(n,1);
% d = randn;

% one solved LP example
n = 2;
A = [1, 2;
     1, 2];
b = [10; 10];
c = [3; 2];
d = randn;

cvx_begin

variable x(n)
dual variables y z
% minimize( x(1) )
subject to
y : A * x == b;
z : x >= 0;

cvx_end

echo off

cvx_status
% one of 'Solved', 'Unbounded',...
% @see http://web.cvxr.com/cvx/doc/solver.html

cvx_optval
cvx_optbnd

cvx_cputime

end
