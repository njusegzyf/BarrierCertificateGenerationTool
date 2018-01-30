function errorIfWrongType(arg, argTypeStr, argName)
% Raises an error if `arg` is not of type `argTypeStr`.
%
% Usage:
% import  lp4util.SolveResBase;
% import  lp4util.CvxSolveRes;
% 
% x = CvxSolveRes(1, 2, 3, 'Solved', 5);
%
% % `x` is of type `lp4util.SolveResBase`, do nothing
% errorIfWrongType(x, 'lp4util.SolveResBase', 'x');
%
% % `x` is not of type `lp4util.SolveResBase`, throw an error
% errorIfWrongType(1, 'lp4util.SolveResBase', 'x');
% ¥ÌŒÛ π”√ errorIfWrongType (line 4)
% Error. x must be a lp4util.SolveResBase, not a double.

if ~isa(arg, argTypeStr)
    error('Error. %s must be a %s, not a %s.', argName, argTypeStr, class(arg));
end

end