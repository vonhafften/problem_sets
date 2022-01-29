function T = static_resid_tt(T, y, x, params)
% function T = static_resid_tt(T, y, x, params)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%
% Output:
%   T         [#temp variables by 1]  double   vector of temporary terms
%

assert(length(T) >= 6);

T(1) = y(1)^(params(1)-1);
T(2) = params(3)^2;
T(3) = y(1)^2;
T(4) = params(4)*(y(3)^2-T(2)*T(3));
T(5) = y(1)^params(1);
T(6) = params(4)*(y(3)-params(3)*y(1))^2;

end
