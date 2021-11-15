function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
% function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   g1
%

if T_flag
    T = noadjustment.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(6, 12);
g1(1,5)=1;
g1(2,4)=(-(y(11)*params(1)*y(10)*T(5)));
g1(2,10)=(-(y(11)*params(1)*T(1)));
g1(2,11)=(-(params(1)*y(10)*T(1)+1-params(3)));
g1(3,1)=(-(1-params(3)));
g1(3,4)=1;
g1(3,6)=(-1);
g1(4,2)=(-params(4));
g1(4,7)=1;
g1(4,12)=(-1);
g1(5,8)=1;
g1(5,3)=(-(T(2)*T(6)*T(7)));
g1(5,9)=(-(T(2)*T(7)*1/y(3)));
g1(6,1)=(-(y(7)*T(4)));
g1(6,6)=1;
g1(6,7)=(-T(3));
g1(6,9)=1;

end
