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
    T = adjustment_zero_gamma.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(6, 14);
g1(1,1)=(-((y(1)*params(4)*(-params(3))-params(4)*(y(6)-params(3)*y(1)))/(y(1)*y(1))));
g1(1,5)=1;
g1(1,6)=(-(params(4)/y(1)));
g1(2,4)=(-(y(13)*T(10)));
g1(2,5)=1;
g1(2,10)=(-(y(13)*(1-params(3))));
g1(2,11)=(-(y(13)*params(4)*2*y(11)/(2*T(3))));
g1(2,12)=(-(y(13)*params(1)*T(1)));
g1(2,13)=(-(params(1)*y(12)*T(1)+T(4)/(2*T(3))+y(10)*(1-params(3))));
g1(3,1)=(-(1-params(3)));
g1(3,4)=1;
g1(3,6)=(-1);
g1(4,2)=(-params(5));
g1(4,7)=1;
g1(4,14)=(-1);
g1(5,8)=1;
g1(5,3)=(-(T(5)*T(11)*T(12)));
g1(5,9)=(-(T(5)*T(12)*1/y(3)));
g1(6,1)=(-(y(7)*T(8)-(2*y(1)*params(4)*(-params(3))*2*(y(6)-params(3)*y(1))-2*T(7))/(2*y(1)*2*y(1))));
g1(6,6)=(-((-(params(4)*2*(y(6)-params(3)*y(1))/(2*y(1))))-1));
g1(6,7)=(-T(6));
g1(6,9)=1;

end
