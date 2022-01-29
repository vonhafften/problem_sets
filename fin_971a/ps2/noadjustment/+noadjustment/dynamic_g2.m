function g2 = dynamic_g2(T, y, x, params, steady_state, it_, T_flag)
% function g2 = dynamic_g2(T, y, x, params, steady_state, it_, T_flag)
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
%   g2
%

if T_flag
    T = noadjustment.dynamic_g2_tt(T, y, x, params, steady_state, it_);
end
g2_i = zeros(14,1);
g2_j = zeros(14,1);
g2_v = zeros(14,1);

g2_i(1)=2;
g2_i(2)=2;
g2_i(3)=2;
g2_i(4)=2;
g2_i(5)=2;
g2_i(6)=2;
g2_i(7)=2;
g2_i(8)=5;
g2_i(9)=5;
g2_i(10)=5;
g2_i(11)=5;
g2_i(12)=6;
g2_i(13)=6;
g2_i(14)=6;
g2_j(1)=40;
g2_j(2)=46;
g2_j(3)=112;
g2_j(4)=47;
g2_j(5)=124;
g2_j(6)=119;
g2_j(7)=130;
g2_j(8)=27;
g2_j(9)=33;
g2_j(10)=99;
g2_j(11)=105;
g2_j(12)=1;
g2_j(13)=7;
g2_j(14)=73;
g2_v(1)=(-(y(11)*params(1)*y(10)*getPowerDeriv(y(4),params(1)-1,2)));
g2_v(2)=(-(y(11)*params(1)*T(5)));
g2_v(3)=g2_v(2);
g2_v(4)=(-(params(1)*y(10)*T(5)));
g2_v(5)=g2_v(4);
g2_v(6)=(-(params(1)*T(1)));
g2_v(7)=g2_v(6);
g2_v(8)=(-(T(2)*(T(7)*(-((-y(9))*(y(3)+y(3))))/(y(3)*y(3)*y(3)*y(3))+T(6)*T(6)*T(8))));
g2_v(9)=(-(T(2)*(T(7)*(-1)/(y(3)*y(3))+T(6)*1/y(3)*T(8))));
g2_v(10)=g2_v(9);
g2_v(11)=(-(T(2)*1/y(3)*1/y(3)*T(8)));
g2_v(12)=(-(y(7)*getPowerDeriv(y(1),params(1),2)));
g2_v(13)=(-T(4));
g2_v(14)=g2_v(13);
g2 = sparse(g2_i,g2_j,g2_v,6,144);
end
