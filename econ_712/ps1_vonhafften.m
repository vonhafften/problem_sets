% PROGRAM NAME: ps1
%
% This program generates and plots the price dynamics given the first-order 
% difference equation discussed in the problem set given some initial price 
%
% Prepared by Fu Tan
% Last update: 09/08/15
% Modified by Alex von Hafften on 09/08/2020

clear;
clc;

%%%%%%%%%%%%%%%%%%
% Problem #3
%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%
r = 0.01; % interest rate
d = 1; % constant dividend
p0_eq = 100; % initial price
p0_neg = 90; % initial price
p0_pos = 110; % initial price
dim = 99; % terminal period t = 99

%%%%%%%%%%%%%%%%%
% INITIALIZATION
%%%%%%%%%%%%%%%%%
pvector_eq = zeros(dim+1,1); % creating a vector of price from t=0 to t=99
pvector_neg = zeros(dim+1,1);
pvector_pos = zeros(dim+1,1);

tvector = linspace(0,dim,dim+1)'; % creating a vector for time from 0 to 99

pvector_eq(1) = p0_eq; % giving value to the first element of the price 
                       % vector with the initial price
pvector_neg(1) = p0_neg;
pvector_pos(1) = p0_pos;

%%%%%%%%%%%
% DYNAMICS
%%%%%%%%%%%

for n = 2:dim+1 % starting from t = 1 to t = 100 
    pvector_eq(n) = (1+r)*pvector_eq(n-1)-d; % updating the price in the 
                                             % next period with the first-
                                             % order difference equation
    pvector_neg(n) = (1+r)*pvector_neg(n-1)-d;
    pvector_pos(n) = (1+r)*pvector_pos(n-1)-d;
end

%%%%%%%%
% PLOTS
%%%%%%%%
figure();
plot(tvector(:),pvector_pos(:), '--');
hold on;
plot(tvector(:),pvector_eq(:));
plot(tvector(:),pvector_neg(:), '-.');
title('Price Dynamics');
xlabel('Time t'); ylabel('Price P_t');
legend('P_0 = 110','P_0 = 100','P_0 = 90','Location','Northeast');
axis([0 dim 0 180])


%%%%%%%%%%%%%%
% Problem #4
%%%%%%%%%%%%%%

clear;
clc;

%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%
d = 1; % constant dividend
p0_eq = 100; % initial price
p0_eqq = 69.597;
dim = 99; % terminal period t = 99

%%%%%%%%%%%%%%%%%
% INITIALIZATION
%%%%%%%%%%%%%%%%%
x = (1:2)'/100;
rvector = repmat(x, 1, 50)';
rvector = rvector(:);

pvector_eq = zeros(dim+1,1); % creating a vector of price from t=0 to t=99

tvector = linspace(0,dim,dim+1)'; % creating a vector for time from 0 to 99

pvector_eq(1) = p0_eq; % giving value to the first element of the price 
                       % vector with the initial price

%%%%%%%%%%%
% DYNAMICS
%%%%%%%%%%%

for n = 2:dim+1 % starting from t = 1 to t = 100 
    pvector_eq(n) = (1+rvector(n-1))*pvector_eq(n-1)-d; % updating the price in the 
                                             % next period with the first-
                                             % order difference equation
end

%%%%%%%%
% PLOTS
%%%%%%%%
figure();
plot(tvector(:),pvector_eq(:));
hold on;
title('Price Dynamics');
xlabel('Time t'); ylabel('Price P_t');
legend('P_0 = 100','Location','Northeast');
axis([0 dim 0 220])
