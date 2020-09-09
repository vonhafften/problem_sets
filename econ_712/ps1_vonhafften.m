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
p0 = d/r; % initial price
shock = 10;
p0_neg = d/r-shock; % initial price with negative shock
p0_pos = d/r+shock; % initial price with positive shock
dim = 99; % terminal period t = 99

%%%%%%%%%%%%%%%%%
% INITIALIZATION
%%%%%%%%%%%%%%%%%
pvector = zeros(dim+1,1); % creating a vector of price from t=0 to t=99
pvector_neg = zeros(dim+1,1);
pvector_pos = zeros(dim+1,1);

tvector = linspace(0,dim,dim+1)'; % creating a vector for time from 0 to 99

pvector(1) = p0; % giving value to the first element of the price 
                       % vector with the initial price
pvector_neg(1) = p0_neg;
pvector_pos(1) = p0_pos;

%%%%%%%%%%%
% DYNAMICS
%%%%%%%%%%%

for n = 2:dim+1 % starting from t = 1 to t = 100 
    % updating the price in the 
    % next period with the first-
    % order difference equation
    pvector(n) = (1+r)*pvector(n-1)-d;
    pvector_neg(n) = (1+r)*pvector_neg(n-1)-d;
    pvector_pos(n) = (1+r)*pvector_pos(n-1)-d;
end

%%%%%%%%
% PLOTS
%%%%%%%%
figure();
hold on;
plot(tvector(:),pvector_pos(:), '--');
plot(tvector(:),pvector(:));
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
fomc_annc = 20; % period of FOMC announcement
fomc_chng = 50; % period of rate change
dim = 99; % terminal period t = 99
r_pre = 0.01; % risk-free rate before FOMC
r_post = 0.02; % risk-free rate after FOMC


%%%%%%%%%%%%%%%%%
% INITIALIZATION
%%%%%%%%%%%%%%%%%
rvector = [r_pre*ones(1, fomc_chng-1) r_post*ones(1, dim-fomc_chng+1)];

pvector = zeros(dim+1,1); % creating a vector of price from t=0 to t=99

tvector = linspace(0,dim,dim+1)'; % creating a vector for time from 0 to 99

pvector(1) = d/rvector(1); 
                       % giving value to the first element of the price 
                       % vector with the initial price
pvector(fomc_chng) = d/rvector(fomc_chng);
                       % giving the steady state price on the day of the 
                       % rate change

%%%%%%%%%%%
% DYNAMICS
%%%%%%%%%%%

% Steady state before FOMC announcement
for n = 2:fomc_annc-1
    pvector(n) = (1+rvector(n-1))*pvector(n-1)-d;
end

% Transition period between announcement and rate change
for n = fomc_chng-1:-1:fomc_annc
    pvector(n) = (d+pvector(n+1))/(1+rvector(n));
end

% Steady state after rate change
for n = fomc_chng:dim+1
    pvector(n) = (1+rvector(n-1))*pvector(n-1)-d;
end


%%%%%%%%
% PLOTS
%%%%%%%%
figure();
plot(tvector(:),pvector(:));
hold on;
title('Price Dynamics');
xlabel('Time t'); ylabel('Price P_t');
legend('P_0 = 100, P_{50}=50','Location','Northeast');
axis([0 dim 0 120])
