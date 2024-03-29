clear
clc

% Sample Solution code of Problem Set 2, ECON712
% Written by Anson Zhou
% Modified by Alex von Hafften

%% Problem 1
% Parameter setup
z = 1;
alpha = 0.3;
delta = 0.1;
beta = 0.97;

% steady state
k_bar = (((beta^-1)-1+delta)/(z*alpha))^(1/(alpha-1)); 
c_bar = z * k_bar^alpha - delta * k_bar;

% compute the Jacobian matrix at steady state
J_11 =  z * alpha * k_bar^(alpha-1) + 1 - delta;
J_12 =  -1;
J_21 =  z * alpha * beta * c_bar * (alpha - 1)*(z * k_bar^alpha+(1-delta)*k_bar-c_bar)^(alpha - 2)*(z * alpha * k_bar^(alpha-1) + 1 - delta);
J_22 =  (1 - delta)*beta + z * alpha * beta * ((z*k_bar^alpha + (1 - delta)*k_bar-c_bar)^(alpha-1)-c_bar*(alpha - 1)*(z*k_bar^alpha + (1 - delta)*k_bar - c_bar)^(alpha-2));
J = [J_11, J_12; J_21, J_22];

% eigenvalues and eigenvectors
[V,D]=eig(J);

% new steady state
z = z + 0.1;
k_bar2 = (((beta^-1)-1+delta)/(z*alpha))^(1/(alpha-1));
c_bar2 = z * k_bar2^alpha - delta * k_bar2;

% Jacobian matrix at new steady state
J2_11 = z * alpha * k_bar2^(alpha-1) + 1 - delta;
J2_12 = -1;
J2_21 = z * alpha * beta * c_bar2 * (alpha - 1)*(z * k_bar2^alpha+(1-delta)*k_bar2-c_bar2)^(alpha - 2)*(z * alpha * k_bar2^(alpha-1) + 1 - delta);
J2_22 = (1 - delta)*beta + z * alpha * beta * ((z*k_bar2^alpha + (1 - delta)*k_bar2-c_bar2)^(alpha-1)-c_bar2*(alpha - 1)*(z*k_bar2^alpha + (1 - delta)*k_bar2 - c_bar2)^(alpha-2));
J2 = [J2_11, J2_12; J2_21, J2_22];

% eigenvalues and eigenvectors
[V2,D2]=eig(J2);


%% determine coefficients for specific solution
t0 = 5;
m2 = (k_bar - k_bar2)/(V2(1,2)*D2(2,2)^t0);

% consumption response at t0
c_t0 = c_bar2 + V2(2,2)*m2*(D2(2,2).^t0);

% trajectory using linearized saddle path
Prd = 20-t0+1;  % number of periods we need to calculate

% Hereafter I use Trj (a 2 by Prd matrix) to store the values of capital
% and consumption for the transition period. First row is capital, second
% row is consumption
Trj1 = zeros(2,Prd);
Trj1(1,1) = k_bar;   % At t0, agent are stuck with k_bar
Trj1(2,1) = c_t0;     % They can choose consumption so that they will fall onto the new saddle path
for i = 2:Prd
    Trj1(:,i) = J2*(Trj1(:,i-1) - [k_bar2, c_bar2]')+ [k_bar2, c_bar2]'; 
end


% trajectory using original equation and the "jump" of c_5 calculated under
% linearization 
Trj2 = zeros(2,Prd);
Trj2(1,1)=k_bar;
Trj2(2,1)=c_t0;
for i = 2:Prd
    k = Trj2(1,i-1);
    c = Trj2(2,i-1);
    Trj2(1,i) = z*k^alpha + (1-delta)*k-c;
    Trj2(2,i) = beta*c*(1-delta+z*alpha*(z*k^alpha+(1-delta)*k-c)^(alpha-1));
end

% The matrix "before" is simply storing the system values up to date t0, it
% is used to graph the whole system from date 0 to date 20
before = repmat([k_bar; c_bar],1, t0-1);
Trj1_whole = [before, Trj1];   % Transition under linear hypothesis
Trj2_whole = [before, Trj2];   % Actual transition using "jump" from linearization

time = 1:20;

% Linear plot
figure
subplot(2,1,1)       % add first plot in 2 x 1 grid
plot(time,Trj1_whole(1,:),'-*')
xlabel('Time')
title('k_t')

subplot(2,1,2)       % add first plot in 2 x 1 grid
plot(time,Trj1_whole(2,:),'-*')
xlabel('Time')
title('c_t')


% Showing linear value under actual dynamics does not converge
figure
subplot(2,1,1)       % add first plot in 2 x 1 grid
plot(time,Trj2_whole(1,:),'-*')
xlabel('Time')
title('k_t')

subplot(2,1,2)       % add first plot in 2 x 1 grid
plot(time,Trj2_whole(2,:),'-*')
xlabel('Time')
title('c_t')

% This part finds the actual c_t0 needed using "shooting method". The idea
% is to try different values of c_t0 and choose the one such that the
% system will converge to the new steady state after a long period of time

M = 10000;  % Picking M number of candidates from c_t0
C_range = linspace(c_bar, c_bar2, M);  % Choosing candidates of c_t0
Distance = zeros(1,M);  % The matrix that stores the distance between the new steady state
% and where the system is at after a long period of time following
% particular c_t0
for i = 1:M
    Distance(1,i)=ps2_calib( C_range(i), alpha, beta, delta, z, k_bar, k_bar2, c_bar2);
    % The ps2_calib.m provided is a function that calibrates where the system
    % is at after 50 periods
end
[DD,I] = min(Distance);     % Picking the particular c_t0 that minimizes the distance
c_star = C_range(I);

% Calibrating the actual nonlinear transition using c_star chosen above.
Trj3 = zeros(2,Prd);
Trj3(1,1)=k_bar;
Trj3(2,1)=c_star;
for i = 2:Prd
    k = Trj3(1,i-1);
    c = Trj3(2,i-1);
    Trj3(1,i) = z*k^alpha + (1-delta)*k-c;
    Trj3(2,i) = beta*c*(1-delta+z*alpha*(z*k^alpha+(1-delta)*k-c)^(alpha-1));
end

Trj3_whole = [before, Trj3];  

% Actual transition using "shooting algorithm"
% Actual transition, using c_5 calculated from "shooting algorithm"
figure
subplot(2,1,1)       % add first plot in 2 x 1 grid
plot(time,Trj3_whole(1,:),'-*')
xlabel('Time')
title('k_t')

subplot(2,1,2)       % add first plot in 2 x 1 grid
plot(time,Trj3_whole(2,:),'-*')
xlabel('Time')
title('c_t')

