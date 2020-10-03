% Pavel Brendler, Anton Babkin, 10/9/15.
% Modified by Jason Choi, 10/18/18.

% Modified by Alex von Hafften (10/3/2020)

% This program evaluates macroeconomic consequences of eliminating social security in the U.S.
% The model is a simplified version of the model by Conesa and Krueger (1999).

% You are asked to understand the logic of the code, complete the missing parts
% and use it to run a policy experiment.
% Look for these comments: % <------------------- YOUR CODE GOES HERE

% This code is not written very efficiently, so it is not fast.
% Instead, it is very explicit, so that you can follow and understand it easier.


% Clear the memory
clear; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Demographics
J=66;                       % life-span
JR=46;                      % age of retirement
tR=J-JR+1;                  % length of retirement
tW=JR-1;                    % length of working life
n=0.011;                    % Population growth

% Preferences
beta=0.97;                  % discount factor
sigma=2;                    % coefficient of relative risk aversion
gamma=0.42;                 % weight on consumption

% Production
alpha=0.36;                 % production elasticity of capital
delta=0.06;                 % rate of depreciation

% Social Security tax rate
tau=0.11; % Use this for model with SS
%tau=0;   % Use this for model without SS


% Measure of each generation
mass=ones(J,1);
for ik0=2:J
    mass(ik0)=mass(ik0-1)/(1+n);
end

% Normalized measure of each generation (sums up to 1)
mass=mass/sum(mass);

% Age-efficiency profile
e = zeros(tW,1);
e(1:11) = linspace(0.6, 1, 11);
e(11:21) = linspace(1, 1.08, 11);
e(21:31) = linspace(1.08, 1.12, 11);
e(31:41) = linspace(1.12, 1.06, 11);
e(41:45) = linspace(1.06, 1.02, 5);

%  Capital grid
%  Be careful when setting maxkap across experiments! This value should not be binding!
%  To see if it's binding, check the optimal capital decision kp(k) of
%  every cohort.
%  If you have problems with convergence, increasing nk might help.

maxkap = 14;                                % maximum value of capital grid
minkap = 0.01;                             % minimum value of capital grid
nk=180;                                    % number of grid points
inckap=(maxkap-minkap)/(nk-1).^2;          % distance between points
aux=1:nk;
kap= minkap+inckap*(aux-1).^2;             % capital grid
% This formula makes a non-uniform grid with more density at lower range.

neg=-1e10;                              % very small number


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial guesses for interest rate, wages and pension benefits
% Comment out the appropriate lines!

% Social Security
K0=3.1392;
L0=0.3496;

% No Social Security
%K0=3.9288;
%L0= 0.3663;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop over capital and labor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tolk=1e-3;              % Numerical tolerance for capital
tollab=1e-3;            % Numerical tolerance for labor
nq=100;                 % Max number of iterations

q=0;                    % Counter for iterations
K1=K0+10;
L1=L0+10;

fprintf('\nComputing equilibrium prices... \n');
while q<nq && (abs(K1-K0)>tolk || abs(L1-L0)>tollab)

    q=q+1;

    fprintf('\nIteration %g out of %g \n',q,nq);

    % Prices
    r0 =  alpha * (L0/K0)^(1-alpha) - delta;
    w0 = (1-alpha) * (K0/L0)^alpha;

    % Pension benefit
    b =  (tau * w0 * L0) / sum(mass(JR:J));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BACKWARD INDUCTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Initialization
    v=zeros(nk,J);                          % value function of agents
    kapopt=ones(nk,J);                      % optimal savings of agents
    % (store INDEX of k' in capital grid, not k' itself!)
    labopt=ones(nk,tW);                     % optimal labor supply


    % Retired households

    % Last period utility
    cons=(1+r0)*kap+b;                    % last period consumption (vector!)
    util=(cons.^(1-sigma))/(1-sigma);       % last period utility (vector!)
    v(:,J)=util;                            % last period indirect utility (vector!)

    for j=J-1:-1:tW+1 % age
        for ik0=1:nk        % assets today

            % Initialize right-hand side of Bellman equation
            vmin=neg;
            ik1=0;

            % Loop over all k's in the capital grid to find the value,
            % which gives max of the right-hand side of Bellman equation

            while ik1<nk  	% assets tomorrow
                ik1=ik1+1;
                kap0=kap(ik0); % current asset holdings
                kap1=kap(ik1); % future asset holdings

                % Instantaneous utility
                cons =  (1 + r0) * kap0 + b - kap1;

                if cons<=0
                    util=neg;
                else
                    util = (cons^(1-sigma)) / (1 - sigma);
                end

                % Right-hand side of Bellman equation
                v0 =  util + beta * v(ik1,j+1);

                % Store indirect utility and optimal saving
                if v0>vmin
                    v(ik0,j)=v0;
                    kapopt(ik0,j)=ik1;
                    vmin=v0;
                end
            end
        end
    end

    % Working households
    for j=tW:-1:1           % age
        for ik0=1:nk          % assets today

            % Initialize right-hand side of Bellman equation
            vmin=neg;
            ik1=0;

            % Loop over all k's in the capital grid to find the value,
            % which gives max of the right-hand side of Bellman equation

            while ik1 < nk  	% assets tomorrow
                ik1=ik1+1;

                kap0=kap(ik0); % current asset holdings
                kap1=kap(ik1); % future asset holdings

                % Optimal labor supply
                lab = (gamma*(1-tau)*e(j)*w0 - (1-gamma)*((1+r0)*kap0 - kap1))/((1-tau)*e(j)*w0);

                % Check feasibility of labor supply
                if lab>1
                    lab=1;
                elseif lab<0
                    lab=0;
                end

                % Instantaneous utility
                cons = (1-tau)*w0*e(j)*lab + (1+r0)*kap0 - kap1;
                
                if cons<=0
                    util = neg;
                else
                    util =  (((cons^gamma)*(1-lab)^(1-gamma))^(1-sigma))/(1 - sigma);
                end

                % Right-hand side of Bellman equation
                v0 =  util + beta * v(ik1,j+1);

                % Store indirect utility, optimal saving and labor
                if v0>vmin
                    v(ik0,j)=v0;
                    kapopt(ik0,j)=ik1;
                    labopt(ik0,j)=lab;
                    vmin=v0;
                end
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Aggregate capital stock and employment                                  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Initializations
    kgen=zeros(J,1);        % capital supply k(j+1) for each generation j
    ikgen=zeros(J,1);       % grid index of k(j+1)
    labgen=zeros(tW,1);     % labor supply l(j) for each generation j

    % Use decision rules to iteratively fill in kgen and labgen
    ik0 = 1;                % starting capital of j = 1, kap(ik0) = 0
    for j=1:J               % iterations over cohort
        % capital decision kp(k)
        ik1 = kapopt(ik0, j);
        ikgen(j) = ik1;
        kgen(j) = kap(ik1);

        % labor decision l(k)
        if j<=tW
            labgen(j) = labopt(ik0, j);
        end

        % update k = kp
        ik0 = ik1;
    end

    K1 = kgen' * mass;              % dot product of vectors
    L1 = (labgen .* e)' * mass(1:tW);      % dot product of vectors

    % Update the guess on capital and labor
    K0=0.9*K0+0.1*K1;
    L0=0.9*L0+0.1*L1;

    % Display results
    disp('  capital     labor   pension');
    disp([K0, L0, b]);
    disp('deviation-capital deviation-labor       ');
    disp([abs(K1-K0),  abs(L1-L0)]);
end

% Display equilibrium results
disp('      K0         L0       w         r         b   ');
disp([K0, L0, w0, r0, b]);

% Check if any cohort wants to save on the upper bound of capital grid
if sum(kgen == kap(end)) > 0
    disp('Capital decision on upper bound, increase!');
end


%% Plots
% Value function for a retired agent
figure(1)
age=50;
plot(kap,v(:,age));
xlabel('asset holdings, k','FontSize',14);
ylabel('value function, V_{50}(k)','FontSize',14);
title(['value function of a retired agent at age ', num2str(age)],'FontSize',14)


% Savings of a working agent
figure(2)
age=20;
plot(kap,kap(kapopt(:,age)),'k-',kap,kap,'r--');
xlabel('asset holdings, k','FontSize',14);
ylabel('saving, k''','FontSize',14);
legend('saving','45 degree line','Location','East');
title(['saving k'' of a working agent at age ', num2str(age)],'FontSize',14)


%% Saving Model Results
if tau ~= 0
% Solve the model with social security and
save('ss.mat');
elseif tau ==0
% Solve the model without social security and
save('no_ss.mat');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compare Model Results %% RUN THIS PART OF THE CODE AFTER YOU'VE SOLVED BOTH SS AND WITHOUT SS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

load('ss.mat','kgen','J');
kgen_ss=kgen;

load('no_ss.mat','kgen','J');
kgen_noss=kgen;

figure(3)
plot1=plot(20+1:20+J,kgen_ss,'b-',20+1:20+J,kgen_noss,'ro');
set(plot1,'LineWidth',1.5)
xlabel('(real-life) age','FontSize',14);
ylabel('wealth','FontSize',14);
AX=legend('with Social Security','without Social Security','Location','NorthWest');
LEG=findobj(AX,'type','text');
set(LEG,'FontSize',14);

% save figure to pdf
set(3, 'PaperSize', [5 5]);
set(3, 'PaperPositionMode', 'manual');
set(3, 'PaperPosition', [0 0 5 5]);
print(3, 'fig_wealth', '-dpdf');

%% Welfare comparison
clear all;

load('ss.mat','v','mass','ikgen','J');

% Newborn generation
V1_ss = v(1,1)

% Aggregate welfare
V_ss = zeros(J,1);
V_ss(1) = v(1,1);
for j = 2:J
    ik0 = ikgen(j-1);
    V_ss(j) = v(ik0,j);
end
W_ss = V_ss' * mass


load('no_ss.mat','v','mass','ikgen','J');

% Newborn generation
V1_noss = v(1,1)

% Aggregate welfare
V_noss = zeros(J,1);
V_noss(1) = v(1,1);
for j = 2:J
    ik0 = ikgen(j-1);
    V_noss(j) = v(ik0,j);
end
W_noss = V_noss' * mass
