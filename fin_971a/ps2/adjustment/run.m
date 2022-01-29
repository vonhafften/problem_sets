clear all;

addpath /Applications/Dynare/4.8-2021-11-09-1819/matlab

% Run dynare to generate the simulated data.
dynare adjustment.mod

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_irf = oo_.irfs.i_eps_z ./ oo_.steady_state(3) .* 100;
k_irf = oo_.irfs.k_eps_z ./ oo_.steady_state(1) .* 100;

plot(i_irf)
title("Investment IRF")
xlabel("Periods from shock")
ylabel("Percent deviation from steady state")
saveas(gcf,"figures/p2_i.png")

plot(k_irf)
title("Capital IRF")
xlabel("Periods from shock")
ylabel("Percent deviation from steady state")
saveas(gcf,"figures/p2_k.png")

% Impulse response if gamma is zero.
dynare adjustment_zero_gamma.mod

i_irf_zero_gamma = oo_.irfs.i_eps_z ./ oo_.steady_state(3) .* 100;
k_irf_zero_gamma = oo_.irfs.k_eps_z ./ oo_.steady_state(1) .* 100;

plot(i_irf_zero_gamma)
title("Investment IRF")
xlabel("Periods from shock")
ylabel("Percent deviation from steady state")
saveas(gcf,"figures/p2_i_zero_gamma.png")

plot(k_irf_zero_gamma)
title("Capital IRF")
xlabel("Periods from shock")
ylabel("Percent deviation from steady state")
saveas(gcf,"figures/p2_k_zero_gamma.png")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load("adjustment/Output/adjustment_results.mat")

% Pull simulation variables
k   = oo_.endo_simul(1,:)'; 
q   = oo_.endo_simul(2,:)'; 
i   = oo_.endo_simul(3,:)';
z   = oo_.endo_simul(4,:)';

% get profit
pi  = z .* k .^ 0.7;

% create x and y
y = i./k;
x = [ones(length(q),1) q pi./k];

% lag x
y = y(2:end);
x = x(1:end-1, :);

% OLS
betas = x\y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Obtain simulated results.
k   = oo_.endo_simul(1,:)'; 
q   = oo_.endo_simul(2,:)'; 
i   = oo_.endo_simul(3,:)';
z   = oo_.endo_simul(4,:)';
sdf = oo_.endo_simul(5,:)'; 
c   = oo_.endo_simul(6,:)'; 

% Add measurement error
TW = length(oo_.endo_simul(2,:));
e  = randn(TW,1) * 0.01;
qob = q + e;

% save the data
save('data_mat.mat','k','qob')

% estimate the parameters using the simulated data.
dynare adjustment_est.mod

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% estimate the parameters using the real data.
dynare adjustment_est_D.mod

