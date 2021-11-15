clear all;

addpath /Applications/Dynare/4.8-2021-11-09-1819/matlab

% Run dynare to generate the simulated data.
dynare noadjustment.mod

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
dynare noadjustment_est.mod
