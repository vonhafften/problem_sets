// This code solves Hayashi's investment model.
// Written by Akio Ino
// Modified on Nov. 14, 2021 by Alex von Hafften

var k, q, i, z, sdf, c, qob, kob; // Even though z is exogenous shock, it is determined by z'=(1-rho) + rho z + eps, so it is determined inside the model and we should include z in this section.

varexo eps_z, measure_error;


parameters 
	theta r delta psi rho gamma;


// Start describing the model. For our purpose, 7 and 8 are redundant.
model;

    // NOTE: in dynare notation, k = k_{t+1}, and k(-1) = k_t.

    // 1. FOC wrt I_t (equation (4) in lecture note)
    q = 1 + psi * (i - delta * k(-1))/k(-1);

    // 2. FOC wrt k_{t+1} (equation (5))
    q = sdf(+1) * (theta * z(+1) * k^(theta-1) + psi * (i(+1)^2 - delta^2 * k^2) / (2*k^2) + q(+1) * (1-delta) );

    // 3. Law of motion for capital-output
    k = i + (1-delta) * k(-1);

    // 4. law of motion for z
    z = rho * z(-1)+ (1-rho) + eps_z;
    
    // 5. Stochastic discount factor
    sdf = (1/(1+r)) * (c/c(-1))^(-gamma);
	
    // 6. Market clearing for consumption goods.
    c = z * k(-1)^theta - psi * (i - delta*k(-1))^2/(2 *k(-1)) - i ;

    // 7. Observation equation for q
    qob = q + measure_error;
    
    // 8. Observation equation for k. 
    // Since we use the detrended capital whose mean is zero, we subtract the steady state value from k.
    kob = k - (theta/(r + delta))^(1/(1-theta));


end;

// This section gives dynare the initial value to solve for the steady state.
initval;
	z = 1;
	q = 1;
	k = 70;
	i = 10;
	sdf = 0.96;
	c = 9;
    kob = 0;

end;


// Set the correlation btw eps and measurement error to be zero.
shocks;
corr eps_z, measure_error = 0;

end;

// Calibrated parameter
r     = 0.04;   // discount rate
delta = 0.15;   // depreciation rate
rho   = 0.7;    // persistency of productivity shock
gamma = 2.0;

// Observable variable
varobs k, qob;

// Parameters to be estimated
// I computed the standard error of kob and qob, and use it as a prior mean for sigma_z and sigma_q.
estimated_params;
theta, uniform_pdf, , ,0, 1;
psi, uniform_pdf, , ,0, 1;
stderr eps_z, inv_gamma_pdf, 78.1330, inf;
stderr measure_error,inv_gamma_pdf, 0.1624, inf;
end;

// with mode_compute = 6, dynare will not use the hessian at the posterior mode in MCMC part.
// Hence even if the hessian is not positive definite, dynare will not stop the computation,
// although it is more inefficient and takes more time.
// mh_jscale is used to scale the var-cov matrix in MCMC part. 
// It is recommended that we choose mh_jscale so that the acceptance ratio is close to 0.234.
// See, for example, Roberts, Gelman and Gilks(‎1997)"Weak convergence and optimal scaling of random walk Metropolis algorithms".

estimation(datafile=data_ps2,mode_compute = 6, mh_jscale=1.8) k qob;

