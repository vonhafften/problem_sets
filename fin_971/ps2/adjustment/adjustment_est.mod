// This code solves Hayashi's investment model.
// Written by Akio Ino
// Modified by Alex von Hafften on Nov. 14, 2021

var k, q, i, z, sdf, c, qob; // Even though z is exogenous shock, it is determined by z'=(1-rho) + rho z + eps, so it is determined inside the model and we should include z in this section.

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


end;

// This section gives dynare the initial value to solve for the steady state.
initval;
	z = 1;
	q = 1;
	k = 70;
	i = 10;
	sdf = 0.96;
	c = 9;

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
estimated_params;
theta, uniform_pdf, , ,0, 1;
psi, uniform_pdf, , ,0, 1;
stderr eps_z, inv_gamma_pdf, 0.007, inf;
stderr measure_error,inv_gamma_pdf, 0.007, inf;
end;

estimation(datafile=data_mat,mh_jscale=1.8) k qob;

