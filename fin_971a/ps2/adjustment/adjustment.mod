// This code solves Hayashi's investment model.
// Written by Akio Ino
// Modified on Nov. 14, 2021 by Alex von Hafften

var k, q, i, z, sdf, c; // Even though z is exogenous shock, it is determined by z'=(1-rho) + rho z + eps, so it is determined inside the model and we should include z in this section.

varexo eps_z;


parameters 
	theta r delta psi rho sigma gamma;

// set parameter values

theta = 0.7;    // curvature of profit function
r     = 0.04;   // discount rate
delta = 0.15;   // depreciation rate
psi   = 0.01;   // adjustment cost
rho   = 0.7;    // persistency of productivity shock
sigma = 0.01;   // std(eps)
gamma = 2.0;    // CRRA utility parameter


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


end;

// This section gives dynare the initial value to solve for the steady state.
initval;
	z = 1;
	q = 1;
	k = (theta/(r + delta))^(1/(1-theta));
	i = delta * k;
	sdf = 1/(1+r);
	c = k^theta - i;

end;
// Check if the initial value given above satisfies the steady state condition.
resid(1);

// Determine the size of shock. Here I set sigma as a standard error of eps.
shocks;
    var eps_z; stderr sigma; 
end;

// use Dynare capabilities to generate TeX-files of the dynamic and static model
write_latex_static_model;
write_latex_dynamic_model;

// solve for steady state and use it as an initial value to compute the impulse reponse.
//steady;

// check if the Blanchard-Kahn condition is satisfied or not.
//check;

// compute the impulse reponse with linear approximation (order = 1) for 40 periods (irf = 40).
// Also simulate the model to generate simulated data for 20000 periods (periods = 20000). Drop first 200 periods to avoid the effect of initial condition.
stoch_simul(order=2,irf=40, periods = 200, nodisplay);
