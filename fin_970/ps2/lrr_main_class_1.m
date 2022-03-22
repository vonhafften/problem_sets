%% Solution of the Long Run risks Model 
% Produced by Ivan Shaliastovich
% !!! For Teaching Purposes Only !!!

clear all
close all


freq = 12;      % number of periods per year

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calibration:
% ------------ Preferences:
delta   =0.999;      
psi=  1.5;          
ga=  10;

% ------------ CONSUMPTION:
mu_c= 0.0015;
rho= 0.978; 
varphi_e = 0.037;

%------------ Volatility:
sigma_w=0*0.0000028*1;
nu= 0.999;
sigma_c0=  0.0072;    

% ------------ DIVIDENDS (MARKET Portfolio)
mu_d=  0.0015;     
phi= 2.5;
varphi_d =5.96;
pi_d = 2.6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ------------ Model Solution: Consumption Asset

theta = (1-ga)/(1-1/psi);

% Solve for kappa1
dif = 100;
kappa1 = delta;
count = 0;

while dif>1.0e-006
    
    
    A_x = (1-1/psi)/(1-kappa1*rho);  

    A_s = 0.5*((1-ga)*(1-1/psi) / (1-kappa1*nu))*( 1 +varphi_e^2* (kappa1/ (1-kappa1*rho))^2);
    
    kappa1_NEW = delta*exp( (1-1/psi)* mu_c + A_s*(1-kappa1*nu)*sigma_c0^2 ...
        +0.5*theta*kappa1^2*(A_s )^2*sigma_w^2 );
                
    dif = abs(kappa1-kappa1_NEW);
    kappa1 = kappa1_NEW;
    
    count = count+1;
    
    if count>500 
        kappa1=2; 
        break;
    end
end

if (kappa1>1 || kappa1<0.0) 
    error('Error: No Equilibrium Solution');    
end


% SDF Solution:
m_x = -1/psi;
m_s = (1-theta)*A_s*(1-kappa1*nu);
m_c0 = theta*log(delta) - (theta-1)*log(kappa1) - ga*mu_c - m_s*sigma_c0^2;

Lambda_g = ga;
Lambda_x = (1-theta)*kappa1*A_x;
Lambda_s = (1-theta)*kappa1*A_s;

% Parameters of the Real Yield curve
B_x(1) = -m_x;
B_s(1) = -m_s -0.5*(ga^2 +  varphi_e^2*Lambda_x^2);
B_c0(1) = -m_c0 - 0.5*Lambda_s^2*sigma_w^2;

fprintf(' Average annualized log return on consumption asset is  %2.5f  \n',  (-log(kappa1) + mu_c)*100*freq);
fprintf(' Average annualized one-period real risk-free  rate is  %2.5f  \n',  (B_c0 + B_s*sigma_c0^2)*100*freq);


