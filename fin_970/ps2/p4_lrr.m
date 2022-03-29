%% Solution of the Long Run risks Model 
% Produced by Ivan Shaliastovich
% !!! For Teaching Purposes Only !!!
% Modified by Alex von Hafften on March 27, 2022

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
sigma_w=0.0000028;
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

while dif>1.0e-10
    
    
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
T = 12*5+1;
B_c0 = zeros(T,1);
B_x = zeros(T,1);
B_s = zeros(T,1);

B_c0(1) = m_c0 + 0.5*Lambda_s^2*sigma_w^2;
B_x(1) = m_x;
B_s(1) = m_s + 0.5*(ga^2 +  varphi_e^2*Lambda_x^2);

disp('Question 1: ')
fprintf(' Average annualized log return on consumption asset is  %2.5f  \n',  (-log(kappa1) + mu_c)*100*freq);
fprintf(' Average annualized one-period real risk-free  rate is  %2.5f  \n',  (-B_c0(1) - B_s(1)*sigma_c0^2)*100*freq);

% SDF volatility decomposition
c_vol = Lambda_g^2*sigma_c0^2;
x_vol = Lambda_x^2*varphi_e^2*sigma_c0^2;
s_vol = sigma_w^2*Lambda_s^2;
sdf_vol = c_vol + x_vol + s_vol;

disp('Question 2: ')
fprintf(' Average conditional volatility of SDF %2.5f  \n',  sdf_vol*100*freq);
fprintf(' Fraction from short-run news %2.5f  percent \n',  c_vol/sdf_vol*100);
fprintf(' Fraction from long-run news %2.5f  percent \n',  x_vol/sdf_vol*100);
fprintf(' Fraction from volatility news %2.5f  percent \n',  s_vol/sdf_vol*100);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------ Term Structure

for t = 2:T
    B_c0(t) = - m_c0 + B_c0(t-1) + B_s(t-1)*sigma_c0^2*(1+nu) ...
        - 1/2*(Lambda_s*sigma_w - B_s(t-1)*sigma_w)^2;
    B_x(t) = - m_x + B_x(t-1)*rho;
    B_s(t) = - m_s + B_s(t-1)*nu - 1/2*(Lambda_g + (Lambda_x*varphi_e + B_x(t-1)*varphi_e)^2);
end

yield_curve = (B_c0 + B_s*sigma_c0^2)./((1:T)')*100;

% for question 4
plot((1:T)'./freq,yield_curve)
title("Yield Curve")
xlabel("Years to Maturity")
ylabel("Interest Rate")
saveas(gcf,"p4_output_4.png")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------ Model solution: Dividend-Paying Asset


% Solve for kappa1d
dif = 100;
kappa1d = delta;
count = 0;

while dif>1.0e-10
    
    H_x = (m_x + phi)/(1 - kappa1d*rho);

    H_s = (m_s + 0.5 * (Lambda_g - pi_d)^2 + 0.5 * (Lambda_x*varphi_e - kappa1d*H_x*varphi_e)^2 + 0.5 * varphi_d^2)/(1-kappa1d*nu);

    kappa1d_NEW = exp(m_c0 + mu_d - H_s*(kappa1d*nu - 1)*sigma_c0^2 + 0.5*(Lambda_s*sigma_w - kappa1d*H_s*sigma_w)^2);
                
    dif = abs(kappa1d-kappa1d_NEW);
    kappa1d = kappa1d_NEW;
    
    count = count+1;
    
    if count>500 
        kappa1d=2; 
        break;
    end
end

if (kappa1d>1 || kappa1d<0.0) 
    error('Error: No Equilibrium Solution');    
end

% Solve for expected price-consumption and price-dividend ratios
E_pc = log(kappa1/(1 - kappa1));
E_pd = log(kappa1d/(1 - kappa1d));

disp('Question 5: ')
fprintf('Average price-consumption ratio is %2.3f  \n',  E_pc);
fprintf('Average price-dividend ratio is %2.3f  \n',  E_pd);
fprintf('Average price-dividend ratio in data is %2.3f  \n',  3.443);

%% Average equity premia

% consumption asset
c_vol = Lambda_g*sigma_c0^2;
x_vol = (kappa1*A_x*varphi_e)*(Lambda_x*varphi_e)*sigma_c0^2;
s_vol = (kappa1*A_s*sigma_w)*(sigma_w*Lambda_s);
rp_c_avg = c_vol + x_vol + s_vol;

disp('Question 6: ')
fprintf('Average equity premium on consumption asset is %2.5f \n',  rp_c_avg*freq*100);
fprintf('Fraction that comes from short-run news is %2.5f percent \n',  c_vol/rp_c_avg*100);
fprintf('Fraction that comes from long-run news is %2.5f percent \n',  x_vol/rp_c_avg*100);
fprintf('Fraction that comes from volatility is %2.5f percent \n',  s_vol/rp_c_avg*100);
disp(' ')

% dividend-paying asset
c_vol = pi_d*Lambda_g*sigma_c0^2;
x_vol = (kappa1d*H_x*varphi_e)*(Lambda_x*varphi_e)*sigma_c0^2;
s_vol = (kappa1d*H_s*sigma_w)*(sigma_w*Lambda_s);
rp_d_avg = c_vol + x_vol + s_vol;

fprintf('Average equity premium on dividend-paying asset is %2.5f \n',  rp_d_avg*freq*100);
fprintf('Fraction that comes from short-run news is %2.5f percent \n',  c_vol/rp_d_avg*100);
fprintf('Fraction that comes from long-run news is %2.5f percent \n',  x_vol/rp_d_avg*100);
fprintf('Fraction that comes from volatility is %2.5f percent \n',  s_vol/rp_d_avg*100);
disp(' ')

% Volatility of equity premia
Sig2V = sigma_w^2/(1-nu^2)^2;
EV2 = sigma_c0^2;
SigV = Sig2V/(4*EV2); 
SigX = varphi_e^2*SigV/(1-rho^2); 

rp_c_vol = (Lambda_g  + (kappa1*A_x*varphi_e)*(Lambda_x*varphi_e))^2*Sig2V;
rp_d_vol = (pi_d*Lambda_g + (kappa1d*H_x*varphi_e)*(Lambda_x*varphi_e))^2*Sig2V;
rf_vol = (m_x)^2*SigX + (m_s + 0.5*(Lambda_g^2 + Lambda_x^2*varphi_e^2))^2*Sig2V;

disp('Question 7: ')
fprintf('Volatility of equity premium on consumption asset is %2.5f \n',  rp_c_vol*freq*100);
fprintf('Volatility of equity premium on dividend-paying asset is %2.5f \n',  rp_d_vol*freq*100);
fprintf('Volatility of risk free rate is %2.5f \n',  rf_vol*freq*100);