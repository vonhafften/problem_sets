% Modified by Alex von Hafften on Sept 14, 2021

% PROGRAM NAME: 387vfigrowth.M
% This program generates the value function and decision rules for
% a nonstochastic growth model.
% Date: 2/17/03

% PARAMETERS
b=.99; %discount factor 
d=0.025; %depreciation rate
a=.36; %capital share
z_low=.2;
z_high=1.25;
pi_low_low=0.926;
pi_low_high=0.074;
pi_high_low=0.023;
pi_high_high=0.977;

% ASSET VECTOR
klb=0.01; %lower bound of grid points
inc=0.025; %increments
kub=45;%upper bound of grid points
k=[klb:inc:kub];% asset (row) vector
N=size(k);
N=N(1,2);
c1=ones(N,1); % column vector of ones
K=c1*k; % rows are k, cols are k'

% TABULATE CURRENT RETURN (UTILITY) FUNCTION FOR LOW PRODUCTIVITY STATE
cs_low=zeros(N,N); %rows are k, cols are k'
cs_low=z_low*K'.^a-(K-K'*(1-d));%cons
is_low=find(cs_low<0);
cs_low(is_low)=0;
us_low=log(cs_low);
t_low=isinf(us_low);
j_low=find(t_low==1);
us(j_low)=-realmax;

% TABULATE CURRENT RETURN (UTILITY) FUNCTION FOR HIGH PRODUCTIVITY STATE
cs_high=zeros(N,N); %rows are k, cols are k'
cs_high=z_high*K'.^a-(K-K'*(1-d));%cons
is_high=find(cs_high<0);
cs_high(is_high)=0;
us_high=log(cs_high);
t_high=isinf(us_high);
j_high=find(t_high==1);
us(j_high)=-realmax;

% TABULATE INITIAL VALUE FUNCTION GUESS
visr_low=us_low(:,1)'; %choose utility associated with k'=0 
visr_high=us_high(:,1)'; %choose utility associated with k'=0 

pcntol=1; %tolerance for value function iteration
n=1; %if want to run vfi for a set number of iterations

while pcntol >.0001;
   vis_low=c1*visr_low; %generates future value function matrix from above row vector
   vis_high=c1*visr_high; %generates future value function matrix from above row vector
   
   %CONSTRUCT TOTAL RETURN FUNCTION
   wis_low=us_low+pi_low_low*b*vis_low+pi_low_high*b*vis_high; 
   wis_high=us_high+pi_high_low*b*vis_low+pi_high_high*b*vis_high;
   
   %CHOOSE HIGHEST VALUE (ASSOCIATED WITH k' CHOICE)
   [vsr_low,I_low]=max(wis_low'); %since max gives highest element in each column of a matrix
   [vsr_high,I_high]=max(wis_high'); %since max gives highest element in each column of a matrix
   n=n+1;
   tol=max([abs(vsr_low-visr_low), abs(vsr_high-visr_high)]); %use sup norm for tolerance
   pcntol=tol/abs(vsr_low(1,N));
   visr_low=vsr_low; %update value functions
   visr_high=vsr_high; %update value functions
end;

t = tiledlayout(3,1);
nexttile
plot(k,[vsr_low', vsr_high']') % plot value function
title('Value Function')
nexttile
plot(k,[k([I_low])', k([I_high])']) % plot policy function
title('Policy Function')
nexttile
plot(k,[k([I_low])'-k', k([I_high])'-k']')% plot change in the decision rule
title('Policy Function Changes')
exportgraphics(t,'ps1_matlab_figures.png')
