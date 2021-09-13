% PROGRAM NAME: 387vfigrowth.M
% This program generates the value function and decision rules for
% a nonstochastic growth model.
% Date: 2/17/03

% PARAMETERS
b=.99; %discount factor 
d=0.025; %depreciation rate
a=.36; %capital share

% ASSET VECTOR
klb=0.01; %lower bound of grid points
inc=0.025; %increments
kub=45;%upper bound of grid points
k=[klb:inc:kub];% asset (row) vector
N=size(k);
N=N(1,2);
c1=ones(N,1); % column vector of ones
K=c1*k; % rows are k, cols are k'

% TABULATE CURRENT RETURN (UTILITY) FUNCTION
cs=zeros(N,N); %rows are k, cols are k'
cs=K'.^a-(K-K'*(1-d));%cons
is=find(cs<0);
cs(is)=0;
us=log(cs);
t=isinf(us);
j=find(t==1);
us(j)=-realmax;

% TABULATE INITIAL VALUE FUNCTION GUESS
visr=us(:,1)'; %choose utility associated with k'=0 

pcntol=1; %tolerance for value function iteration
n=1; %if want to run vfi for a set number of iterations
while pcntol >.0001;
   vis=c1*visr; %generates future value function matrix from above row vector
   
   %CONSTRUCT TOTAL RETURN FUNCTION
   wis=us+b*vis;
   
   %CHOOSE HIGHEST VALUE (ASSOCIATED WITH k' CHOICE)
   [vsr,I]=max(wis'); %since max gives highest element in each column of a matrix
   n=n+1;
   tol=max(abs(vsr-visr)); %use sup norm for tolerance
   pcntol=tol/abs(vsr(1,N));
   visr=vsr; %update value functions
   end;
save 387vdr vsr I k;
save 387parm b a d N inc klb kub;

plot(k,k([I])-k)% plot change in the decision rule
%plot(k,vsr) % plot value function
