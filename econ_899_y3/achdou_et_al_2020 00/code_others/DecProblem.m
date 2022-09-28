% Code written by Eirik Eylands Brandsaas for ECON899 (Computational
% Economics at UW-Madison). Code is based on the code contained at 
% http://www.princeton.edu/~moll/HACTproject.htm, just somewhat simplified

%% Set parameters
clc
clear 
lambda = 0.05;
gamma = 2.5;

rho = 0.018;
r = 0.01;

% grids
Na = 200;
amin=0.0;
amax=2.0;
agrid = linspace(amin,amax,Na)';
Da = agrid(2)-agrid(1);

w =0.01;
z = [1 3];
ygrid = w*z;
Ny=length(ygrid);

% create matrices for the states to make some stuff easier
aa = agrid * ones(1,Ny);
yy = ones(Na,1)*ygrid;

% set the utility functions and it's derivative
util = @(c,gamma) c.^(1-gamma)/(1-gamma);
uprime = @(c,gamma) c.^(-gamma);
uprimeinv = @(dV,gamma) dV.^(-1/gamma);
bc = @(c,y,a) y - c + r*a ;


% numerical parameters
maxit = 15;
crit = 10^(-4);
Delta = 1000;

% preallocate some variables
dVf = zeros(Na,Ny);
dVb = zeros(Na,Ny);
dV0 = zeros(Na,Ny);
cf = zeros(Na,Ny);
cb = zeros(Na,Ny);
c0 = zeros(Na,Ny);
adotf = zeros(Na,Ny);
adotb = zeros(Na,Ny);
If = false(Na,Ny);
Ib = false(Na,Ny);
I0 = false(Na,Ny);

% initial guess (present value of staying put forever)
V0 =util(r.*aa + yy,gamma)/rho;
Vnew = V0;

%%
for n=1:maxit
    %disp(sprintf('Starting iteration %d',n))
    V = Vnew;
    for iy=1:Ny
        dVf(1:Na-1,iy) = (V(2:Na,iy) - V(1:Na-1,iy))/Da;
        dVb(2:Na,iy) = (V(2:Na,iy) - V(1:Na-1,iy))/Da;
        
        % End point corrections, only the first is important
        dVb(1,iy) = uprime(ygrid(iy) + r*amin,gamma);
        dVf(Na,iy) = uprime(ygrid(iy) + r*amax,gamma); 

        cf(:,iy) = uprimeinv(dVf(:,iy),gamma) ;
        cb(:,iy) = uprimeinv(dVb(:,iy),gamma) ;
        adotf(:,iy) = bc(cf(:,iy),ygrid(iy),agrid);
        adotb(:,iy) = bc(cb(:,iy),ygrid(iy),agrid);
        
        c0(:,iy) = ygrid(iy) + r*agrid;
        dV0(:,iy) = uprime(c0(:,iy),gamma);
        
        If(:,iy) = adotf(:,iy)>0; %10^(-6);
        Ib(:,iy) = adotb(:,iy)<0; %-10^(-6);
        % I0(:,iy) = 1 - If(:,iy) - Ib(:,iy); (Using this line you instead
        % of the line outside the loop can lead to imaginary numbers, no
        % idea why...
    end
    I0 = (1-If-Ib);
    
   
    dVupwind = (dVf.*If + dVb.*Ib + dV0.*I0);
    c = uprimeinv(dVupwind,gamma);
    adot = bc(c,yy,aa);
    u = util(c,gamma);
    
    % Construct the A matrixIb
    Xvec = - Ib.*adotb/Da;
    Yvec = Ib.*adotb/Da - If.*adotf/Da - lambda;
    Zvec = If.*adotf/Da;
    
    
    A1block = spdiags(Yvec(:,1),0,Na,Na) + spdiags(Xvec(2:Na,1),-1,Na,Na) + spdiags([0;Zvec(1:Na-1,1)],1,Na,Na);
    A2block = spdiags(Yvec(:,2),0,Na,Na) + spdiags(Xvec(2:Na,2),-1,Na,Na) + spdiags([0;Zvec(1:Na-1,2)],1,Na,Na);
    lambdablock = lambda*speye(Na,Na);
    A = [A1block,lambdablock; lambdablock, A2block];
    
    B = (rho + 1/Delta)*speye(2*Na) - A;
    
    ustack = [u(:,1); u(:,2)];
    Vstack = [V(:,1); V(:,2)];
    
    b = ustack + Vstack/Delta;
    Vstack = B\b ;
    Vnew = [Vstack(1:Na), Vstack(Na+1:2*Na)];
    
    diff = max(max(abs(Vnew - V)));
    if diff<crit
        fprintf('Value function converged on iteration %d, distance %f \n',n,diff);
        diff
        break
    end
end
V=Vnew;

It = If + Ib;



%% Distribution
AT= A';
tempvec = zeros(Na*2,1);

% Need to hack one value so that it's not singular
ihack = 1;
tempvec(ihack) = 0.1;
row = zeros(1,Na*2);
row(ihack) = 1;
AT(ihack,:) = row;

gstack = AT\tempvec;
gmass = ones(1,2*Na)*gstack*Da;
gstack = gstack/gmass;

g = [gstack(1:Na), gstack(Na+1:2*Na)];
%% Next lets plot the results

figure(1)
subplot(1,3,1)
plot(agrid,V)
legend('Low Income','High Income','Location','SE')
title("Value function")

subplot(1,3,2)
plot(agrid,adot)
title("Savings adot")
hold on
plot(agrid,zeros(Na),'k')
hold off

subplot(1,3,3)
line(agrid,g)
title("Conditional Distribution")

% 
% subplot(2,2,4)
% bar(agrid,sum(g,2))
% title("Wealth Distribution")
fig.PaperPositionMode = 'auto';
print("tabfig/benchmarksol",'-depsc2','-r0')



    
    
    
