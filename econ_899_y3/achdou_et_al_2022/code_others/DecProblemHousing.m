% Code written by Eirik Eylands Brandsaas for ECON899 (Computational
% Economics at UW-Madison). Code is based on the code contained at 
% http://www.princeton.edu/~moll/HACTproject.htm, just somewhat simplified


%% Set parameters
clc
clear 
lambda = 0.5;
gamma = 3;
eta = 0.3;
p = 7;
rho = 0.031;
r = 0.013;
d = 0.3;
hmin = 0.23;
hmax = 1.8;
alpha = 0.5;
% grids
Na = 300;
amin=-0.5;
amax=2;
agrid = linspace(amin,amax,Na)';
Da = agrid(2)-agrid(1);

w =0.01;
z = [1 6];
ygrid = w*z;

%% 
Ny=length(ygrid);

% create matrices for the states to make some stuff easier
aa = agrid * ones(1,Ny);
yy = ones(Na,1)*ygrid;

% set the utility functions and it's derivative
util = @(c,gamma) c.^(1-gamma)/(1-gamma);
uprime = @(c,gamma) c.^(-gamma);
uprimeinv = @(dV,gamma) dV.^(-1/gamma);

func = @(h,eta) -alpha*exp(-eta*h) + exp(0) - r*p*h;
fprime =@(h,eta) alpha*eta*exp(-eta*h);
fprimeinv =@(dV,eta) -1/eta * log(r*p/(eta*alpha));

h = min(fprimeinv(r*p,eta), aa./(d*p));
h = h.*(h>hmin);
h = min(h,hmax);
bc = @(c,f,y,a) y + f +r*a - c;


%%
% numerical parameters
maxit = 30;
crit = 10^(-8);
Delta = 100;

% preallocate some variables
dVf = zeros(Na,Ny);
dVb = zeros(Na,Ny);
dV0 = zeros(Na,Ny);
cf = zeros(Na,Ny);
cb = zeros(Na,Ny);

adotf = zeros(Na,Ny);
adotb = zeros(Na,Ny);
If = false(Na,Ny);
Ib = false(Na,Ny);
I0 = false(Na,Ny);

% initial guess (present value of staying put forever)
V0 = util(yy + func(h,eta) + r.*aa,gamma)/rho;
Vnew = V0;
tic
%%
for n=1:maxit
    %disp(sprintf('Starting iteration %d',n))
    V = Vnew;

    dVf(1:Na-1,:) = (V(2:Na,:) - V(1:Na-1,:))/Da;
    dVb(2:Na,:) = (V(2:Na,:) - V(1:Na-1,:))/Da;

    % End point corrections, only the first is important
    dVb(1,:) = uprime(ygrid + r*amin + func(h(1,:),eta),gamma);
    dVf(Na,:) = uprime(ygrid + r*amax + func(h(Na,:),eta),gamma); 
    
    cf = uprimeinv(dVf,gamma) ;
    cb = uprimeinv(dVb,gamma) ;
        
    adotf = bc(cf,func(h,eta),yy,aa);
    adotb = bc(cb,func(h,eta),yy,aa);
 
    Hf = util(cf,gamma)  + dVf.*adotf;
    Hb = util(cb,gamma)  + dVb.*adotb;
    
    c0 = yy +func(h,eta) +r.*aa;
    Ineither = (1-(adotf>0)) .* (1-(adotb<0));
    Iunique = (adotb<0).*(1-(adotf>0)) + (1-(adotb<0)).*(adotf>0);
    Iboth = (adotb<0).*(adotf>0);
    Ib = Iunique.*(adotb<0) + Iboth.*(Hb>=Hf);
    If = Iunique.*(adotf>0) + Iboth.*(Hf>=Hb);
    
    Ib(1,:) = false;
    I0 = Ineither;
    I0 = (1-If-Ib);
    
    
    dV0 = uprime(c0,gamma) ;%+ fprime(h0,eta);
      
    dVupwind = (dVf.*If + dVb.*Ib + dV0.*I0);
    c = uprimeinv(dVupwind,gamma);
    c = cf.*If + cb.*Ib + c0.*I0;
    adot = bc(c,func(h,eta),yy,aa);
    u = util(c,gamma) ;
    
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
        break
    end
end
%V=Vnew;

AT=transpose(A);


%% Distribution - iterative
% start with uniform
gstack = ones(2*Na,1);
gmass = ones(1,2*Na)*gstack*Da;
gstack = gstack./gmass;
gnew = gmass;
N = 500;
dt = 10;
for i=1:N
    gnew= (speye(2*Na) - AT*dt)\gstack;
    dist = max(abs(gnew-gstack));
    gstack = gnew;
    if dist < crit
        break
    end
end

g2 = [gstack(1:Na), gstack(Na+1:2*Na)];

%% Distribution - iterative 2
% start with uniform
gstack = zeros(2*Na,1);
gstack(1) = 1.0;
gmass = ones(1,2*Na)*gstack*Da;
gstack = gstack./gmass;
gnew = gmass;
N = 500;
dt = 10;
for i=1:N
    gnew= (speye(2*Na) - AT*dt)\gstack;
    dist = max(abs(gnew-gstack));
    gstack = gnew;
    if dist < crit
        break
    end
end

g1 = [gstack(1:Na), gstack(Na+1:2*Na)];

%% Next lets plot the results

astar = p*hmin*d % First time you can afford the house


figure(1)
subplot(2,3,1)

plot(agrid,V)
title("Value function")

hold on;
line([astar astar],[-65,50],'Color','k','LineStyle','--')
hold off;
legend('Low Income','High Income','hmin feasible','Location','SE')
ylim([-65 -45])


subplot(2,3,2)
plot(agrid,adot)
title("Savings adot")

hold on;
plot(agrid,zeros(Na),'k')
line([astar astar],[-0.1,0.1],'Color','k','LineStyle','--')
hold off;
%ylim([-65 -45])

subplot(2,3,3)
plot(agrid,c)
title("Total Consumption (x)")
hold on;
line([astar astar],[0.45,0.65],'Color','k','LineStyle','--')
hold off;

subplot(2,3,4)
plot(agrid,h)
title("Housing")
hold on;
line([astar astar],[0.,1],'Color','k','LineStyle','--')
hold off;


subplot(2,3,5)
plot(agrid,g1)
title("Stationary Distribution #1")
hold on;
line([astar astar],[0.,40],'Color','k','LineStyle','--')
hold off;

subplot(2,3,6)
plot(agrid,g2)
title("Stationary Distribution #2")
hold on;
line([astar astar],[0.,10],'Color','k','LineStyle','--')
hold off;

fig.PaperPositionMode = 'auto';
print("tabfig/benchmarksol_housing",'-depsc2','-r0')
%print("tabfig/benchmarksol_housing",'-depsc2','-r0')

    
    
%% 
%% These parameters work!
% clc
% clear 
% lambda = 0.3;
% gamma = 3.5;
% eta = 1.1;
% p = 8.0;
% rho = 0.045;
% r = 0.015;
% d = 0.38;
% hmin = 0.19;
% alpha = 0.2;
% % grids
% Na = 200;
% amin=0.0;
% amax=1.5;
% agrid = linspace(amin,amax,Na)';
% Da = agrid(2)-agrid(1);
% 
% w =0.011;
% z = [1 7];
% ygrid = w*z;
