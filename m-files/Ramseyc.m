%% This is a fixed point iteration loop to solve the Hamilton Jacobi BellmanPDE
% for the Neoclassical Growth Model
% Written by Tiago Cavalcanti but based on Ben Moll HJB_NGM.m,
% which is available on Ben Moll' website

clear all; clc;

tic;
%% Parameters
s = 2;          % CRRA coefficient
a = 0.3;        % Capital share in income
d = 0.05;       % Depreciation rate
r = 0.05;       % Subective discount factor
A = 1;          % TFP
%% Grid for the capital stock
kss = (a*A/(r+d))^(1/(1-a));     % Steady-state level of capital

I=150;              % Grid points for capital stock
kmin = 0.001*kss;   % Minimum level of capital stock
kmax = 2*kss;       % Maximum level of capital stock
k = linspace(kmin,kmax,I)'; % Column vector for capital
dk = (kmax-kmin)/(I-1);     % Change in the capital stock in the grid
%% Parameters of the algorithm
maxit=10000;
crit = 10^(-6);

dVf = zeros(I,1);  % 
dVb = zeros(I,1);
c = zeros(I,1);

%INITIAL GUESS
v0 = (A.*k.^a).^(1-s)/(1-s)/r;
v = v0;

for n=1:maxit
    V = v;
    % forward difference
    dVf(1:I-1) = (V(2:I)-V(1:I-1))/dk;
    dVf(I) = 0; %will never be used
    % backward difference
    dVb(2:I) = (V(2:I)-V(1:I-1))/dk;
    dVb(1) = 0; %will never be used
    
    I_concave = dVb > dVf; %indicator whether value function is concave (problems arise if this is not the case)
    
    %consumption and savings with forward difference
    cf = dVf.^(-1/s);
    muf = A.*k.^a - d.*k - cf;
    %consumption and savings with backward difference
    cb = dVb.^(-1/s);
    mub = A.*k.^a - d.*k - cb;
    %consumption and derivative of value function at steady state
    c0 = A.*k.^a - d.*k;
    dV0 = c0.^(-s);
    
    % dV_upwind makes a choice of forward or backward differences based on
    % the sign of the drift    
    If = muf > 0; %below steady state
    Ib = mub < 0; %above steady state
    I0 = (1-If-Ib); %at steady state
    %make sure the right approximations are used at the boundaries
    Ib(1) = 0; If(1) = 1; Ib(I) = 1; If(I) = 0;
    dV_Upwind = dVf.*If + dVb.*Ib + dV0.*I0; %important to include third term
    
    c = dV_Upwind.^(-1/s);
    Vchange = c.^(1-s)/(1-s) + dV_Upwind.*(A.*k.^a - d.*k - c) - r.*V;
       
    %% This is the update
    % the following CFL condition seems to work well in practice
    Delta = .9*dk/max(A.*k.^a - d.*k);
    v = v + Delta*Vchange;
    
    dist(n) = max(abs(Vchange));
    if dist(n)<crit
        disp('Value Function Converged, Iteration = ')
        disp(n)
        break
    end
end
toc;

% Graphs
figure(1)
set(gca,'FontSize',14)
plot(dist,'LineWidth',2)
grid
xlabel('Iteration')
ylabel('||V^{n+1} - V^n||')

Verr = c.^(1-s)/(1-s) + dV_Upwind.*(A.*k.^a - d.*k - c) - r.*V;

set(gca,'FontSize',14)
figure(2)
plot(k,Verr,'LineWidth',2)
grid
xlabel('k')
ylabel('Error in HJB Equation')
xlim([kmin kmax])

kdot = A.*k.^a - d.*k - c;

figure(3)
set(gca,'FontSize',12)
plot(k,V,'LineWidth',2)
grid
xlabel('k')
ylabel('V(k)')
xlim([kmin kmax])

figure(4)
set(gca,'FontSize',14)
plot(k,c,'LineWidth',2)
grid
xlabel('k')
ylabel('c(k)')
xlim([kmin kmax])

figure(5)
set(gca,'FontSize',14)
plot(k,kdot,k,zeros(1,I),'--','LineWidth',2)
grid
xlabel('$k$','FontSize',16,'interpreter','latex')
ylabel('$s(k)$','FontSize',16,'interpreter','latex')
xlim([kmin kmax])

%% Simulate the economy
% Start with a k0
%
T=100;
k_path=zeros(T,1);
y_path=zeros(T,1);
c_path=zeros(T,1);
kdot_path=zeros(T,1);
ky_path=zeros(T,1);
k_path(1)=k(10);			% Initial value for k
%
for t=1:T,
   % Find the indices in the capita grid 
    index1 = max(find(k <= k_path(t)));
    index2 = min(find(k >= k_path(t)));
    diff1=k_path(t)-k(index1);
    diff2=k(index2)-k_path(t);
    if diff1>=diff2
       index=index2;
    end
    if diff1<diff2
       index=index1;
    end 
   y_path(t)=A*k_path(t)^a;
   c_path(t)=c(index);
   ky_path(t)=k_path(t)/y_path(t);
   kdot_path(t)=y_path(t)-d*k_path(t)-c_path(t);
   k_path(t+1)=k_path(t)+kdot_path(t);   
end

figure(6)
subplot(2,2,1)
plot(k_path(1:T-1))
title('Time path for capital')
xlabel('t')
ylabel('k')
%
subplot(2,2,2)
plot(y_path)
title('Time path for output')
xlabel('t')
ylabel('y')
%
subplot(2,2,3)
plot(c_path)
title('Time path for consumption')
xlabel('t')
ylabel('c')
%
subplot(2,2,4)
plot(ky_path)
title('Time path for capital to output ratio')
xlabel('t')
ylabel('k/y')
% End


