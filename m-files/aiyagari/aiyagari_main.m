% Code to solve the Aiyagari model
% Heterogeneous Agents Model with No Aggregate Uncertainty
% Based on Aiyagari (1994, QJE)
% Solution Methods: 
%   (1) Simulation Algorithm
%   (2) Iteration on the distribution
%   
% Tiago Cavalcanti, August 2016

% Clear MATLAB workspace and empty the command window
clear; clc;

%% Declare Parameters
% =======================================================================
% Define global objects, such as parameters and computational objects
% This saves you defining them as function inputs later on
global betta rtp mu theta delta l_grid lnl_grid Pi Na Nl Ni T Tburn a_grid b

% PARAMETERS
% =======================================================================
betta   = 0.96;     % Discount Factor
rtp     = (1/betta) - 1;    % Rate of Time Preference
mu      = 3;        % Coefficient of Relative Risk Aversion
theta   = 0.36;     % Capital Income Share
delta   = 0.08;     % Depreciation Rate

%% Approximate Labour Endowment Shocks with Nl State Markov Chain
% =======================================================================
% Approximate the AR(1) Process
% ln(l(t+1)) = rho*ln(l(t)) + sigma*(1-rho^2)^(1/2)*eps(t)
%   with eps ~ N(0,1)
% The discretised grid should have Nl, equally-spaced nodes.

% The associated NlxNl transition matrix Pi has element pi_{i,j} denoting
% the probability of moving from state i to state j

Nl      = 7;        % Number of Labour Grid Points
rho     = 0.6;
sigmain = 0.4;
sigma   = sigmain*sqrt(1-rho^2);

% Markov Approximation
[Pi,lnl_grid,invdist,alambda,asigma] = markovappr(rho,sigma,3,Nl);

% Convert the Discrete Grid for ln(l) to a corresponding discrete grid for
% the level of labour l - this is important, because household decisions
% are written in terms of levels (not logarithms)
l_grid = exp(lnl_grid);

% Compute the Stationary Distribution for Labour Supply
% =======================================================================
labour = l_grid*invdist;

%% To Convince You It is Really Chain, Let's Simulate It
% =======================================================================
[eg_chain,eg_state] = markovsimul(Pi,lnl_grid,100);

figure('name','AR(1) and Markov Approximation')
plot([1:1:100],exp(eg_chain),'k','LineWidth',2)
xlabel('Simulation of Labour Grid','Interpreter','latex')

clear eg_chain eg_state

%% Some Useful Functions
% =======================================================================
% Use MATLAB's function handler
% Equilibrium Capital Demand as a function of interest rate and labour,
% although the latter is time-invariant, pinned down by the discrete grid
% for labour
KapD = @(r) labour*(theta./(r+delta)).^(1/(1-theta));

%% (1) Solve by Simulation
% =======================================================================
% Set Some Preliminaries for the Solution Method
Na      = 200;      % Number of Grid Points for Individual Asset Holdings
Ni      = 1000;     % Number of Households in Cross-Section
T       = 10000;    % Number of Time Periods in Simulation
Tburn   = 1000;     % Number of Time Periods Discarded as Burn-In

% Set the Debt Limit
b       = 0;        % From Question, Cannot Borrow

% Create a bracket for the Real Interest Rate within which I bisect
minrate = -delta + 0.04;    % I could pick -delta, but I want to save time
maxrate = rtp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REDEFINE THIS NAIVE BRACKET TO SAVE TIME
% minrate = 0.03525; maxrate = 0.03526;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the tolerance and some initial, arbitrarily large error to start
% the bisection
tol     = 1e-6;
err     = 1;
it      = 1;    % Counter

% Start Solving
while err > tol
    disp(['=====================']);
    disp(['Iteration Number ',num2str(it)]);
    % Define an intermediate interest rate in the bracket
    r_g = (minrate + maxrate)/2;
    
    % Calculate Implied Capital Demand
    Kd  = KapD(r_g);
    
    % Solve for the Implied Capital Supply by Simulation
    Ks  = aiyagari_simul(r_g);
    
    % Now Do the Bisection
    if Kd > Ks          % Excess Demand
        minrate = r_g;  % Search for a higher real interest rate
    elseif Kd < Ks      % Excess Supply
        maxrate = r_g;
    end
    err = abs(maxrate - minrate);
    it  = it + 1;
    
    % Get MATLAB to break the bad news
    disp(['Error ',num2str(err),' for interest rate ',num2str(r_g)]);
end

disp(['Equilibrium Interest Rate: ',num2str(r_g*100)])
disp(['Equilibrium Capital Stock: ',num2str(Ks)]);

%% Plot the Policy Functions at this Interest Rate
% =======================================================================
PF = plot_pol(r_g);

%% Simulate the Asset Holdings to Marvel at Histogram
% =======================================================================
Kd_sol = KapD(r_g);
simul_asset(r_g,Kd_sol,PF);

%% (2) Solve by Iterating on Asset Holding Distribution
% =======================================================================
% Set Some Preliminaries for the Solution Method
Na      = 200;      % Number of Grid Points for Individual Asset Holdings

% Set the Debt Limit
b = 0;

% Create a bracket for the Real Interest Rate within which I bisect
minrate = -delta + 0.04;    % I could pick -delta, but I want to save time
maxrate = rtp;

% Define the tolerance and some initial, arbitrarily large error to start
% the bisection
tol     = 1e-6;
err     = 1;
it      = 1;    % Counter

% Start Solving
while err > tol
    disp(['=====================']);
    disp(['Iteration Number ',num2str(it)]);
    
    % Define an intermediate interest rate in the bracket
    r_g = (minrate + maxrate)/2;
    
    % Calculate Implied Capital Demand
    Kd  = KapD(r_g);
    
    % Solve for the Implied Capital Supply by Iterating in the Distribution
    % of Assets
    Ks  = aiyagari_statdist(r_g,0);
    
    % Now Do the Bisection
    if Kd > Ks          % Excess Demand
        minrate = r_g;  % Search for a higher real interest rate
    elseif Kd < Ks      % Excess Supply
        maxrate = r_g;
    end
    err = abs(maxrate - minrate);
    it  = it + 1;
    
    % Get MATLAB to break the bad news
    disp(['Error ',num2str(err),' for interest rate ',num2str(r_g)]);
end

disp(['Equilibrium Interest Rate: ',num2str(r_g*100)])
disp(['Equilibrium Capital Stock: ',num2str(Ks)]);

%% Plot Capital Demand and Supply Curves
% =======================================================================
b = 0;
minrate = -delta+0.04;
maxrate = rtp;
rate = minrate:(maxrate-minrate)/19:maxrate;
for q = 1:length(rate)
    k1(q) = aiyagari_statdist(rate(q),0);
end

b = 5;
for q = 1:length(rate)
    k2(q) = aiyagari_statdist(rate(q),0);
end

rgrid = minrate:0.005:maxrate;
Kd = labour*(theta./(rgrid+delta)).^(1/(1-theta));

horaxis = [-5:1:35];
H = length(horaxis);

figure('name','Capital Supply and Demand')
hax = axes;
plot(k1,rate,'b','LineWidth',2); hold on;
plot(k2,rate,'r--','LineWidth',2); hold on;
plot(Kd,rgrid,'k','LineWidth',2); 
legend('Capital Supply, b=0','Capital Supply, b=5','Capital Demand',...
    'FontSize',8,'Location','East','Orientation','Vertical','Interpreter','latex'); 
hold on;
plot(horaxis,rtp*ones(H),'k--','LineWidth',0.5); hold on;
plot(horaxis,zeros(H),'k','LineWidth',0.5); hold on;
line([0 0],get(hax,'YLim'),'Color','k','LineWidth',0.5); hold off;
title('Aiyagari Model','FontSize',12,'Interpreter','latex')
xlabel('Capital Quantity','FontSize',10,'Interpreter','latex')
ylabel('Interest Rate','FontSize',10,'Interpreter','latex')
