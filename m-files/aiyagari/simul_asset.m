function simul_asset(r,Kd,PF)
% Function to simulate household asset holdings for the interest rate r
% with corresponding capital demand Kd
% Use policy function PF from earlier VFI
global lnl_grid Pi Na a_grid


Ni      = 10000;     % Number of Agents
T       = 1000;   % Length of Simulation
Tburn   = 100;    % Burn In

% Use Simulation with Ni individuals for T periods (Tburn discarded as
% burn-in)

% Simulate a stochastic labour endowment path for each of the Ni
% individuals over T periods
% The empty matrices:
L_POS = nan(Ni,T);  % This stores the position of the shock on the grid
L_R   = nan(Ni,T);  % This stores the quantitative value of the shock
% Fill them:
for i = 1:Ni
    % Simulate the chain for T periods for each individual
    [value,position] = markovsimul(Pi,lnl_grid,T);
    L_POS(i,:) = position';     % Numerical Position
    L_R(i,:)   = exp(value');   % Value of Labour Endowment
    clear value position
end

% Initialise the Cross Section of Households with Initial Asset
% Holdings
% Construct Empty Matrices
A_POS = nan(Ni,T);  % This stores the position on the asset grid
A_R   = nan(Ni,T);  % This stores the quantitative value of the assets
% Fill in the first column (t=1) randomly
A_POS(:,1) = randi(Na,Ni,1);    % Na random integers in an Ni x 1 vector
A_R(:,1)   = a_grid(A_POS(:,1));

% Now iterate over time to calculate how household's asset holdings
% evolve
for t = 2:T             % Loop over time
    for i = 1:Ni        % Loop over households
    % Compute Next Period Asset Holdings
    A_POS(i,t) = PF(A_POS(i,t-1),L_POS(i,t-1));
    A_R(i,t)   = a_grid(A_POS(i,t));
    end
end

figure('name','Histogram of Asset Holdings');
histogram(A_R(:,end)); hold on;
histogram(Kd*ones(100));hold off;
%plot(Kd,10000,'r','LineWidth',2); hold off;