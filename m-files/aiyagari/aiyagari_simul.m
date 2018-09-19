function [Ks] = aiyagari_simul(r)
% Function to solve for the implied capital supply Ks for a given real
% interest rate r in the Aiyagari (1994) model.
% Solution Method: Simulation

% Redefine the Global Variables from the Master File
global betta rtp mu theta delta l_grid lnl_grid Pi Na Nl Ni T Tburn a_grid b

%% Calculate the Implied Wage for a Given Interest Rate
% =======================================================================
% Because I only solve for this once, I don't bother defining a function
wage = (1-theta)*(theta/(r+delta)^(theta/(1-theta)));

%% Define the Debt Limit
% =======================================================================
% Use the b defined before the function is called
phi = -b;

% If you want a natural borrowing limit:
%if r <= 0
%    phi = b;
%else
%    phi = min(b,wage*l_grid(1)/r);
%end
% -phi is borrowing limit, b is ad hoc, while second term is the natural
% limit

%% Form the Household Asset Holding Grid
% =======================================================================
a_max   = 20;       % Maximum Capital Holdings (needed for computer to solve)
a_min   = phi;      % Minumum Capital Holdings (mathematically defined in problem)

% TRICK: Define the wealth grid logarithmically to put more grid points at
% areas of the state space that are most non-linear (i.e., around the
% borrowing constraint)
a_grid  = (exp(linspace(log(a_min+1-a_min),log(a_max+1-a_min),Na))-1+a_min)'; % Grid for wealth (endogenous state)

% If, instead, you want a vanilla equally-spaced grid
% a_grid = linspace(a_min,a_max,Na);

%% Calculate the Household Decision Rules
% =======================================================================
% I use Value Function Iteration (VFI) to do this.
% For a complete exposition of VFI, see my notes and code for Ph.D. 21,
% Problem Set 3 on the Faculty Website.

% I carry out the VFI in two stages:
% (1) I build a 3-dimensional grid for contemporaneous utility, where
% dimension 1 refers to today's asset holdings, dimension 2 to
% tomorrow's asset holdings and dimension 3 to today's labour endowment
% shock.
% (2) Using the output from (1), I run a vanilla value function
% iteration.

% This approach to VFI is advantageous as complex inequality
% constraints can be easily applied in stage (1) by setting
% contemporaneous utility to arbitrarily large negative values at grid
% points that violate constraints.

%% Build a 3-Dimensional Contemporaneous Utility Grid
% ===================================================================
Ut = zeros(Na,Na,Nl);   % Dimension 1: Today's assets (a); 
                        % Dimension 2: Tomorrow's assets (a');
                        % Dimension 3: Today's labour endowment (l);

% Household Utility
u = @(c) (c^(1-mu)-1)/(1-mu);

for kk = 1:Nl           % Loop Over Labour Supply Endowment Today
    for ii = 1:Na       % Loop Over Savings Today
        for jj = 1:Na   % Loop Over Savings Tomorrow
            l  = l_grid(kk);     % Labour Supply Today
            a  = a_grid(ii);     % Savings Today
            ap = a_grid(jj);     % Savings Tomorrow
            % Solve for Consumption at Each Point
            c = wage*l + (1+r)*a - ap;

            % Apply Constraints
            if (ap < phi)||(c < 0)
                % If Tomorrow's Savings Violate Debt Limit or Today's
                % Consumption is Negative, set Utility to be arbitrarily
                % negative
                Ut(ii,jj,kk) = -99999999999999;
            else
                % If constraints satisfied, calculate contemporaneous
                % utility
                Ut(ii,jj,kk) = u(c);
            end
        end
    end
end

%% Vanilla Value Function Iteration
% ===================================================================
% Initial Guess of Value Function
V0 = kron(l_grid,ones(Na,1));  % l_grid is Nl x 1; ones(Na,1) is Na x 1
                                % V0 is Na x Nl

% Calculate the Guess of the Expected Value Function
for kk = 1:Nl
 EVF(:,kk) = V0*Pi(kk,:)'; 
end

% VFI Tolerance and Convergence Criteria
vfitol  = 0.0001;
vfierr  = 2;
vfiiter = 0;

while vfierr > vfitol
    for kk = 1:Nl   % Loop Over Today's Labour Endowment Shock
        % Here, I use Kronecker products to save a loop over asset holdings 
        % for today and tomorrow 
        Objgrid = Ut(:,:,kk) + betta*kron(ones(Na,1),EVF(:,kk)'); % Kron gives NaxNa
        for ii = 1:Na
            % Loop Over Today's Capital Stock to Find Max of Value Function
            [V1(ii,kk),PF(ii,kk)] = max(Objgrid(ii,:));
        end
    end

    for kk = 1:Nl
        % Loop Over Today's Labour Endowment to Solve for Expected Value
        % Function
        EVF(:,kk) = V1*Pi(kk,:)';
        % To be used as the next guess for the value function
    end

    vfiiter = vfiiter + 1;
    vfierr  = norm(V1(:) - V0(:));
    iter100 = mod(vfiiter,100);
    if iter100 == 0
        display(['Number of VF Iterations ',num2str(vfiiter)]);
    end
    V0   = V1;
end 

%% Build Household Policy Functions for Saving and Consumption
% ===================================================================
% Policy Function for Assets
AF  = zeros(size(PF));
% Policy Function for Consumption
CF  = zeros(size(PF));

for kk = 1:Nl
    for  ii = 1:Na
        a  = a_grid(ii);        % Savings Brought into Today
        ap = a_grid(PF(ii,kk)); % Savings Made Today for Tomorrow
        AF(ii,kk) = ap;
        l  = l_grid(kk);
        CF(ii,kk) = wage*l + (1+r)*a - ap;
    end
end

%% Compute the Stationary Asset Distribution
% ===================================================================
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

%% Calculate Capital Supply in Economy and Update Guess
% ===================================================================
Ks = mean(mean(A_R(:,(Tburn+1):end)));

