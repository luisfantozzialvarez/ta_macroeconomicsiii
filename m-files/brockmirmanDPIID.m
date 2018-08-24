%  
% 
%  SOLVES THE BROCK-MIRMAN MODEL WITH UNCERTAINTY
%  WITH VALUE FUNCTION ITERATIONS
%
%  The model is: 
%
%  U = max E(SUM(b^t*log(c_t)))
%  s.t. c_t = A*z*k^a-k'
%  k_0 given
% 
%  The shocks are iid and defined as:
%  z(t) = {0.9, 1.0, 1.1} 
%  with probability vector:
%  P = [0.2, 0.6, 0.2]
%
%  The code asks the user for the number of grid points for the 
%  state space discretization, as well as the upper and lower 
%  grid points, as percentages of the steady state k_bar.
%
%  It produces figures of the value function and simulates for
%  T = 120, then calculates average statistics after N replications
%  of the simulation (N is inserted by user).
% 
%  NOTE: The code uses the function IID to generate the shocks.
%  Make sure that this function exists in the file path.
%

clear all; disp('  ')

% SECTION 1 - DEFINE PARAMETERS AND VARIABLES

% Parameters
alpha  = 0.4;
A      = 5;
beta   = 0.9888;

% Steady states
k_bar  = (1/(alpha*beta*A))^(1/(alpha-1));
y_bar  = A*k_bar^alpha;
c_bar  = y_bar - k_bar;


% Shock values and probability vector
z      = [0.9, 1.0, 1.1]; 
pr     = [0.2, 0.6, 0.2];

% Number of periods for simulation
T      = 120;  

disp(' ********************* VALUE FUNCTION ITERATIONS **************************** ')
disp(' ')

% Number of grid points
m      = input(' Enter number of grid points for k                               :  ');
m      = round(abs(m));
disp('  ');

% Grid for capital

k_l    = input(' Enter percent of k_bar for lowest grid point (between 0 and 1)  :  ')*k_bar;
disp('  ');
if k_l < 0 | k_l >= k_bar
    error(' The first input for the grid needs to be between zero and one.')
end

k_h    = input(' Enter percent of k_bar for highest grid point (between 1 and 2) :  ')*k_bar;
disp('  ');
if k_h < k_bar | k_h >= 2*k_bar
    error(' The first input for the grid needs to be between one and two.')
end

k      = (k_l:((k_h - k_l)/(m-1)):k_h)'; 

% Grid for consumption matrices (with vectorization)
k1     = k;
c1     = ones(m,1)*(A*z(1)*k.^alpha)' - k1*ones(1,m);
c2     = ones(m,1)*(A*z(2)*k.^alpha)' - k1*ones(1,m);
c3     = ones(m,1)*(A*z(3)*k.^alpha)' - k1*ones(1,m);
% *****************************************************
% NOTE: The consumption matrices can also be defined
%       as a 3D array (to economize in typing), but
%       3D arrays are a bit hard to handle in loops
%       like the value function iteration.
% *****************************************************

% the following sets all non-positive c's equal to zero
for i = 1:m
   for j = 1:m
      if c1(i,j) <= 0
         c1(i,j) = 0;              
      end
      if c2(i,j) <=0
         c2(i,j) = 0;
      end
      if c3(i,j) <=0
         c3(i,j) = 0;
      end
   end
end


% Grid for Utility
U1 = log(c1);
U2 = log(c2);
U3 = log(c3);

% SECTION 2 - VALUE FUNCTION ITERATIONS

check  = [1, 1, 1];                 % Initializes condition for stopping rule
iter = 0;                           % Initializes number of iterations
o = ones(1,m);                      % Auxiliary vector 

V1 = zeros(m,1);                    % Initializes value function
V2 = zeros(m,1);
V3 = zeros(m,1);

TV1 = (max(U1 + beta*V1*o))';        % Finds the first step iteration for TV
TV2 = (max(U2 + beta*V2*o))';
TV3 = (max(U2 + beta*V3*o))';

while min(check) > 0.0001
   iter = iter + 1;
   V1  = TV1;                                  
   V2  = TV2;
   V3  = TV3;
   aux = pr(1)*V1+pr(2)*V2+pr(3)*V3;
   TV1 = (max(U1 + beta*aux*o))';
   TV2 = (max(U2 + beta*aux*o))';
   TV3 = (max(U3 + beta*aux*o))';
   check  = [norm(TV1-V1)/norm(V1), norm(TV2-V2)/norm(V2), norm(TV3-V3)/norm(V3)];  
end

disp(' The number of iterations for the Bellman equation were:'); disp(iter);
disp('  ')

% SECTION 3 - OPTIMAL POLICY FUNCTION

% the following uses the last found V (i.e. V*) to find indeces for policy
[TV1,j1] = max(U1 + beta*aux*o);
[TV2,j2] = max(U2 + beta*aux*o);
[TV3,j3] = max(U3 + beta*aux*o);

policy_k  = [k(j1), k(j2), k(j3)];
policy_c  = (A*z'*(k.^alpha)')' - policy_k;
policy_y  = policy_k + policy_c;

% Plot policy functions

disp(' ********************* PLOTTING POLICY FUNCTIONS **************************** ')
disp(' ') 

figure(1)
plot(k1,policy_k)
hold on
plot(k,k,':')
legend('k_1','k_2','k_3','45^0 line')
title('POLICY FUNCTIONS FOR CAPITAL')
hold off

disp(' Inspect figure -- press ENTER to continue')
pause

figure(2)
plot(k1,policy_c)
hold on
plot(k,k,':')
legend('c_1','c_2','c_3','45^0 line')
title('POLICY FUNCTIONS FOR CONSUMPTION')
hold off

disp(' Inspect figure -- press ENTER to continue')
pause

figure(3)
plot(k1,policy_y)
hold on
plot(k,k,':')
legend('y_1','y_2','y_3','45^0 line')
title('POLICY FUNCTIONS FOR OUTPUT')
hold off
   
disp(' Inspect figure -- press ENTER to continue') 
pause
   
% SECTION 4 - SIMULATE TIME PATH FOR CAPITAL

k_path = zeros(T,1);
c_path = zeros(T,1);
y_path = zeros(T,1);
z_path = zeros(T,1);

S = iid(T,pr);        % produces shock realization using function IID
i = round(ceil(m/2)); % produces the initial value for k (taking the middle point of the grid, or the one next to)

for t = 1:T
    for n = 1:length(z)
        if S(t) == n
            k_path(t+1) = policy_k(i,n);
            c_path(t)   = policy_c(i,n);
            y_path(t)   = policy_y(i,n);
            i           = find(k == k_path(t+1));
        end
    end
end

disp('  ')
disp(' ********************* PLOTTING SIMULATED TIME PATH ************************* ')
disp('  ')
% Plot a simulation
figure(5)
plot([k_path(2:T), c_path(2:T), y_path(2:T)])
hold on
plot([k_bar*ones(T-1,1), c_bar*ones(T-1,1), y_bar*ones(T-1,1)],':')
title('Simulation of time path')
legend('k_t','c_t','y_t')
hold off

disp(' Inspect figure -- press ENTER to continue')
disp('  ')
pause

% SECTION 5 - REPEATING THE SIMULATIONS TO GET STATISTICS

%warning off  % this avoids warning messages for the repetition of simulations

disp(' ********************* GETTING STATISTICS *********************************** ')
disp('  ')
N  = input(' Enter number for replications of the simulation                  :  ');
N  = round(abs(N));

sd = zeros(N,3);
rd = zeros(N,3);

for q = 1:N
    k_path_s = zeros(T,1);
    c_path_s = zeros(T,1);
    y_path_s = zeros(T,1);
    z_path_s = zeros(T,1);

    S        = iid(T,pr);        % produces shock realisation using function IID
    i        = round(ceil(m/2)); % produces the initial value for k (taking the middle point of the grid)
                                 % or the one next to it

    for t = 1:T
        for n = 1:length(z)
            if S(t) == n
                k_path_s(t+1) = policy_k(i,n);
                c_path_s(t)   = policy_c(i,n);
                y_path_s(t)   = policy_y(i,n);
                i             = find(k == k_path_s(t+1));
            end
        end
    end
    
    sd(q,:) = std([k_path_s(2:T), c_path_s(2:T), y_path_s(2:T)]); % remove initial value
    rd(q,:) = sd(q,:)/sd(q,3);
end

average_sd = mean(sd);
average_rd = mean(rd);

disp('  ')   
disp(' The average standard deviations for the three variables are:')
disp('        k        c        y                                  ') 
disp(average_sd)
disp('  ')
disp(' The relative standard deviations for the three variables are:')
disp('        k        c        y                                   ') 
disp(average_rd)

%warning on  % this turns on the warnings after completing the simulations