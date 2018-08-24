% MATLAB CODE FOR SOLVING THE BROCK-MIRMAN MODEL WITH VALUE FUNCTION
% ITERACTIONS
%
% This is an example of how to solve deterministic models using dynamic
% programming. The model is U=max(sum(b^t*log(c_t))) s.t c_t+k'=Ak^a, k_0
% given
%
% The code asks the user for the number of grid points for the  state space
% discretization, as well as the upper and lower grid points, as
% percentages of the steady-state k_bar
%
%
clear all
%
%
%% SECTION I: Define Parameters and Variables
% 
% Parameters
alpha=0.4;
A=1;
beta=0.98;
% 
% Number of periods for simulation
%
T=30;
%
%% SECTION: Discretise the State Space
% Number of grid points
%
m=input('Enter number of grid points for k (odd number): ');
disp(' ');
k_bar=(1/(alpha*beta*A))^(1/(alpha-1));
% 
% Grid for capital
%
k_l=input('Enter percent of k_bar for lowest grid point (between 0 and 1): ')*k_bar;
disp(' ');
k_h=input('Enter percent of k_bar for highest grid point: ')*k_bar;
disp(' ');
k=(k_l:((k_h-k_l)/(m-1)):k_h)';
%
% Grid for consumption
%
k1=k;
c=ones(m,1)*(A*k.^alpha)'-k1*ones(1,m);
%
% The following sets all non-positive cs equal to zero
%
for i=1:m
   for j=1:m
      if c<=0
         c(i,j)=0.001;
      end
   end
end
%
% Grid for utility
U=log(c);
% 
%% SECTION III: Value Function Iterations
%
check=1;					%	Initializes condition for stopping rule
iter=0;						%	Initializes number of iterations
o=ones(1,m);				%   Auxiliary vector
%
V=zeros(m,1);				%	Initializes value function
TV=(max(U+beta*V*o))';	%	Finds the first step iteration for TV
%
while check > 0.0001
   iter=iter+1
   V=TV;								% Sets V to be the last TV we found
   TV=(max(U+beta*V*o))';		% Finds the new TV 
   check=norm(TV-V)/norm(V);  % Sets the new numerical value for the stopping rule
end
V_star=TV;	% Defines the optimal value function
%
disp('The number of interations for the Bellman equation were: ');
disp('iter');
%
%% SECTION IV: Optimal Policy Function
%
[TV,j]= max(U+beta*V_star*o); % The vector j record the indexes where the
% value function is maximized
%
policy_k=k(j);
% 
% Plot policy function
figure(1)
plot(k, policy_k)
xlabel('k_t')
ylabel('k_t+1')
hold on
%% SECTION V: Compare to the Analytical Solution
%
analytic_pol_k=(alpha*beta*A*(k.^alpha));
%
% Plot both policies in the same figure
%
plot(k, analytic_pol_k,'--')
hold on
% 
plot(k,k,':')
legend('Numerical policy function', 'Analytical policy function', '45 Degree line')
hold off
%
%% SECTION VI: Simulate Time Path for Capital and Other Endogenous Variables
%
k_path=zeros(T,1);
k_path(1)=k_l;			% Initial value for k
%
for t=1:T
   % The following finds the index for which k_t is equal to some  element
   % of k, i.e., finds where k_t is in the grid
   i= find(k==k_path(t));
   % For the given element of the grid, indexed by i, the next k_t+1 will
   % be given by the optimal policy function at (i)
   k_path(t+1)=policy_k(i);
   y_path(t)=A*k_path(t)^alpha;
   c_path(t)=A*k_path(t)^alpha-k_path(t+1);
   ky_path(t)=k_path(t)/y_path(t);
end
%
figure(2)
subplot(2,2,1)
plot(k_path)
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



                                                         
                              





