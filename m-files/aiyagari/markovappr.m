function [Tran,s,probst,alambda,asigmay] = markovappr(lambda,sigma,m,N)
% Approximating first-order autoregressive process with Markov chain
% y_t = lambda * y_(t-1) + u_t

% u_t is a Gaussian white noise process with standard deviation sigma.

% m determines the width of discretised state space, Tauchen uses m=3
% ymax=m*vary,ymin=-m*vary, ymax and ymin are two boundary points

% N is the number of possible states chosen to approximate
% the y_t process, usually N=9 should be fine

% Tran is the transition matrix of the Markov chain

% s is the discretised state space of y_t

% alambda is the theoretical first order autoregression coefficient 
% for Markov chain

% asigma is the theoretical standard deviation for Markov chain Y_t

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Discretise the state space
stvy = sqrt(sigma^2/(1-lambda^2)); % standard deviation of y_t

ymax = m*stvy;                     % upper boundary of state space
ymin = -ymax;                      % lower boundary of state space

w = (ymax-ymin)/(N-1);             % length of interval 

s = ymin:w:ymax;                   % the discretised state space


% Calculate the transition matrix - see Tauchen
for j=1:N;
   for k=2:N-1;
      Tran(j,k)= normcdf(s(k)-lambda*s(j)+w/2,0,sigma)...
         - normcdf(s(k)-lambda*s(j)-w/2,0,sigma);
   end
    % Fill in end points
   Tran(j,1) = normcdf(s(1)-lambda*s(j)+w/2,0,sigma);
   Tran(j,N) = 1 - normcdf(s(N)-lambda*s(j)-w/2,0,sigma);
end

if sum(Tran') ~= ones(1,N)
   str = find(Tran'-ones(1,N));  % find rows not adding up to one
   disp('error in transition matrix');
   disp(['rows ',num2str(str),' does not sum to one']);
end


% Calculate the invariant distribution of Markov chain by simulating the
% chain to reach a long-run level
Trans= Tran';
probst = (1/N)*ones(N,1); % initial distribution of states
test = 1;

while test > 10^(-8);
	probst1 = Trans*probst;
	test=max(abs(probst1-probst));
	probst = probst1;   
end
   
% Calculate Properties of Invariant Distribution
meanm = s*probst;               % mean of invariant distribution
varm = ((s-meanm).^2)*probst;   % variance of invariant distribution
    
midaut1 = (s-meanm)'*(s-meanm); % cross product of deviation from the
                                % mean of y_t and y_t-1
                                   
probmat = probst*ones(1,N);     % each column is invariant distribution   
   
midaut2 = Tran.*probmat.*midaut1; % product of the first two terms is 
                                  % the joint distribution of (Y_t-1,Y_t)
                                                                      
autcov1 = sum(sum(midaut2));      % first-order auto-covariance
   
% Calculate the Asymptotic second moments of Markov chain
   
alambda = autcov1/varm;           % persistence of discrete process
disp(['Persistence of original process: ',num2str(lambda)])
disp(['Persistence of discrete Markov chain: ',num2str(alambda)]) 
disp('')
   
asigmay = sqrt(varm);            % s.d. of discrete process
disp(['S.D. of original process: ',num2str(stvy)])
disp(['S.D. of discrete Markov chain: ',num2str(asigmay)])