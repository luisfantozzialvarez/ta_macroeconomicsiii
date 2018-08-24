function S = markov(T,prob,S0);
%  FUNCTION MARKOV
%  
% 
%  This function creates an realization of a random variable
%  that comes from a markov chain.
%
%  MARKOV(T,PROB,S0) creates T realizations of a r.v. 
%  with transition matrix PROB. By default, the states are 1,2,...,m
%  where m is the size of PROB
%
%  MARKOV(T,PROB) creates T realizations by generating
%  randomly an initial state.
%
%  See also IID
%


% 1. The following take cares of possible mistakes in the inputs of the function

if nargin == 1; % detects the number of inputs, need 2
    error('The function needs at least two inputs') 
else
    [r1,c1] = size(T);
    [r2,c2] = size(prob);
    
    % (a) checking the first input
    if r1 == 1 & c1 == 1 % confirms that first input is a scalar
        % if the number of realizations is not a natural number, 
        % the following takes the abs. value and rounds to the nearest natural number.
        T = round(abs(T));
    else % checks if T is a scalar or not
        error('The first input must be a scalar') 
    end
            
    % (b) checking the second input
    if (r2 == 1 | c2 == 1) | (r2 ~= c2) % checks if the second input is a square matrix or not
        error('The second input must be a square matrix')
    end
    
    prob = abs(prob); % changes the sign if some probabilities are negative
    
    % checks if the probabilities sum to one, if not, it normalizes
    for i = 1:r2
        if sum(prob(i,:)) ~= 1
           warning(' The probabilities don`t sum to 1. Normalizing probabilities...')
           prob(i,:) = prob(i,:)/sum(prob(i,:));
        end
     end
    
    % (c) checking third input/defining initial state
    S    = zeros(T,1);
    if nargin == 3; 
       [r3,c3] = size(S0);
       if (r3 == 1 & c3 == 1) & (S0 >= 1 & S0 <= r2)
          S(1) = round(abs(S0));
       else
          error('The third argument must be a scalar between 1 and the size of the 2nd argument')
       end
    else 
       S(1) = ceil(rand*r2); % sets a default initial state if not determined by the user   
    end
    
end


% ----------------------------------------------------------------------------------
% 2. Creating the shock realizations

seed = sum(100*clock);
P = cumsum(prob');            % Creates a matrix with cumulated sum over rows, (prob') 

for t = 2:T;
  j      = find(rand < P(:,S(t-1)));
  S(t,1) = j(1);
end
