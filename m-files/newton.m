function x = newton(func,x0,param,crit,maxit);

% Newton.m Program to solve a system of equations
% x=newton(func, x0, param, crit, maxit) 

for i=1:maxit
  	[f,J]=feval(func,x0,param);
  	x=x0-inv(J)*f;
  	if norm(x-x0)<crit; 
    	   break 
  	end
  	x0=x;
end

if i>=maxit 
   sprintf('WARNING: Maximum number of %g iterations reached',maxit)
end

