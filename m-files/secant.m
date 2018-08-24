function x=secant(func, x0, param, crit, maxit)
%   secant.m solves a system of equations
%   f(z1, z2,...,zn)=0 with the secant method
%   where x=[z1, z2,...,zn] is the solution vector.
%   function 'func', which is a string.
%   param corresponds to additonal parameters of the function 'func'.
%   x0, crit, and maxit
del=diag(max(abs(x0)*1e-4, 1e-8));
n=length(x0);
    for i=1:maxit
        f=feval(func,x0,param);
        for j=1:n
            J(:,j)=(f-feval(func,x0-del(:,j),param))/del(j,j);
        end
        x=x0-inv(J)*f;
        if norm(x-x0)<crit; break;  end
        x0=x;
    end
    if i>=maxit
        sprintf('maximum number of iterations was reached')
    end
   
        