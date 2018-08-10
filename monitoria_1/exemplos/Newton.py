#Univariate Newton method
def newton_univariate(f, df, start, tol = 1e-5):
    iteration = start
    while abs(f(iteration)) > tol:
        iteration = iteration - f(iteration)/df(iteration)     
    return iteration

        