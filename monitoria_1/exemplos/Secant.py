import numpy as np 
#computes partial derivative of ith equation at jth variable 
def partial(f,x, i=0, j=0, prec= 1e-5):
    dev, x = np.array(x), np.array(x)
    dev[j] = dev[j]+ prec
    return (f(dev)[i] - f(x)[i])/prec

#computes Jacobian of Function
def jacobian(f, x):
    x= np.array(x)
    jacob = np.empty(shape= (len(f(x)), len(x)))
    for i in range(len(f(x))):
        for j in range(len(x)):
            jacob[i,j] = partial(f,x, i=i, j=j)
    return(jacob)

#secant Method: uses Jacobian approximation
def secant(f, start, tol=1e-5):
    iteration = np.array(start)
    while np.linalg.norm(f(iteration), ord=2)> tol:
        iteration = iteration - np.dot(np.linalg.inv(jacobian(f,iteration)),f(iteration))
    return iteration