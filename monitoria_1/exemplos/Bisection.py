# Bisection method
def bisection(f, a, b, tol=1e-5):
    lower, upper = a, b 
    if f(a)*f(b)>=0:
        print("Method won't work. Needs opposite signs.")
    else:
        middle = 0.5*(upper + lower)
        while abs(f(middle)) > tol:
            if f(middle)*f(upper)<0:
                lower, upper = middle, upper
            else:
                lower, upper = lower, middle
            middle = 0.5*(lower+upper)
        return middle
            
