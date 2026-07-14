from numpy import sinc, exp, sin, cos, pi

def sinc_fit(x, a,b,c,d):
    return a*sinc(2*c*(x-b) / pi)**2+d

def gaussian_fit(x,a,b,c):
    return a*exp(-1/(2*c**2)*(x-b)**2)

def exp_decay_fit(x, a,b,c):
    return a*exp(-x/b)+c

def parabola_fit(x, a,b,c):
    return a*(x-b)**2 + c

def exp_sin_fit(x, a,b,c, d, e):
    return a*exp(-b*x)*sin(c*x+d)+e

def cos_fit(x, a, b, c, d):
    return a*cos(b*x+c)+d