import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.interpolate import interp1d
import odestep as step

# Constants
sigma = 10
b = 8/3
tmax = 100
t = np.arange(0.0, tmax, 0.1)
def lorenzsystem(t, X, dx, **kwargs):
    r = float(kwargs['rho'])
    dydx = np.zeros(6)
    x1, x2, x3, y1, y2, y3 = X
    dydx[0] = sigma * (x2 - x1)
    dydx[1] = -x1 * x3 + r*x1 - x2
    dydx[2] = x1 * x2 - b*x3
    if 'd' in kwargs:
        smooth = kwargs['smooth']
        d = smooth(t)
        x1+=d

    dydx[3] = sigma * (y2 - y1)
    dydx[4] = -x1 * y3 + r*x1 - y2
    dydx[5] = x1 * y2 - b*y3
    return dydx
    
def ode_init(stepper, driver):
    if (stepper == 'rk4'):
        fORD = step.rk4
    elif (stepper == 'rk45'):
        fORD = step.rk45
    else:
        raise Exception('[ode_init]: invalid stepper value: %s' % (stepper))
    fBVP = 0 # default is IVP, but see below.
    fINT = ode_ivp
    if driver:
        fRHS = lorenzsystem
    else:
        fRHS = lorenzsystem
    return fINT,fORD,fRHS,fBVP

def ode_ivp(fRHS,fORD,fBVP,x0,y0,x1,nstep, driv, **kwargs):
    for key in kwargs:
        if (key=='driving'):
            driv = True
            driving = kwargs['driving']
    r = float(kwargs['rho'])

    nvar    = y0.size                      # number of ODEs
    x       = np.linspace(x0,x1,nstep+1)   # generates equal-distant support points
    y       = np.zeros((nvar,nstep+1))     # result array 
    y[:,0]  = y0                           # set initial condition
    dx      = x[1]-x[0]                    # step size
    it      = np.zeros(nstep+1)
    if driv:
        for k in range(1,nstep+1):
            finterp = interp1d(x,driving,fill_value='extrapolate')
            y[:,k],it[k] = fORD(fRHS,x[k-1],y[:,k-1],dx,d = driving, rho = r, smooth =finterp)
    else:
        for k in range(1,nstep+1):        
            y[:,k],it[k] = fORD(fRHS,x[k-1],y[:,k-1],dx, rho = r)
    return x,y,it

if __name__ == "__main__":
    ics = np.array([2000, 20, 30, 15, 20, 30])
    stepper = 'rk45'
    n = 10**4
    rho = 40

    fINT,fORD,fRHS,fBVP = ode_init(stepper,False)
    x,y,it = fINT(fRHS,fORD,fBVP,t[0],ics,t[t.size-1],n-1, False, rho = rho)

    driving = np.zeros(n)
    x,y,it = fINT(fRHS,fORD,fBVP,t[0],ics,t[t.size-1],n-1, True, driving = driving, rho = rho)  

    plt.plot(x,y[2]-y[-1])
    plt.show()