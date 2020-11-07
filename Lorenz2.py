import numpy as np
import time as time
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.interpolate import interp1d
import audio as ad
import odestep as step
import multiprocessing as mp
from itertools import repeat
# Constants
sigma = 10.0
beta = 8/3
tlen = 100

#If speed = True --> rk45 integrator else odeint integrator
speed = True


def solve(ic,rho,t):
    sol = odeint(lorenz, ic, t, args=(sigma, rho, beta), rtol=1e-6, atol=1e-6)
    return sol

def solved(ic,rho,driving,t):
    finterp = interp1d(t,driving,fill_value='extrapolate')
    sol = odeint(lorenznew, ic, t, args=( sigma, rho, beta,driving,finterp), rtol=1e-6, atol=1e-6)
    return sol

def lorenz(q, t, sigma, rho, beta):
    x, y, z = q
    return [sigma*(y - x), x*(float(rho) - z) - y, x*y - beta*z]

def lorenznew(q, t, sigma, rho, beta, driving,finterp):
    x, y, z = q
    d = finterp(t)
    return [sigma*(y - x), d*(float(rho) - z) - y, d*y - beta*z]

def lorenzi(x,y,dx,**kwargs):
    dydx    = np.zeros(3)
    r = float(kwargs['rho'])
    # ????? from here
    dydx[0] = sigma*(y[1]-y[0])
    dydx[1] = r*y[0]-y[1]-y[0]*y[2]
    dydx[2] = y[0]*y[1]-beta*y[2]
    # ????? to here
    return dydx

def lorenzr(x,y,dx,**kwargs):
    dydx    = np.zeros(3)
 
    d = kwargs['d']
    smooth = kwargs['smooth']
    d = smooth(x)
    r = float(kwargs['rho'])



    # ????? from here
    dydx[0] = sigma*(y[1]-y[0])
    dydx[1] = r*d-y[1]-d*y[2]
    dydx[2] = d*y[1]-beta*y[2]
    # ????? to here
    return dydx
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
        finterp = interp1d(x,driving,fill_value='extrapolate')
        for k in range(1,nstep+1):
            y[:,k],it[k] = fORD(fRHS,x[k-1],y[:,k-1],dx,d = driving, rho = r, smooth =finterp)
    else:
        for k in range(1,nstep+1):        
            y[:,k],it[k] = fORD(fRHS,x[k-1],y[:,k-1],dx, rho = r)
    return x,y,it



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
        fRHS = lorenzr
    else:
        fRHS = lorenzi
    return fINT,fORD,fRHS,fBVP


def Unperturbed(sol1,rho, ic,t, stepper):
    tlen = t.size
    signal = np.zeros(sol1[:,0].size)
    driving = signal + sol1[:,0]

    if speed == False:
        tstart = time.time()

        fINT,fORD,fRHS,fBVP = ode_init(stepper, True)
        x,yrec,it = fINT(fRHS,fORD,fBVP,t[0],ic,t[t.size-1],signal.size-1,True, driving = driving, rho = rho)  

        tend = time.time()
        print(f"It took approximatly {tend-tstart} seconds")
    else:
        tstart = time.time()
        sol2 = solved(ic,rho, driving,t)
        yrec = np.transpose(sol2)
        tend = time.time()
        print(f"It took approximatly {tend-tstart} seconds")
    y = np.transpose(sol1)
    return y,yrec

def Problem2(n,rho):
    stepper = 'rk45'
   # ics = [rho**(1/2), rho**(1/2), rho+2]

    t = np.linspace(0, tlen, n)
    #ics = np.random.rand(3)
    ics = np.ones(3)*50
    if speed:
        ics = list(ics)

        sol1 = solve(ics,rho,t)

    else:
        ics = np.array(ics)
        fINT,fORD,fRHS,fBVP = ode_init(stepper,False)
        x,y,it = fINT(fRHS,fORD,fBVP,t[0],ics,t[t.size-1],n-1, False, rho = rho)  
        sol1 = np.transpose(y)
    ics += np.ones(3)*10
    y, yrec = Unperturbed(sol1,rho, ics,t, stepper)
    error = np.abs(y[0]-yrec[0])

    ss = y[0][n//2:-1].size
    m = np.min(np.abs(y[0][n//2:-1]-yrec[0][n//2:-1]))
    
    location = list(np.abs(y[0][n//2:-1]-yrec[0][n//2:-1])).index(m)
    location += ss

    m = np.average(error[location:-1])

    return n,m, location, error

def Problem3(rho,n,ic):
    stepper = 'rk45'
   # ics = [rho**(1/2), rho**(1/2), rho+2]

    t = np.linspace(0, tlen, n)
    #ics = np.random.rand(3)
    if speed:
        ics = list(ic)

        sol1 = solve(ics,rho,t)

    else:
        ics = np.array(ic)
        fINT,fORD,fRHS,fBVP = ode_init(stepper,False)
        x,y,it = fINT(fRHS,fORD,fBVP,t[0],ics,t[t.size-1],n-1, False, rho = rho)  
        sol1 = np.transpose(y)
    ics += np.ones(3)*10
    y, yrec = Unperturbed(sol1,rho, ics,t, stepper)
    error = np.abs(y[0]-yrec[0])

    m = np.average(error)

    return rho,m, error
if __name__ == "__main__":
    ndependance = True
    if ndependance:
        rho = 60
        n = [100,1000,5000,10000,50000,100000,1000000,5000000]#,9999999]

        print("multiprocessing:", end='')
        tstart = time.time()
        num_processes = 6 
        p = mp.Pool(num_processes)
        mp_solutions = p.starmap(Problem2, zip(n, repeat(rho)))
        #serial_solutions = [Problem2(ic) for ic in n]
        tend = time.time()
        tmp = tend - tstart
        print(" %8.3f seconds" % tmp)
        narray = []
        minarray = []
        for sol1 in mp_solutions:
            n, minim, location, error = sol1
            narray.append(n)
            minarray.append(minim)
            print(n,location/n*100, minim)
            if minim < 10**(-6):
                t = np.linspace(0, tlen, n)

                plt.plot(t[location:-1],error[location:-1],color='r')
                plt.title(f'{n} error')
                plt.show()
        plt.plot(np.log(narray),minarray)
        stri = f'Parameters:\n Sigma = {sigma}\n beta = {round(beta,3)}\n r = {rho}'
        plt.text(13, 13,stri, style='italic', bbox = {'facecolor': 'white'})
        plt.title('Performance of system for various Steps')
        plt.xlabel('log(Step)')
        plt.ylabel('Minimal Error')
        plt.show()
    else:
        ics = [[0.75,0.75,0.75],[10,10,10],[5,5,5],[100,100,100],[50,50,50]]
        n = 100000
        for ic in ics:
            print(f"IC in progress:{ic}")
            rho = np.linspace(25,60,15)

            ic2 = tuple(ic)
            print("multiprocessing:", end='')
            tstart = time.time()
            num_processes = 6 
            p = mp.Pool(num_processes)
            mp_solutions = p.starmap(Problem3, zip(rho, repeat(n), repeat(ic2)))
            #serial_solutions = [Problem2(ic) for ic in n]
            tend = time.time()
            tmp = tend - tstart
            print(" %8.3f seconds" % tmp)
            rhoarray = []
            minarray = []
            for sol1 in mp_solutions:
                rho, minim, error = sol1
                rhoarray.append(rho)
                minarray.append(minim)

            plt.plot(rhoarray,minarray, label=f'IC = {ic}')
        plt.legend()
        stri = f'Parameters:\n Sigma = {sigma}\n beta = {round(beta,3)}'
        plt.text(25, 0.028,stri, style='italic', bbox = {'facecolor': 'white'})
        plt.title('Performance of system for various R')
        plt.xlabel('rho')
        plt.ylabel('Drive Function Average Error')
        plt.show()