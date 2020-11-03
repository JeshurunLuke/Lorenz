from __future__ import division, print_function

import sys
import time
import multiprocessing as mp
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import audio as ad 


n = 1000000
tlen = 1000
rho = 40.0
timer = 5


binary = False

Name  = 'Truest'


def param(x,**kwargs):
    sigma = 10.0
    beta = 8/3

    # ????? to here
    return sigma,beta,rho
def lorenz(q, t, sigma, rho, beta):
    x, y, z = q
    return [sigma*(y - x), x*(rho - z) - y, x*y - beta*z]

def lorenzr(x,y,dx,**kwargs):
    dydx    = np.zeros(3)
    sigma, b, r = param(x,**kwargs)
    d = float(kwargs['d'])

    # ????? from here
    dydx[0] = sigma*(y[1]-y[0])
    dydx[1] = r*d-y[1]-d*y[2]
    dydx[2] = d*y[1]-b*y[2]
    # ????? to here
    return dydx

def solve(ic):
    t = np.linspace(0, tlen, n)
    sigma = 10.0
    beta = 8/3
    sol = odeint(lorenz, ic, t, args=(sigma, rho, beta), rtol=1e-10, atol=1e-12)
    return sol

def solved(ic,driving):
    sigma = 10.0
    rho = 28.0
    beta = 8/3
    set_step()
    sol = odeint(lorenzr, ic, t, args=(sigma, rho, beta,driving), rtol=1e-10, atol=1e-12)
    return sol

def rk4(fRHS,x0,y0,dx,**kwargs):
    #???????? from here
    k1 = fRHS(x0,y0,dx,**kwargs)*dx
    k2 = fRHS(x0+dx/2,y0+k1/2,dx,**kwargs)*dx
    k3 = fRHS(x0+dx/2,y0+k2/2,dx,**kwargs)*dx
    k4 = fRHS(x0+dx,y0+k3,dx,**kwargs)*dx
    y = y0 + (k1+k4)/6 + (k2+k3)/3
    #???????? to here
    return y,1

def ode_ivp(fRHS,fORD,fBVP,x0,y0,x1,nstep, **kwargs):
    driv = False
    for key in kwargs:
        if (key=='driving'):
            driv = True
            driving = kwargs['driving']                      

    nvar    = y0.size                      # number of ODEs
    x       = np.linspace(x0,x1,nstep+1)   # generates equal-distant support points
    y       = np.zeros((nvar,nstep+1))     # result array 
    y[:,0]  = y0                           # set initial condition
    dx      = x[1]-x[0]                    # step size
    it      = np.zeros(nstep+1)
    if driv:
        for k in range(1,nstep+1):
        
            currentx = driving[k]
            y[:,k],it[k] = fORD(fRHS,x[k-1],y[:,k-1],dx,d = currentx)
    else:
        for k in range(1,nstep+1):        
            currentx = driving[k]
            y[:,k],it[k] = fORD(fRHS,x[k-1],y[:,k-1],dx,d = currentx)
    return x,y,it



def ode_init():

    fBVP = 0 # default is IVP, but see below.
    fINT = ode_ivp
    fORD = rk4
    fRHS = lorenzr
    return fINT,fORD,fRHS,fBVP



if __name__ == "__main__":
    t = np.linspace(0, tlen, n)

    ics = np.random.randn(1, 3)
    print(ics)

    print("multiprocessing:", end='')
    tstart = time.time()
    num_processes = 5
    p = mp.Pool(num_processes)
    mp_solutions = p.map(solve, ics)
    tend = time.time()
    tmp = tend - tstart
    print(" %8.3f seconds" % tmp)

    n = 0 
    for sol1 in mp_solutions:
        sample_rate = sol1[:,0].size/timer
        cycpersec = 4 #Cycles/second


        #Generates binary signal according to frequency and sample rate
        if binary:
            mode = 1

            s = ad.Setter(Name, mode,T = cycpersec,time = timer, A = 30000,rate = sample_rate)
        else:
            mode = 2
            s = ad.Setter(Name, mode,T = timer, rate = sample_rate)
            plt.figure(figsize=(30, 4))
            plt.plot(t, s[1])
            plt.title("Initial Signal")
            plt.draw()

        signal = s[1] #np.zeros(sol1[:,0].size)
        driving = signal/max(signal) + sol1[:,0]



        print(driving)


        fINT,fORD,fRHS,fBVP = ode_init()
        x,y,it = fINT(fRHS,fORD,fBVP,t[0],ics[n],t[t.size-1],signal.size-1, driving = driving)        
        recovered = (driving - y[0])*max(signal)
        if binary == False:
            print("masked signal")
            try:
                ad.play(driving*500, Name, 'masked', sample_rate, time)
            except:
                print('Signal too strong')
                pass
            print("recovered signal")
            try:
                ad.play(recovered, Name, 'recovered', sample_rate, time)
            except:
                print('Signal too strong')
                pass


        plt.figure(figsize=(30, 4))
        plt.plot(t, recovered)
        plt.title("Recovered Signal")
        plt.show()

        n+=1
        

'''
    print("serial:         ", end='')
    sys.stdout.flush()
    tstart = time.time()
    serial_solutions = [solve(ic) for ic in ics]
    tend = time.time()
    tserial = tend - tstart
    print(" %8.3f seconds" % tserial)

    print("num_processes = %i, speedup = %.2f" % (num_processes, tserial/tmp))

    check = [(sol1 == sol2).all()
             for sol1, sol2 in zip(serial_solutions, mp_solutions)]
    if not all(check):
        print("There was at least one discrepancy in the solutions.")
'''