from __future__ import division, print_function

import sys
import time
import multiprocessing as mp
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.interpolate import interp1d

import audio as ad 
import odestep as step

'''
argument: synchronization and masking depends on dt. You need a reasonable dt such that each iteration from the next has a noticable difference in the value of the 
evaluated function. Therefore decreasing dt to much will result in a slow "continous" function that fails to properly mask as efficently as a more "wild" function.
Nonetheless, this comes at a cost, because of this larger dt convergence and syncronization is not as efficent.

Conclusion: You can't win with both an amazingly masked function and an extremely pure recovered signal
Ways to optimize: Because of this conclusion our nstep/tlen has to be a constant, so the best plan of action is to optimize the magnitude of the hidden signal, such that 
the level of synchronization (for example with 0 signal and nste/tlen ratio given here syncrhonizes to (+/-0.005)), is a small fraction of the hidden signal.

Sample Rate max it because chaotic signal has to best frequency?
'''
n = 100000
tlen = 10


rho = 40.0

timer = 5
times = 5

binary = True
speed = True

Name  = 'Truest'

tnew = np.linspace(0, tlen, n)

def param(x,**kwargs):
    sigma = 10.0
    beta = 8/3

    # ????? to here
    return sigma,beta,rho
def lorenz(q, t, sigma, rho, beta):
    x, y, z = q
    return [sigma*(y - x), x*(rho - z) - y, x*y - beta*z]

def lorenznew(q, t, sigma, rho, beta, driving,finterp):
    x, y, z = q
    d = finterp(t)
    return [sigma*(y - x), d*(rho - z) - y, d*y - beta*z]


def lorenzr(x,y,dx,**kwargs):
    dydx    = np.zeros(3)
    d = float(kwargs['d'])
    sigma = 10.0
    b = 8/3
    r = rho

    # ????? from here
    dydx[0] = sigma*(y[1]-y[0])
    dydx[1] = r*d-y[1]-d*y[2]
    dydx[2] = d*y[1]-b*y[2]
    # ????? to here
    return dydx

def solve(ic):
    t = np.linspace(0, tlen, n)
    print(t)
    sigma = 10.0
    beta = 8/3

    sol = odeint(lorenz, ic, t, args=(sigma, rho, beta), rtol=1e-10, atol=1e-12)
    return sol

def solved(ic,driving):
    t = np.linspace(0, tlen, n)

    sigma = 10.0
    beta = 8/3

    finterp = interp1d(t,driving,fill_value='extrapolate')
    sol = odeint(lorenznew, ic, t, args=(sigma, rho, beta,driving,finterp), rtol=1e-6, atol=1e-6)
    return sol



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
        
            currentx = driving[k-1]
            y[:,k],it[k] = fORD(fRHS,x[k-1],y[:,k-1],dx,d = currentx)
    else:
        for k in range(1,nstep+1):        
            y[:,k],it[k] = fORD(fRHS,x[k-1],y[:,k-1],dx)
    return x,y,it



def ode_init():

    fBVP = 0 # default is IVP, but see below.
    fINT = ode_ivp
    fORD = step.rk45
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

    c = 0 
    for sol1 in mp_solutions:
        sample_rate = sol1[:,0].size/timer
        print(sample_rate)
        cycpersec = 4 #Cycles/second


        #Generates binary signal according to frequency and sample rate
        if binary:
            mode = 1
            
            s = ad.Setter(Name, mode,T = cycpersec,time = timer, A = 30000,rate = sample_rate)
            sig = s[1]
        else:
            mode = 2
            s = ad.Setter(Name, mode,T = times, rate = sample_rate)
            sig = np.zeros(sol1[:,0].size)
            sig[0:s[1].size] = s[1]

            plt.figure(figsize=(30, 4))
            plt.plot(t[0:s[1].size], sig[0:s[1].size])
            plt.title("Initial Signal")
            plt.show()

        signal = sig
        driving = signal/max(signal) + sol1[:,0]


        if speed == False:
            tstart = time.time()

            fINT,fORD,fRHS,fBVP = ode_init()
            x,y,it = fINT(fRHS,fORD,fBVP,t[0],ics[c],t[t.size-1],signal.size-1, driving = driving)  

            tend = time.time()
            print(f"It took approximatly {tend-tstart} seconds")
        else:
            tstart = time.time()
            sol2 = solved(ics[c],driving)
            y = np.transpose(sol2)
            tend = time.time()
            print(f"It took approximatly {tend-tstart} seconds")



        recovered = (driving - y[0])*max(signal)
        plt.figure(figsize=(30, 4))
        plt.plot(t[0:s[1].size], recovered[0:s[1].size])
        plt.title("Recovered Signal")
        plt.show()

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



        c+=1
        

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
