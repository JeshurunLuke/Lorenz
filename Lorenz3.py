from __future__ import division, print_function

import sys
import time
import multiprocessing as mp
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.interpolate import interp1d
import multiprocessing as mp
import audio as ad 
import argparse  # allows us to deal with arguments to main()

from argparse import RawTextHelpFormatter
import odestep as step
'''
Purpose of Lorenz3: It can perturb a signal with either a binary or recording and additionally checks if
for the value you chose if there is sychrnozation
Work needed: Implementation of Basic rk45 for 1st prob
             rk45 interpolate
'''
'''
argument: synchronization and masking depends on dt. You need a reasonable dt such that each iteration from the next has a noticable difference in the value of the 
evaluated function. Therefore decreasing dt to much will result in a slow "continous" function that fails to properly mask as efficently as a more "wild" function.
Nonetheless, this comes at a cost, because of this larger dt convergence and syncronization is not as efficent.

Conclusion: You can't win with both an amazingly masked function and an extremely pure recovered signal
Ways to optimize: Because of this conclusion our nstep/tlen has to be a constant, so the best plan of action is to optimize the magnitude of the hidden signal, such that 
the level of synchronization (for example with 0 signal and nste/tlen ratio given here syncrhonizes to (+/-0.005)), is a small fraction of the hidden signal.

Sample Rate max it because chaotic signal has to best frequency?
'''
sigma = 10.0
beta = 8/3


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

def solve(ic,rho,t):
    sol = odeint(lorenz, ic, t, args=(sigma, rho, beta), rtol=1e-6, atol=1e-6)
    return sol

def solved(ic,rho,driving,t):
    finterp = interp1d(t,driving,fill_value='extrapolate')
    sol = odeint(lorenznew, ic, t, args=(sigma, rho, beta,driving,finterp), rtol=1e-6, atol=1e-6)
    return sol



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


def Perturbed(sol1, rho, ic, binary,t, stepper,timer,speed,Name,times):
    tlen = t.size
    sample_rate = sol1[:,0].size/timer


    #Generates binary signal according to frequency and sample rate
    if binary:
        mode = 1
        cycpersec = 4 #Cycles/second

        s = ad.Setter(Name, mode,T = cycpersec,time = times, A = 30000,rate = sample_rate)
        sig = np.zeros(sol1[:,0].size)
        sig[0:s[1].size] = s[1]

    else:
        mode = 2
        s = ad.Setter(Name, mode,T = times, rate = sample_rate)
        sig = np.zeros(sol1[:,0].size)
        sig[0:s[1].size] = s[1]

    #Plots Singla, chaos and mask
    signal = sig
    driving = signal/max(signal) + sol1[:,0]

    figure, axs = plt.subplots(3, 1, figsize=(30,4))
    plt.subplots_adjust(hspace=1)
    axs[0].plot(s[0], signal[0:s[1].size])
    axs[0].set_title("Signal")
    axs[0].set_xlabel("time")
    axs[1].plot(s[0], sol1[0:s[1].size,0])
    axs[1].set_title("Chaos")
    axs[1].set_xlabel("time")
    axs[2].plot(s[0],driving[0:s[1].size])
    axs[2].set_title("Mask")
    axs[2].set_xlabel("time")
    plt.show()

    #Fourier Analysis on the Sounds:
    if binary == False:
        #Fourier Spectrum Plots
        figure, axes = plt.subplots(2, 1, figsize=(30,4))
        plt.subplots_adjust(hspace=1)

        frequencies, fourierTransform = fft(signal[0:s[1].size]/max(signal[0:s[1].size]), times, sample_rate) #Signal FFT
        #First subplot deals with graphing signal and chaotic fft
        i = list(frequencies).index(float(3000))
        axes[0].plot(frequencies[0:i], abs(fourierTransform[0:i]), label = 'Signal')

        frequencies, fourierTransform = fft(sol1[0:s[1].size,0]/max(sol1[0:s[1].size,0]), times, sample_rate) #Chaotic FFT
        axes[0].plot(frequencies[0:i], abs(fourierTransform[0:i]), label = 'Chaos')                       
        axes[0].set_xlabel('Frequency')
        axes[0].set_ylabel('Amplitude')
        axes[0].set_title('Signal')

        #2nd subplot graphs mask fft
        frequencies, fourierTransform = fft(driving[0:s[1].size], times, sample_rate) #Signal + Chaotic FFT
        axes[1].plot(frequencies[0:i], abs(fourierTransform[0:i]))
        axes[1].set_xlabel('Frequency')
        axes[1].set_ylabel('Amplitude')
        axes[1].set_title('Sig+Chaos')
        plt.legend()

    plt.show()        

    if speed == False:
        tstart = time.time()

        fINT,fORD,fRHS,fBVP = ode_init(stepper, True)
        x,y,it = fINT(fRHS,fORD,fBVP,t[0],ic,t[t.size-1],signal.size-1, True, driving = driving, rho = rho)  

        tend = time.time()
        print(f"It took approximatly {tend-tstart} seconds")
    else:
        tstart = time.time()
        sol2 = solved(ic,rho, driving,t)
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

def error(x, y, yrec):
    e = np.zeros((3,y[0].size))
    e[0] = y[0]-yrec[0]
    e[1] = y[1]-yrec[1]
    e[2] = y[2]-yrec[2]
    
    figure, axs = plt.subplots(3, 1, figsize=(30,4))

    plt.subplots_adjust(hspace=1)

    Ldot = -(e[0]-0.5*e[1])**2-0.75*(e[1]**2)-4*beta*(e[2]**2) #Lyaponev Derivative
    axs[0].plot(x, Ldot)
    axs[0].set_title('Lyaponev Derivative')
    L = 0.5*(1/sigma)*(e[0]**2) + e[1]**2 + e[2]**2 #Lyaponev Function
    axs[1].plot(x,L)
    axs[1].set_title('Lyaponev Function')
    axs[2].plot(x,e[0])
    axs[2].set_title('Error of Drive Function')

    plt.show()
    return L

def Unperturbed(sol1,rho, ic,t, stepper,timer,speed):
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
    error(t, y, yrec)

def fft(filename, time, sample_rate):
    fftfile = filename
    times  = np.arange(0, time, 1/sample_rate)
    fourierTransform = np.fft.fft(fftfile)/len(fftfile)
    fourierTransform = fourierTransform[range(int(len(fftfile)/2))]
    tpCount     = len(fftfile)
    values      = np.arange(int(tpCount/2))
    timePeriod  = tpCount/sample_rate
    frequencies = values/timePeriod
    return frequencies, fourierTransform


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument("stepper",type=str,default='euler',
                        help="stepping function:\n"
                             "   odeint: Production grade integrator\n"
                             "   rk4  : Runge-Kutta 2nd order\n"
                             "   rk45  : Runge-Kutta 4th order\n")

    parser.add_argument("namefolder",type=str,default=1,
                        help="paramaters tested")
    parser.add_argument("r",type=str,default=1,
                        help="r tested")
    parser.add_argument("sound",type=str,default=1,
                        help="  record: recording\n"
                          "  binary: binary")
    parser.add_argument("time",type=int,default=1,
                        help="Recording or Binary Time")    

    args   = parser.parse_args()
    stepper= args.stepper
    Name = args.namefolder
    global rho
    rho = args.r
    sound = args.sound
    times = args.time


    '''
    3 varables that need attention
    '''
    n = 100000 #Needs to be atleast above 100,000 for decent recovery
    tlen = 1000 #Need atleast above 100 for decent masking and a certain ratio between n/tlen has to be maintained for convergence dt<0
    #But if you try tlen = 10 with binary it breaks horribly for r>24.8 why?

    sample_rate = 20000 #How does sample rate influence?


    timer = n/sample_rate #Selects t such that sample rate is maxed to 200,000

    if times>timer:
        raise NameError(f'Choose a time less than {timer} seconds')
    if sound == 'record':
        binary = False
    else:
        binary = True

    if stepper == 'odeint':
        speed = True
    else:
        speed = False



    
    t = np.linspace(0, tlen, n)
    ics = np.random.rand(3)


    tstart = time.time()
    if speed:
        sol1 = solve(ics,rho,t)
    else:
        fINT,fORD,fRHS,fBVP = ode_init(stepper,False)
        x,y,it = fINT(fRHS,fORD,fBVP,t[0],ics,t[t.size-1],n-1, False, rho = rho)  
        sol1 = np.transpose(y)
    tend = time.time()
    tmp = tend - tstart
    print(f"Initialization took approximatly {tend-tstart} seconds")


    p2 = mp.Process(target=Perturbed, args=(sol1,rho,ics,binary,t, stepper,timer,speed,Name,times)) #Process for pertubration
    p1 = mp.Process(target=Unperturbed, args=(sol1,rho, ics,t, stepper,timer,speed)) #Process for synchronization check
    p1.start()
    p2.start()

