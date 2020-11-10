#Setup: conda install pyaudio

'''
What does this Function do?: The lorenz3.py primarily functions as a tool to 
encrypt messages whether they be binary or recorded. It not only runs a perturbed
synchronization to recover your message, but also runs an unperturbed synchronizaiton
to see if the values you chose lead to convergence.

Arguments:
    Stepper: What integrator is used?
        odeint: scipy integrator (100x faster)
        rk45: our adaptive integrator (100x slower but provides 100x more accurate results)
        rk4: Shouldn't be used if your looking for speed use odeint instead
        euler: included for completion sake
        rk2: included for completion sake
    Folder Name:
        This info is the name of the folder in which your sound files and data files are stored
        Note the folder is stored in the Data and Sound Folder respoectively
    r: This is the value for your r lorenz parameter
    sound:
        Record: Recording 
        Binary: Binary Waveform
    time:
        The time length of your recording or the time of the binary signal


Output:
The .wav files and waveform pictures are in the Sound folder

'''
from __future__ import division, print_function
import os 
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

beta = 8/3
sigma = 10

'''
IMPORTANT IF YOU HAVE A MAC OR UNIX CHANGE ALL THE \\ to / to properly save the file
ALSO INSTEAD of
    import audio as ad
USE
    import audiolinux as ad
'''

#Gets derivative setup for scipy integrator
def solve(ic,rho,t):
    sol = odeint(lorenz, ic, t, args=(sigma, rho, beta), rtol=1e-6, atol=1e-6)
    return sol

def solvrec(ic,rho,driving,t):
    finterp = interp1d(t,driving,fill_value='extrapolate')
    sol = odeint(lorenzrec, ic, t, args=(sigma, rho, beta,driving,finterp), rtol=1e-6, atol=1e-6)
    return sol
#Scipy integrator setup
def lorenz(q, t, sigma, rho, beta):
    x, y, z = q
    return [sigma*(y - x), x*(float(rho) - z) - y, x*y - beta*z]
def lorenzrec(q, t, sigma, rho, beta, driving,finterp):
    x, y, z = q
    d = finterp(t)
    return [sigma*(y - x), d*(float(rho) - z) - y, d*y - beta*z]



#Derivative evaluation for our integrator: sender
def ourlorenz(x,y,dx,**kwargs):
    dydx    = np.zeros(3)
    r = float(kwargs['rho'])
    # ????? from here
    dydx[0] = sigma*(y[1]-y[0])
    dydx[1] = r*y[0]-y[1]-y[0]*y[2]
    dydx[2] = y[0]*y[1]-beta*y[2]
    # ????? to here
    return dydx
#Derivative evaluation for: Receiver
def lorenzr(x,y,dx,**kwargs):
    dydx    = np.zeros(3)
 
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
        finterp = interp1d(x,driving,fill_value='extrapolate') #Uses basic linear extrapolation to get values within two support poitns
        for k in range(1,nstep+1):
            y[:,k],it[k] = fORD(fRHS,x[k-1],y[:,k-1],dx, rho = r, smooth =finterp)
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
        fRHS = ourlorenz
    return fINT,fORD,fRHS,fBVP

#Perturbed Synchronization 
def Perturbed(sol1, rho, ic, binary,t, stepper,timer,speed,Name,times, random):
    sample_rate = sol1[:,0].size/timer
    frac = 0.01*(t[-1]/100) #Higher t's require better fraction for better synchronization

    #Generates binary signal according to frequency and sample rate
    if binary:
        mode = 1
        cycpersec = 2 #Cycles/second of binary wave

        s = ad.Setter(Name, mode,T = cycpersec,time = times, A = 100,rate = sample_rate, rand = random) #Generates wave of length specified in optional args
        sig = np.zeros(sol1[:,0].size) #Makes signal length consistant with chaos length
        sig[0:s[1].size] = s[1]#Generates wave of length specified in optional args

    else:
        mode = 2
        s = ad.Setter(Name, mode,T = times, rate = sample_rate)
        sig = np.zeros(sol1[:,0].size)
        sig[0:s[1].size] = s[1]
        ad.play(s[1], Name, 'initial', sample_rate, time)


    #Plots Singla, chaos and mask
    signal = sig
    driving = frac*signal/max(signal) + sol1[:,0]

    plt.figure(figsize=(30,4))
    plt.plot(s[0], signal[0:s[1].size])
    plt.title("Signal")
    plt.xlabel("time")
    plt.show()
    plt.figure(figsize=(30,4))

    plt.plot(s[0],driving[0:s[1].size])
    plt.title("Mask")
    plt.xlabel("time")
    plt.show()

    #Fourier Analysis on the Sounds:
    if binary == False or binary == True:
        #Fourier Spectrum Plots
        figure, axes = plt.subplots(2, 1, figsize=(15,10))
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
        axes[0].legend()

        #2nd subplot graphs mask fft
        frequencies, fourierTransform = fft(driving[0:s[1].size], times, sample_rate) #Signal + Chaotic FFT
        axes[1].plot(frequencies[0:i], abs(fourierTransform[0:i]))
        axes[1].set_xlabel('Frequency')
        axes[1].set_ylabel('Amplitude')
        axes[1].set_title('Sig+Chaos')
        plt.legend()

    plt.show()        

    if speed == False: #Receiver differential equation recovery with our program integrator
        tstart = time.time()

        fINT,fORD,fRHS,fBVP = ode_init(stepper, True)
        x,y,it = fINT(fRHS,fORD,fBVP,t[0],ic,t[t.size-1],signal.size-1, True, driving = driving, rho = rho)  

        tend = time.time()
        print(f"Recovery: It took approximatly {tend-tstart} seconds")
    else: #Receiver differntial equation recovery with scipy integrator
        tstart = time.time()
        sol2 = solvrec(ic,rho, driving,t)
        y = np.transpose(sol2)
        tend = time.time()
        print(f"Recovery: It took approximatly {tend-tstart} seconds")


    #Gets and plots recovered signal
    recovered = (driving - y[0])*max(signal)/frac
    plt.figure(figsize=(30, 4))
    plt.plot(t[0:s[1].size], recovered[0:s[1].size])
    plt.title("Recovered Signal")

    plt.show()

    #Plays recovered and masked signal
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
            ad.play(recovered/2.0, Name, 'recovered', sample_rate, time)

            print('Signal too strong but was able to solve by halfing intensity')
            pass

#plots the lyaponov function, derivative and drive system error as a visual way to 
#quantify convergence
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

###
# Unperturbed function
# Performance syncriozation with signal = 0 to check for convergence
def Unperturbed(sol1,rho, ic,t, stepper,speed):
    signal = np.zeros(sol1[:,0].size)
    driving = signal + sol1[:,0]

    if speed == False: #Receiver differential equations with our integrator
        tstart = time.time()

        fINT,fORD,fRHS,fBVP = ode_init(stepper, True)
        x,yrec,it = fINT(fRHS,fORD,fBVP,t[0],ic,t[t.size-1],signal.size-1,True, driving = driving, rho = rho)  

        tend = time.time()
        print(f"{rho}: It took approximatly {tend-tstart} seconds")
    else: #Receiver differential equation with scipy integrator
        tstart = time.time()
        sol2 = solvrec(ic,rho, driving,t)
        yrec = np.transpose(sol2)
        tend = time.time()
        print(f"{rho}: It took approximatly {tend-tstart} seconds")

    y = np.transpose(sol1) 
    error(t, y, yrec) #Lyaponov evalutions and plot

    return y,yrec

##
# Fourier Transform Function:
# uses numpy's fourier transform to take the array "filename" into fourier space
# Code Based of: https://pythontic.com/visualization/signals/fouriertransform_fft
#   
#  
def fft(filename, time, sample_rate):
    fftfile = filename
    fourierTransform = np.fft.fft(fftfile)/len(fftfile)
    fourierTransform = fourierTransform[range(int(len(fftfile)/2))]
    tpCount     = len(fftfile)   
    frequencies = np.arange(int(tpCount/2))/(tpCount/sample_rate)
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
    parser.add_argument("r",type=float,default=1,
                        help="r tested")
    parser.add_argument("sound",type=str,default=1,
                        help="  record: recording\n"
                          "  binary: binary")
    parser.add_argument("time",type=float,default=1,
                        help="Recording or Binary Time")    

    args   = parser.parse_args()
    stepper= args.stepper
    Name = args.namefolder
    rho = float(args.r)
    sound = args.sound
    times = args.time


    '''
    3 varables that need attention
    '''
    #Keep n at 1 mil any lower and you can't test t = 1000 as well as t = 100
    n = 1000000 #Number of Steps: Needs to be atleast above 1 mil for good recovery or can be 100,000 if t= 100
    tlen = 1000 #Lorenz Integration Time: Need atleast above 100 for decent masking and a certain ratio between n/tlen has to be maintained for convergence dt<0

    sample_rate = 200000 #Sample Rate: Keep at 200000 and only lower if n goes down

    random = True #Random binary True = Yes False = no




    timer = n/sample_rate #For sample rate and steps whats the maximum time of recording

    if times>timer:
        raise NameError(f'Choose a time less than {timer} seconds or lower sample rate of {sample_rate}')


    if sound == 'record':
        binary = False
    else:
        binary = True

    if stepper == 'odeint':
        speed = True
    else:
        speed = False


    pather = os.getcwd() + '\\Data'
    try:
        os.mkdir(pather)
    except:
        pass
    pather = os.getcwd() + '\\Sound'
    try:
        os.mkdir(pather)
    except:
        pass

    
    t = np.linspace(0, tlen, n)
    ics = np.array([50.0,50.0,50.0])


    tstart = time.time()
    print("Initialization Started") #Takes about 100-200 seconds
    if speed: #uses scipy integrator
        sol1 = solve(ics,rho,t)
    else: #Uses our implemented integrator
        fINT,fORD,fRHS,fBVP = ode_init(stepper,False)
        x,y,it = fINT(fRHS,fORD,fBVP,t[0],ics,t[t.size-1],n-1, False, rho = rho)  
        sol1 = np.transpose(y)
    tend = time.time()
    tmp = tend - tstart
    print(f"Initialization took approximatly {tmp} seconds")
    if speed == False:
        print(f"Synchronization Prediction: {4.5*tmp} seconds")


    #uses parallel processing to check for convergence and for recovery
    p2 = mp.Process(target=Perturbed, args=(sol1,rho,ics,binary,t, stepper,timer,speed,Name,times, random)) #Process for pertubration
    p1 = mp.Process(target=Unperturbed, args=(sol1,rho, ics,t, stepper,speed)) #Process for synchronization check
    p1.start()
    p2.start()

