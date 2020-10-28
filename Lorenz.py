import numpy as np
import odestep as step
import utilities as util
import argparse	                 # allows us to deal with arguments to main()
from argparse import RawTextHelpFormatter
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from scipy.fftpack import fft
import os
import multiprocessing as mp
import audio as ad



def param(x,**kwargs):
    for key in kwargs:
        if (key=='r'):
            r = float(kwargs[key])
        if (key=='s'):
            sigma = float(kwargs[key])
        if (key=='b'):
            b = float(kwargs[key])
    sigma = 10

    # ????? to here
    return sigma,b,r

#========================================
# fRHS for the (non-)uniform string
# Should return an array dydx containing three
# elements corresponding to y'(x), y''(x), and lambda'(x).
def dydx_lorenz(x,y,dx,**kwargs):
    dydx    = np.zeros(3)
    sigma, b, r = param(x,**kwargs)
    # ????? from here
    dydx[0] = sigma*(y[1]-y[0])
    dydx[1] = r*y[0]-y[1]-y[0]*y[2]
    dydx[2] = y[0]*y[1]-b*y[2]
    # ????? to here
    return dydx
def dydx_rlorenz(x,y,dx, **kwargs):
    for key in kwargs:
        if key == 'driving':
            driving = kwargs[key]

    dydx    = np.zeros(3)
    sigma, b, r = param(x,**kwargs)

    currentx = driving[steps]
  

    # ????? from here
    dydx[0] = sigma*(y[1]-y[0])
    dydx[1] = r*currentx-y[1]-currentx*y[2]
    dydx[2] = currentx*y[1]-b*y[2]
    return dydx
def dydx_rylorenz(x,y,dx, **kwargs):
    for key in kwargs:
        if key == 'driving':
            driving = kwargs[key]

    dydx    = np.zeros(3)
    sigma, b, r = param(x,**kwargs)

    currenty = driving[steps]
  

    # ????? from here
    dydx[0] = sigma*(currenty-y[0])
    dydx[1] = r*y[0]-y[1]-y[0]*y[2]
    dydx[2] = y[0]*currenty-b*y[2]
    return dydx

    # ????? to her
def get_step():
    global steps
    steps +=1
    return steps



def ode_init(stepper,nstep, **kwargs):
    
    ver = int(kwargs['version'])
    time = int(kwargs['time'])
    var = int(kwargs['var'])

    fBVP = 0 # default is IVP, but see below.
    if (stepper == 'euler'):
        fORD = step.euler
    elif (stepper == 'rk2'):
        fORD = step.rk2
    elif (stepper == 'rk4'):
        fORD = step.rk4
    elif (stepper == 'rk45'):
        fORD = step.rk45
    else:
        raise Exception('[ode_init]: invalid stepper value: %s' % (stepper))
    x1 = time
    x0 = 0
    fINT = ode_ivp
    if ver == 1:
        y0 = np.array([10,10,10])
        fRHS = dydx_lorenz
    elif ver ==2:
        y0 = np.array([10,10,10]) #+ np.array([0.1,0.1,0.1]) #Displacement Vector
        if var == 1:
            fRHS = dydx_rylorenz
        elif var == 0:
            fRHS = dydx_rlorenz
    fBVP = 0
    return fINT,fORD,fRHS,fBVP,x0,y0,x1,nstep

#========================================
# Single rk4 step.
# Already provided. Good to go; 2.5 free points :-)

#========================================
# ODE IVP driver.
# Already provided. Good to go; 2.5 free points :-)
def ode_ivp(fRHS,fORD,fBVP,x0,y0,x1,nstep,**kwargs):
    for key in kwargs:
        if (key=='driving'):
            get_step()  
                   

    nvar    = y0.size                      # number of ODEs
    x       = np.linspace(x0,x1,nstep+1)   # generates equal-distant support points
    y       = np.zeros((nvar,nstep+1))     # result array 
    y[:,0]  = y0                           # set initial condition
    dx      = x[1]-x[0]                    # step size
    it      = np.zeros(nstep+1)
    for k in range(1,nstep+1):
        for key in kwargs:
            if (key=='driving'):
                get_step()         
        y[:,k],it[k] = fORD(fRHS,x[k-1],y[:,k-1],dx,**kwargs)
    return x,y,it


#=======================================
# A single trial shot.
# Sets the initial values (guesses) via fLOA, calculates 
# the corresponding solution via ode_ivp, and returns 
# a "score" via fSCO, i.e. a value for the rootfinder to zero out.

#=======================================
def check(x,y,it,r,Name, **kwargs):
    col = ['black','green','cyan','blue','red','black','black','black','black']
    version = False
    for key in kwargs:
        if (key=='receiving'):
            version = True
            yrec = kwargs[key]
      
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    

    ax.plot(y[0, :], y[1, :], y[2, :], color=col[1])
    if version:
        ax.plot(yrec[0, :], yrec[1, :], yrec[2, :], color=col[4])
    # Setting the axes properties
    ax.set_xlim3d([min(y[0])-1, max(y[0])+1])
    ax.set_xlabel('X')

    ax.set_ylim3d([min(y[1])-1, max(y[1])+1])
    ax.set_ylabel('Y')

    ax.set_zlim3d([min(y[2])-1, max(y[2])+1])
    ax.set_zlabel('Z')
    ax.legend()



def update(num, data,line):
    line.set_data(data[:2, :num])
    line.set_3d_properties(data[2, :num])

def set_step():
    global steps
    steps = -1

def error(x, y, yrec, sigma,b):
    e = np.zeros((3,len(y[0])))
    e[0] = y[0]-yrec[0]
    e[1] = y[1]-yrec[1]
    e[2] = y[2]-yrec[2]
    
    figure, axs = plt.subplots(2, 1, figsize=(30,4))

    plt.subplots_adjust(hspace=1)

    Ldot = -(e[0]-0.5*e[1])**2-0.75*(e[1]**2)-4*b*(e[2]**2) #Lyaponev Derivative
    axs[0].plot(x, Ldot)
    axs[0].set_title('Lyaponev Derivative')
    L = 0.5*(1/sigma)*(e[0]**2) + e[1]**2 + e[2]**2 #Lyaponev Function
    axs[1].plot(x,L)
    axs[1].set_title('Lyaponev Function')
    plt.show()
    return L


def run(stepper, nstep, ver, Name, rparam, sound, drive):
    sample_rate = 50000

    if drive == 'x':
        g = 0
    elif drive == 'y':
        g = 1
    if ver ==2:
        time = (nstep+1)/sample_rate #NEEDS TO BE DIVISABLE
        
        fINT,fORD,fRHS,fBVP,x0,y0,x1,nstep = ode_init(stepper,nstep, version = 1, var = g, time = time)
        x,y,it = fINT(fRHS,fORD,fBVP,x0,y0,x1,nstep,s=10,b=8/3,r=rparam)

        set_step()

        fINT,fORD,fRHS,fBVP,x0,y0,x1,nstep = ode_init(stepper,nstep, version = ver, var = g, time = time)
        x,yrec,it2 = fINT(fRHS,fORD,fBVP,x0,y0,x1,nstep,s=10,b=8/3,r=rparam, driving = y[g])
        L = error(x,y,yrec,10,8/3)
        check(x,y,it, rparam,Name,  receiving = yrec)
    elif ver == 1:
        time = 100
        fINT,fORD,fRHS,fBVP,x0,y0,x1,nstep = ode_init(stepper,nstep, version = 1, var = 0, time = time)
        x,y,it = fINT(fRHS,fORD,fBVP,x0,y0,x1,nstep,s=10,b=8/3,r=rparam)
        check(x,y,it, rparam,Name)
    elif ver == 3:
        if sound == 'record':
            mode = 2
            time = (nstep+1)/sample_rate #NEEDS TO BE DIVISABLE
            print(f'Message Time: {time}')


            fINT,fORD,fRHS,fBVP,x0,y0,x1,nstep = ode_init(stepper,nstep, version = 1, var = g, time = time)
            x,y,it = fINT(fRHS,fORD,fBVP,x0,y0,x1,nstep,s=10,b=8/3,r=rparam)


            s = ad.Setter(Name, mode,T = time, rate = sample_rate)
       
            frac = 0.001
            mag = frac*max(abs(y[g]))/(max(abs(s[1])))
            mask = s[1]*mag + y[g]  
            nstep = len(s[1])-1
            
            #mag = max(abs(y[1]))/(frac*max(abs(y[0])))

            figure, axs = plt.subplots(3, 1, figsize=(30,4))
            plt.subplots_adjust(hspace=1)
            axs[0].plot(s[0], s[1]) 
            axs[0].set_xlim(s[0][0], s[0][-1])
            axs[0].set_xlabel('time (s)')
            axs[0].set_ylabel('amplitude')
            axs[0].set_title('Chaotic x')
            axs[1].plot(s[0], y[g]) 
            axs[1].set_xlim(s[0][0], s[0][-1])
            axs[1].set_xlabel('time (s)')
            axs[1].set_ylabel('amplitude')
            axs[1].set_title('Chaotic x')
            axs[2].plot(s[0], mask) 
            axs[2].set_xlim(s[0][0], s[0][-1])
            axs[2].set_xlabel('time (s)')
            axs[2].set_ylabel('amplitude')
            axs[2].set_title('Chaos + signal')
            #plt.show()
            print("masked signal")
            ad.play(mask*500, Name, 'masked', sample_rate, time)
            #Fourier Spectrium
            figure, axes = plt.subplots(2, 1, figsize=(30,4))

            plt.subplots_adjust(hspace=1)



            frequencies, fourierTransform = fft(s[1]/max(s[1]), time, sample_rate) #Signal FFT

            axes[0].plot(frequencies, abs(fourierTransform), label = 'Signal')
            frequencies, fourierTransform = fft(y[0]/max(y[0]), time, sample_rate) #Chaotic FFT
            axes[0].plot(frequencies, abs(fourierTransform), label = 'Chaos')


                       
            axes[0].set_xlabel('Frequency')
            axes[0].set_ylabel('Amplitude')
            axes[0].set_title('Signal')

            frequencies, fourierTransform = fft(mask, time, sample_rate) #Signal + Chaotic FFT


            axes[1].plot(frequencies, abs(fourierTransform))
            axes[1].set_xlabel('Frequency')
            axes[1].set_ylabel('Amplitude')
            axes[1].set_title('Sig+Chaos')
            plt.legend
            plt.show()

            set_step()

            fINT,fORD,fRHS,fBVP,x0,y0,x1,nstep = ode_init(stepper,nstep, version = 2, var = g, time = time)
            x,yrec,it2 = fINT(fRHS,fORD,fBVP,x0,y0,x1,nstep,s=10,b=8/3,r=rparam, driving = mask)
        elif sound == 'binary':
            mode = 1
            time = (nstep+1)/sample_rate
            cycpersec = 2 #Cycles/second
            s = ad.Setter(Name, mode,T = cycpersec,time = time, A = 30000,rate = sample_rate)
            time = 100
            nstep = len(s[1])-1

            fINT,fORD,fRHS,fBVP,x0,y0,x1,nstep = ode_init(stepper,nstep, version = 1, var = g, time = time)
            x,y,it = fINT(fRHS,fORD,fBVP,x0,y0,x1,nstep,s=10,b=8/3,r=rparam)

            frac = 0.1
            mag = frac*max(abs(y[g]))/(max(abs(s[1])))
            mask = s[1]*mag + y[g]
            figure, axs = plt.subplots(3, 1, figsize=(30,4))
            plt.subplots_adjust(hspace=1)
            axs[0].plot(s[0],s[1])
            axs[0].set_title("Signal")
            axs[0].set_xlabel("time")
            axs[1].plot(s[0],y[g])
            axs[1].set_title("Chaos")
            axs[1].set_xlabel("time")
            axs[2].plot(s[0],mask)
            axs[2].set_title("Mask")
            axs[2].set_xlabel("time")
            plt.show()
            fINT,fORD,fRHS,fBVP,x0,y0,x1,nstep = ode_init(stepper,nstep, version = 2, var = g, time = time)
            set_step()

            x,yrec,it2 = fINT(fRHS,fORD,fBVP,x0,y0,x1,nstep,s=10,b=8/3,r=rparam, driving = mask)
        
        recovered = (mask-yrec[g])/mag
        plt.figure(figsize=(30, 4))
        plt.plot(s[0], recovered)
        plt.show()

        if mode == 2:
            if abs(max(recovered))>32767:
                recovered = recovered*32767/(max(recovered+10))

            print("Recovered Signal")
            ad.play(recovered, Name, 'recovered', sample_rate, time )
        

    plt.show()
    

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

    #Default Run : python Lorenz.py rk4 100000 3 Practice 30 record y
    '''
    Work: what values of r best synchronize the message
          Is there a single value we can use to characterize synchronizaiton?
          How does step size effect the convergence of the synchronization?
          Is there any parallel work or anything else we can use to speed up the integrator?
    '''
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument("stepper",type=str,default='euler',
                        help="stepping function:\n"
                             "   euler: Euler step\n"
                             "   rk2  : Runge-Kutta 2nd order\n"
                             "   rk4  : Runge-Kutta 4th order\n")
    parser.add_argument("steps",type=int,default=100,
                        help="number of steps")
    parser.add_argument("ver",type=int,default=1,
                        help="  1:Lorenz\n" 
                        "  2:Unperturbed\n"
                        "  3: Record/Binary Perturbation")
    parser.add_argument("namefolder",type=str,default=1,
                        help="paramaters tested")
    parser.add_argument("r",type=str,default=1,
                        help="r tested")
    parser.add_argument("sound",type=str,default=1,
                        help="  record: recording\n"
                          "  binary: binary\n" 
                          "  NA: for version 1 and 2" )
    parser.add_argument("drive",type=str,default=1,
                        help="Drive signal?\n"
                          "  x: x is your drive\n" 
                          "  y: y is your drive\n"
                          "  NA: if version 1" )

    args   = parser.parse_args()

    stepper= args.stepper
    nstep = args.steps-1
    ver = args.ver
    Name = args.namefolder
    rparam = args.r
    sound = args.sound
    drive = args.drive


    pather = os.getcwd() + '\\Data' + "\\" + Name
    '''
    try:
        os.mkdir(pather)
    except:
        pass
    '''
    if ver == 3:
       # run(stepper, nstep, ver, Name, rparam, sound, drive)
        
        p2 = mp.Process(target=run, args=(stepper, nstep, ver, Name, rparam, sound, drive))
        p1 = mp.Process(target=run, args=(stepper, nstep, 2, Name, rparam, 'NA', drive))
        p1.start()
        p2.start()
    else: 
        run(stepper, nstep, ver, Name, rparam, sound, drive )    


#=======================================
