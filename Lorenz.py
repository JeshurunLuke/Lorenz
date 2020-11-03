import argparse  # allows us to deal with arguments to main()
import multiprocessing as mp
import os
from argparse import RawTextHelpFormatter

import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import audio as ad
import odestep as step
import utilities as util


def param(x,**kwargs):
    for key in kwargs:
        if (key=='r'):
            r = float(kwargs[key])
        if (key=='s'):
            sigma = float(kwargs[key])
        if (key=='b'):
            b = float(kwargs[key])

    # ????? to here
    return sigma,b,r


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
    try:
        time = int(kwargs['time'])
    except:
        print("Defaulting to: 100 ")
        time = 100
    ver = int(kwargs['version'])
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
    x1 = time*10 #time = 1
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

def error(x, y, yrec, sigma,b,g):
    e = np.zeros((3,len(y[0])))
    e[0] = y[0]-yrec[0]
    e[1] = y[1]-yrec[1]
    e[2] = y[2]-yrec[2]
    
    figure, axs = plt.subplots(3, 1, figsize=(30,4))

    plt.subplots_adjust(hspace=1)

    Ldot = -(e[0]-0.5*e[1])**2-0.75*(e[1]**2)-4*b*(e[2]**2) #Lyaponev Derivative
    axs[0].plot(x, Ldot)
    axs[0].set_title('Lyaponev Derivative')
    L = 0.5*(1/sigma)*(e[0]**2) + e[1]**2 + e[2]**2 #Lyaponev Function
    axs[1].plot(x,L)
    axs[1].set_title('Lyaponev Function')
    axs[2].plot(x,e[g])
    axs[2].set_title('Error of Drive Function')

    plt.show()
    return L


def run(stepper, nstep, ver, Name, rparam, sound, drive):
    sample_rate = 100000

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

        L = error(x,y,yrec,10,8/3,g)
        check(x,y,it, rparam,Name,  receiving = yrec)
    elif ver == 1:
        time = 100
        fINT,fORD,fRHS,fBVP,x0,y0,x1,nstep = ode_init(stepper,nstep, version = 1, var = 0, time = time)
        x,y,it = fINT(fRHS,fORD,fBVP,x0,y0,x1,nstep,s=10,b=8/3,r=rparam)
        check(x,y,it, rparam,Name)
    elif ver == 3:
        if sound == 'record':
            mode = 2 #Record
            time = (nstep+1)/sample_rate #NEEDS TO BE DIVISABLE
            print(f'Message Time: {time}')

            #Initialization of sender generages chaotic signal
            fINT,fORD,fRHS,fBVP,x0,y0,x1,nstep = ode_init(stepper,nstep, version = 1, var = g, time = time)
            x,y,it = fINT(fRHS,fORD,fBVP,x0,y0,x1,nstep,s=10,b=8/3,r=rparam)

            #Obtains signal where s[0] = time and s[1] = your voice
            s = ad.Setter(Name, mode,T = time, rate = sample_rate)
       
            #signal = np.zeros(s[0].size)+1000
            signal = s[1]
            chaos = y[g]


            frac = 0.001
            #Normalization and reduction of signal to chaotic amplitude
            mag = frac*max(abs(y[g]))/(max(abs(s[1])))
            #mag = 1
            #Mask signal
            mask = signal*mag + chaos  
            nstep = len(s[1])-1
            
            #The Plots of signal vs. time, chaos vs. time, mask vs. time
            figure, axs = plt.subplots(3, 1, figsize=(30,4))
            plt.subplots_adjust(hspace=1)
            axs[0].plot(s[0], signal) 
            axs[0].set_xlim(s[0][0], s[0][-1])
            axs[0].set_xlabel('time (s)')
            axs[0].set_ylabel('amplitude')
            axs[0].set_title('Chaotic x')
            axs[1].plot(s[0], chaos) 
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

            #Plays Masked Signal
            print("masked signal")
            try:
                ad.play(mask*500, Name, 'masked', sample_rate, time)
            except:
                pass


            #Fourier Spectrum Plots
            figure, axes = plt.subplots(2, 1, figsize=(30,4))
            plt.subplots_adjust(hspace=1)

            frequencies, fourierTransform = fft(signal/max(signal), time, sample_rate) #Signal FFT
            #First subplot deals with graphing signal and chaotic fft
            axes[0].plot(frequencies, abs(fourierTransform), label = 'Signal')
            frequencies, fourierTransform = fft(chaos/max(chaos), time, sample_rate) #Chaotic FFT
            axes[0].plot(frequencies, abs(fourierTransform), label = 'Chaos')                       
            axes[0].set_xlabel('Frequency')
            axes[0].set_ylabel('Amplitude')
            axes[0].set_title('Signal')

            #2nd subplot graphs mask fft
            frequencies, fourierTransform = fft(mask, time, sample_rate) #Signal + Chaotic FFT
            axes[1].plot(frequencies, abs(fourierTransform))
            axes[1].set_xlabel('Frequency')
            axes[1].set_ylabel('Amplitude')
            axes[1].set_title('Sig+Chaos')
            plt.legend()
            plt.show()

            #This up the index of the drive signal that is required for the ode
            set_step()

            #ODE sync with the receiver system
            fINT,fORD,fRHS,fBVP,x0,y0,x1,nstep = ode_init(stepper,nstep, version = 2, var = g, time = time)
            x,yrec,it2 = fINT(fRHS,fORD,fBVP,x0,y0,x1,nstep,s=10,b=8/3,r=rparam, driving = mask)

        elif sound == 'binary':
            mode = 1 #Binary
            time = (nstep+1)/sample_rate
            cycpersec = 4 #Cycles/second


            #Sender system 
            fINT,fORD,fRHS,fBVP,x0,y0,x1,nstep = ode_init(stepper,nstep, version = 1, var = g, time = time)
            x,y,it = fINT(fRHS,fORD,fBVP,x0,y0,x1,nstep,s=10,b=8/3,r=rparam)

            #Generates binary signal according to frequency and sample rate
            s = ad.Setter(Name, mode,T = cycpersec,time = time, A = 30000,rate = sample_rate)

            frac = 0.1*(100000/nstep)

            signal = s[1]
            chaos = y[g]
            
            
            #Sets up the normalization of the signal respective to the chaos
            mag = frac*max(abs(y[g]))/(max(abs(s[1])))
            #Mask
            mask = signal*mag + chaos
            print(max(signal*mag) )
            nstep = len(s[1])-1

            #Plots Singla, chaos and mask
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
            plt.draw()
            #This up the index of the drive signal that is required for the ode
            set_step()

            #ODE sync with the receiver system
            fINT,fORD,fRHS,fBVP,x0,y0,x1,nstep = ode_init(stepper,nstep, version = 2, var = g, time = time)
            x,yrec,it2 = fINT(fRHS,fORD,fBVP,x0,y0,x1,nstep,s=10,b=8/3,r=rparam, driving = mask)


        
        #Obtains Recoverd and graphs
        recovered = (mask-yrec[g])/mag
        plt.figure(figsize=(30, 4))
        plt.plot(s[0], recovered)
        plt.title("Recovered Signal")
        plt.draw()
        
        #Plays Recovered
        if mode == 2:
            

            print("Recovered Signal")
            try:
                ad.play(recovered, Name, 'recovered', sample_rate, time )
            except:
                pass
        

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

    #Default Run : python Lorenz.py rk4 100000 3 Practice 30 record x
    '''
    Work: what values of r best synchronize the message
          Is there a single value we can use to characterize synchronizaiton?
          How does step size effect the convergence of the synchronization?
          Is there any parallel work or anything else we can use to speed up the integrator?
          Search with a cost function = integral of Lyaponov
    '''
    '''
    Time Dependance of Lyaponov Function and level of Synchronization
    Constant function failed syncrhonization why?
    Low levels marked similar to drive lyaponov difference but higher inten relationship lost 
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
        

        p2 = mp.Process(target=run, args=(stepper, nstep, ver, Name, rparam, sound, drive)) #Process for pertubration
        p1 = mp.Process(target=run, args=(stepper, nstep, 2, Name, rparam, 'NA', drive)) #Process for synchronization check
        p1.start()
        p2.start()
    else: 
        run(stepper, nstep, ver, Name, rparam, sound, drive )    


#=======================================
