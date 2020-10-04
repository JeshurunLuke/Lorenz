import numpy as np
import odestep as step
import utilities as util
import argparse	                 # allows us to deal with arguments to main()
from argparse import RawTextHelpFormatter
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
from matplotlib import pyplot as plt
from scipy.fftpack import fft
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import animation
import os
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
    # ????? to her
def get_step():
    global steps
    steps +=1
    return steps



def ode_init(stepper,nstep, **kwargs):
    ver = int(kwargs['version'])

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
    x1 = 100
    x0 = 0
    fINT = ode_ivp
    if ver == 1:
        y0 = np.array([10,10,10])
        fRHS = dydx_lorenz
    elif ver ==2:
        y0 = np.array([10,10,10]) + np.array([0.1,0.1,0.1])
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

    for key in kwargs:
        if (key=='receiving'):
            version = True
            yrec = kwargs[key]
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    

    ax.plot(y[0, :], y[1, :], y[2, :], '-o', color=col[1])
    ax.plot(yrec[0, :], yrec[1, :], yrec[2, :], '-o', color=col[4])
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
    L = -(e[0]-0.5*e[1])**2-0.75*e[1]**2-4*b*e[2]**2
    plt.plot(x, L)
    return L
def main():


    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument("stepper",type=str,default='euler',
                        help="stepping function:\n"
                             "   euler: Euler step\n"
                             "   rk2  : Runge-Kutta 2nd order\n"
                             "   rk4  : Runge-Kutta 4th order\n")
    parser.add_argument("steps",type=int,default=100,
                        help="number of steps")
    parser.add_argument("ver",type=int,default=1,
                        help="1:Lorenz\n" 
                        "   2:Unperturbed")
    parser.add_argument("namefolder",type=str,default=1,
                        help="paramaters tested")
    parser.add_argument("r",type=str,default=1,
                        help="r tested")
    parser.add_argument("sound",type=str,default=1,
                        help="record: recording\n"
                          "     binary: binary\n" 
                          "     NA: for version 1 and 2" )

    args   = parser.parse_args()

    stepper= args.stepper
    nstep = args.steps-1
    ver = args.ver
    Name = args.namefolder
    rparam = args.r
    sound = args.sound

    pather = os.getcwd() + '\\Data' + "\\" + Name
    '''
    try:
        os.mkdir(pather)
    except:
        pass
    '''




    if ver ==2:
        fINT,fORD,fRHS,fBVP,x0,y0,x1,nstep = ode_init(stepper,nstep, version = 1)
        x,y,it = fINT(fRHS,fORD,fBVP,x0,y0,x1,nstep,s=10,b=8/3,r=rparam)

        set_step()

        fINT,fORD,fRHS,fBVP,x0,y0,x1,nstep = ode_init(stepper,nstep, version = ver)
        x,yrec,it2 = fINT(fRHS,fORD,fBVP,x0,y0,x1,nstep,s=10,b=8/3,r=rparam, driving = y[0])
        L = error(x,y,yrec,10,8/3)
        check(x,y,it, rparam,Name,  receiving = yrec)
    elif ver == 1:
        fINT,fORD,fRHS,fBVP,x0,y0,x1,nstep = ode_init(stepper,nstep, version = 1)
        x,y,it = fINT(fRHS,fORD,fBVP,x0,y0,x1,nstep,s=10,b=8/3,r=rparam)
        check(x,y,it, rparam,Name,  receiving = y)
    elif ver == 3:

        if sound == 'record':
            fINT,fORD,fRHS,fBVP,x0,y0,x1,nstep = ode_init(stepper,nstep, version = 1)
            x,y,it = fINT(fRHS,fORD,fBVP,x0,y0,x1,nstep,s=10,b=8/3,r=rparam)
            mode = 2
            sample_rate = 50000
            time = (nstep+1)/sample_rate

            s = ad.Setter(Name, mode,T = time, rate = sample_rate)
            mask = s[1]/(max(s[1])) + y[0]
            nstep = len(s[1])-1
            plt.plot(s[0], mask)
            plt.show()
            ad.play(mask*1000, Name, 'masked' )
            fft_output = fft(mask)
            plt.plot(mask, np.abs(fft_output))
            plt.show()

            set_step()

            fINT,fORD,fRHS,fBVP,x0,y0,x1,nstep = ode_init(stepper,nstep, version = 2)
            x,yrec,it2 = fINT(fRHS,fORD,fBVP,x0,y0,x1,nstep,s=10,b=8/3,r=rparam, driving = mask)
        elif sound == 'binary':
            mode = 1
            sample_rate = 50000
            time = (nstep+1)/sample_rate
            cycpersec = 2 #Cycles/second
            s = ad.Setter(Name, mode,T = cycpersec,time = time, A = 30000,rate = sample_rate)

            nstep = len(s[1])-1

            fINT,fORD,fRHS,fBVP,x0,y0,x1,nstep = ode_init(stepper,nstep, version = 1)
            x,y,it = fINT(fRHS,fORD,fBVP,x0,y0,x1,nstep,s=10,b=8/3,r=rparam)


            mask = s[1]/(max(s[1])) + y[0]
            plt.plot(s[0], mask)
            plt.show()
            fINT,fORD,fRHS,fBVP,x0,y0,x1,nstep = ode_init(stepper,nstep, version = 2)
            set_step()

            x,yrec,it2 = fINT(fRHS,fORD,fBVP,x0,y0,x1,nstep,s=10,b=8/3,r=rparam, driving = mask)
        
        recovered = (mask-yrec[0])*(max(s[1]))
        if mode == 2:
            ad.play(recovered, Name, 'recovered' )
        plt.plot(s[0], recovered)
        plt.show()

        

    plt.show()
    




#=======================================
main()
