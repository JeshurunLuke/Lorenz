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
from scipy.integrate import odeint


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
    # var = int(kwargs['var'])

    fBVP = 0 # default is IVP, but see below.
    if (stepper == 'euler'):
        fORD = step.euler
    elif (stepper == 'rk2'):
        fORD = step.rk2
    elif (stepper == 'rk4'):
        fORD = step.rk4
    elif (stepper == 'rk45'):
        fORD = step.rk45
    # elif (stepper =='odeint'):
    #     fORD = odeint
    else:
        raise Exception('[ode_init]: invalid stepper value: %s' % (stepper))

    x1 = 63#time = 1

    x0 = 0
    fINT = ode_ivp
    if ver == 1:
        y0 = np.array([1e-16,1e-16,1e-16])
        fRHS = dydx_lorenz
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
    n = len(x)       
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    normt=plt.Normalize(min(x),max(x))(x)
    
    # invariant color plot
    # ax.plot(y[0, :], y[1, :], y[2, :],color='g',linewidth=1)

    # color varying with time plot
    for i in range(n-1):
        ax.plot(y[0,i:i+2], y[1,i:i+2], y[2,i:i+2], color=plt.cm.viridis(normt[i]))

    ax.set_xlim3d([min(y[0])-1, max(y[0])+1])
    ax.set_xlabel('X')

    ax.set_ylim3d([min(y[1])-1, max(y[1])+1])
    ax.set_ylabel('Y')

    ax.set_zlim3d([min(y[2])-1, max(y[2])+1])
    ax.set_zlabel('Z')
    ax.legend()
    
    # Lorenz Map
    plt.figure()
    Mn = []
    Mn1 = []
    i = 0
    j = -1
    while(i<n):
        while(y[2][i]<y[2][i+1]):
            i+= 1
            if i == n-2:
                break
        
        if i == n-2:
            break
    
        if  j == -1:
            Mn.append(y[2][i])
        
        else:
            Mn1.append(y[2][i])
        j *= -1
        while(y[2][i]>y[2][i+1]):
            i+=1
            if i==n-2:
                break
    x = np.linspace(0,300,100)        
    plt.plot(x, x, label='y = x',color='k')
    plt.scatter(Mn[0:len(Mn1)],Mn1,color='g',s =10)
    plt.xlim([0.5*max(Mn), max(Mn)+10])
    plt.ylim([0.7*max(Mn), max(Mn)+10])

    plt.xlabel('z_n')
    plt.ylabel('z_(n+1)')
    plt.legend()

# =============================================================================
#     # 2D presentation of the result
#     fig,axs = plt.subplots(3,1)
#     for i in range(n-1):
#         axs[0].plot(y[0,i:i+2], y[1,i:i+2], color=plt.cm.viridis(normt[i]))
#         axs[0].set_xlabel('x(t)')
#         axs[0].set_ylabel('y(t)')
#         
#         axs[1].plot(y[0,i:i+2],y[2,i:i+2],color = plt.cm.viridis(normt[i]))
#         axs[1].set_xlabel('x(t)')
#         axs[1].set_ylabel('z(t)')
#         
#         axs[2].plot(y[1,i:i+2],y[2,i:i+2],color = plt.cm.viridis(normt[i]))
#         axs[2].set_xlabel('y(t)')
#         axs[2].set_ylabel('z(t)')
# =============================================================================

def set_step():
    global steps
    steps = -1
    
def run(stepper, nstep, ver, Name, rparam, sound, drive):
    time = 100
    fINT,fORD,fRHS,fBVP,x0,y0,x1,nstep = ode_init(stepper,nstep, version = 1, var = 0, time = time)
    x,y,it = fINT(fRHS,fORD,fBVP,x0,y0,x1,nstep,s=10,b=8/3,r=rparam)
    
    if rparam == 14:
        y00=np.array([1e-16,-1e-16,1e-16])
        x1,y1,it = fINT(fRHS,fORD,fBVP,x0,y00,x1,nstep,s=10,b=8/3,r=rparam)
        
        fig = plt.figure()
        ax = fig.gca(projection='3d')
    
        ax.plot(y[0, :], y[1, :], y[2, :],color='g',linewidth=1)
        ax.plot(y1[0, :], y1[1, :], y1[2, :],color='r',linewidth=1)
    
        ax.set_xlim3d([min(y[0])-1, max(y[0])+1])
        ax.set_xlabel('X')
    
        ax.set_ylim3d([min(y[1])-1, max(y[1])+1])
        ax.set_ylabel('Y')
    
        ax.set_zlim3d([min(y[2])-1, max(y[2])+1])
        ax.set_zlabel('Z')
        plt.show()
    
    else:
        check(x, y, it, rparam, Name)
    
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
    run(stepper, nstep, ver, Name, rparam, sound, drive )    