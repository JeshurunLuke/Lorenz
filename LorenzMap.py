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
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.cm as cm


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
    else:
        raise Exception('[ode_init]: invalid stepper value: %s' % (stepper))
    x1 = 100#time = 1
    x0 = 0
    fINT = ode_ivp
    if ver == 1:
        y0 = np.array([10,10,10])
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
      
    # fig = plt.figure()
    # ax = p3.Axes3D(fig)
    n = len(x)       
    # points = np.array([y[0,:],y[1,:],y[2,:]]).T.reshape(-1,1,3)
    # segments = np.concatenate([points[:-1], points[1:]], axis=1)
    
    # cmap=plt.get_cmap('copper')
    # colors=[cmap(float(ii)/(n-1)) for ii in range(n-1)]
    
    #plot
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    
    # for ii in range(n-1):
    #     segii=segments[ii]
    #     lii,=ax.plot(segii[:,0],segii[:,1],segii[:,2],color=colors[ii],linewidth=2)
    #     #lii.set_dash_joinstyle('round')
    #     #lii.set_solid_joinstyle('round')
    #     lii.set_solid_capstyle('round')
    norm = plt.Normalize(min(y[0,:]),max(y[0,:]))(y[0,:])
    for i in range(n-1):
        ax.plot(y[0,i:i+2], y[1,i:i+2], y[2,i:i+2], color=plt.cm.viridis(norm[i]))

    ax.set_xlim3d([min(y[0])-1, max(y[0])+1])
    ax.set_xlabel('X')

    ax.set_ylim3d([min(y[1])-1, max(y[1])+1])
    ax.set_ylabel('Y')

    ax.set_zlim3d([min(y[2])-1, max(y[2])+1])
    ax.set_zlabel('Z')
    ax.legend()
    
    fig2, axs = plt.subplots(1,3)
    normx = plt.Normalize(min(y[0,:]),max(y[0,:]))(y[0,:n-2])
    normy = plt.Normalize(min(y[1,:]),max(y[1,:]))(y[0,:n-2])
    normz = plt.Normalize(min(y[2,:]),max(y[2,:]))(y[0,:n-2])
    
    axs[0].scatter(y[0,0:n-2],y[0,1:n-1],c = plt.cm.viridis(normx))
    axs[0].set_xlabel('x_n')
    axs[0].set_ylabel('x_(n+1)')
        
    axs[1].scatter(y[1,0:n-2],y[1,1:n-1],c = plt.cm.viridis(normy))
    axs[1].set_xlabel('y_n')
    axs[1].set_ylabel('y_(n+1)')
        
    axs[2].scatter(y[2,0:n-2],y[2,1:n-1],c = plt.cm.viridis(normz))
    axs[2].set_xlabel('z_n')
    axs[2].set_ylabel('z_(n+1)')

def set_step():
    global steps
    steps = -1
    
def run(stepper, nstep, ver, Name, rparam, sound, drive):
    time = 100
    fINT,fORD,fRHS,fBVP,x0,y0,x1,nstep = ode_init(stepper,nstep, version = 1, var = 0, time = time)
    x,y,it = fINT(fRHS,fORD,fBVP,x0,y0,x1,nstep,s=10,b=8/3,r=rparam)
    check(x,y,it, rparam,Name)
    plt.show()
    
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