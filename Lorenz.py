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
import numpy as np
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import animation
import os



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
def dydx_lorenz(x,y,z,**kwargs):
    dydx    = np.zeros(3)
    sigma, b, r = param(x,**kwargs)
    # ????? from here
    dydx[0] = sigma*(y[1]-y[0])
    dydx[1] = r*y[0]-y[1]-y[0]*y[2]
    dydx[2] = y[0]*y[1]-b*y[2]
    # ????? to here
    return dydx
def dydx_lorenz6(x,y,z,**kwargs):
    dydx    = np.zeros(6)
    sigma, b, r = param(x,**kwargs)
    # ????? from here
    dydx[0] = sigma*(y[1]-y[0])
    dydx[1] = r*y[0]-y[1]-y[0]*y[2]
    dydx[2] = y[0]*y[1]-b*y[2]
    return dydx
    # ????? to her

#========================================
# fLOA for the string problem
# This is the (boundary) loading function.
# Should return the integration boundaries x0,x1, and
# the initial values y(x=0). Therefore, y must be
# an array with three elements. 
# The function takes an argument v, which in our case
# is just the eigenvalue lambda.
# Note that some of the initial values in y can
# (and should!) depend on lambda.
"""
def load_IC(): #Invalid Function atm
    # ????? from here
    y[0] = 10
    y[1] = 10
    y[2] = 10
    # ????? to here
    return y
"""
#========================================
# fSCO for the string problem
# This is the scoring function. Should return
# the function that needs to be zeroed by the
# root finder. In our case, that's just y(1). 
# The function takes the arguments x and y, where
# y is an array with three elements.
"""
def score_IC(x,y): #Invalid Function atm
    # ????? from here
    score = y[0,y.shape[1]-1] - 0
    # ????? to here
    return score # displacement should be zero.
"""


def ode_init(stepper,nstep):

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
    y0 = np.array([10,10,10])
    x1 = 100
    x0 = 0
    fINT = ode_ivp
    fRHS = dydx_lorenz
    fBVP = 0
    return fINT,fORD,fRHS,fBVP,x0,y0,x1,nstep

#========================================
# Single rk4 step.
# Already provided. Good to go; 2.5 free points :-)

#========================================
# ODE IVP driver.
# Already provided. Good to go; 2.5 free points :-)
def ode_ivp(fRHS,fORD,fBVP,x0,y0,x1,nstep,**kwargs):

    nvar    = y0.size                      # number of ODEs
    x       = np.linspace(x0,x1,nstep+1)   # generates equal-distant support points
    y       = np.zeros((nvar,nstep+1))     # result array 
    y[:,0]  = y0                           # set initial condition
    dx      = x[1]-x[0]                    # step size
    it      = np.zeros(nstep+1)
    for k in range(1,nstep+1):
        y[:,k],it[k] = fORD(fRHS,x[k-1],y[:,k-1],dx,**kwargs)
    return x,y,it


#=======================================
# A single trial shot.
# Sets the initial values (guesses) via fLOA, calculates 
# the corresponding solution via ode_ivp, and returns 
# a "score" via fSCO, i.e. a value for the rootfinder to zero out.
"""
def bvp_shoot(fRHS,fORD,fLOA,fSCO,v,nstep,**kwargs): #Assistor to root finder atm invalid root finder
    # ????? from here
    x0,x1,y = fLOA(v)
    x,y = ode_ivp(fRHS,fORD,x0,x1,y,nstep,**kwargs)

    score = fSCO(x,y)
  
    # ????? to here
    return score # this should be zero, and thus can be directly used.

#=======================================
# The rootfinder.
# The function pointers are problem-specific (see main()). 
# v0 is the initial guess for the eigenvalue (in our case).
# Should return x,y, so that the solution can be plotted.
def bvp_root(fRHS,fORD,fLOA,fSCO,v0,nstep,**kwargs): # Possible Root Finder atm invalid method

    eig = [v0]
    n = np.arange(1,100000,1) #Sets maximum iteration amount
    for elem in n: #Iterates through v1 by multiplying v0*1.05 till there is a sign change
        v1 = v0*(1.05)**elem
        score2 = bvp_shoot(fRHS,fORD,fLOA,fSCO,v0,nstep,**kwargs)
        score3 = bvp_shoot(fRHS,fORD,fLOA,fSCO,v1,nstep,**kwargs)
        if score2*score3<0:
            eig.append(v1)
            break

    cond = True
    v1 = eig[1]

    iter = 1
    errcon  = 1.89e-4 #Desirable Tolarance
    maxit   = 100000 #Max iteration amount (stops infinite loops)
    while cond and (iter<maxit): #Normal Bisection routine -> finds eigenvalue
        v = (v0+v1)/2
        score1 = bvp_shoot(fRHS,fORD,fLOA,fSCO,v,nstep,**kwargs)
        score2 = bvp_shoot(fRHS,fORD,fLOA,fSCO,v0,nstep,**kwargs)
        score3 = bvp_shoot(fRHS,fORD,fLOA,fSCO,v1,nstep,**kwargs)
        if score1*score2 > 0:
            v0 = v
        else:
            v1 = v
        if abs(score1) <= errcon:
            cond = False
        iter +=1
    x0,x1,y = fLOA(v)
    x,y = ode_ivp(fRHS,fORD,x0,x1,y,nstep,**kwargs)

    # ????? to here
    return x,y
"""
#=======================================
def check(x,y,it,r,iter,Name):
    col = ['black','green','cyan','blue','red','black','black','black','black']
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    
    data = y
    N = data.shape[1]-1
    line, = ax.plot(data[0, 0:1], data[1, 0:1], data[2, 0:1], '-o', color=col[iter-1])

    # Setting the axes properties
    ax.set_xlim3d([min(data[0])-1, max(data[0])+1])
    ax.set_xlabel('X')

    ax.set_ylim3d([min(data[1])-1, max(data[1])+1])
    ax.set_ylabel('Y')

    ax.set_zlim3d([min(data[2])-1, max(data[2])+1])
    ax.set_zlabel('Z')
    ax.legend()
    ani = animation.FuncAnimation(fig, update, N, fargs=(data, line), interval=10000/N, blit=False)
    pl = os.getcwd() + "\\Data" + "\\" + str(Name) + "\\r" + str(r)+'.gif'
    ani.save(pl, writer='imagemagick')


def update(num, data,line):
    line.set_data(data[:2, :num])
    line.set_3d_properties(data[2, :num])
    
    
def main():


    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument("stepper",type=str,default='euler',
                        help="stepping function:\n"
                             "   euler: Euler step\n"
                             "   rk2  : Runge-Kutta 2nd order\n"
                             "   rk4  : Runge-Kutta 4th order\n")
    parser.add_argument("steps",type=int,default=100,
                        help="number of steps")
    parser.add_argument("limit",type=int,default=1,
                        help="paramaters tested")
    parser.add_argument("namefolder",type=str,default=1,
                        help="paramaters tested")


    args   = parser.parse_args()

    stepper= args.stepper
    nstep = args.steps
    limit = args.limit
    Name = args.namefolder
    pather = os.getcwd() + '\\Data' + "\\" + Name
    try:
        os.mkdir(pather)
    except:
        pass
    fINT,fORD,fRHS,fBVP,x0,y0,x1,nstep = ode_init(stepper,nstep)
    rpar = np.array([10,22,24.5,100,126.52,150,166.3,212])

    tests = 1
    for rparam in rpar:
        print('Testing paramater: ' + str(rparam))
        print(str((tests)/limit)+"% Done ")
        x,y,it = fINT(fRHS,fORD,fBVP,x0,y0,x1,nstep,s=10,b=8/3,r=rparam)
        check(x,y,it,rparam,tests, Name)
        if tests == limit:
            break
        tests +=1
    plt.show()
    




#=======================================
main()
