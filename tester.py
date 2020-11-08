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
import Lorenz3 as lor
sigma = 10.0
beta = 8/3

def check(x,y):
      
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

    args   = parser.parse_args()
    stepper= args.stepper
    Name = args.namefolder
    global rho
    rho = args.r


    '''
    3 varables that need attention
    '''
    n = 10000 #Needs to be atleast above 100,000 for decent recovery
    tlen = 100 #Need atleast above 100 for decent masking and a certain ratio between n/tlen has to be maintained for convergence dt<0
    #But if you try tlen = 10 with binary it breaks horribly for r>24.8 why?

    if stepper == 'odeint':
        speed = True
    else:
        speed = False



    
    t = np.linspace(0, tlen, n)
    ics = np.array([50.0,50.0,50.0])


    tstart = time.time()
    if speed:
        sol1 = lor.solve(ics,rho,t)
        y = np.transpose(sol1)
    else:
        fINT,fORD,fRHS,fBVP = lor.ode_init(stepper,False)
        x,y,it = fINT(fRHS,fORD,fBVP,t[0],ics,t[t.size-1],n-1, False, rho = rho)  
    tend = time.time()

    tmp = tend - tstart
    print(f"Initialization took approximatly {tend-tstart} seconds")
    check(t,y)
    plt.show()