#Setup: conda install pyaudio

'''
What does this Function do?: The Lorenz2.py was used to generate all the plots for 
the presentations
Plots include:
    Animation of lorenz system
    Chaotic Nature
    Transient Conditions
Arguments:
No command line arguments
Look at ___name___ = ____main____ for arguments



Output:
Problem Dependant
'''


import time
import multiprocessing as mp
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import mpl_toolkits.mplot3d.axes3d as p3
import os as os
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import cnames
from matplotlib import animation
import Lorenz3 as lor
from itertools import repeat
num_processes = 7  # Number of processors you want to use to generate data: has to be less than your # of logical processors

'''
IMPORTANT IF YOU HAVE A MAC OR UNIX CHANGE ALL THE // to \ to properly save the file
'''
def solverk45(ics, rho,t):
    n = t.size
    stepper = 'rk45'
    fINT,fORD,fRHS,fBVP = lor.ode_init(stepper,False)
    x,y,it = fINT(fRHS,fORD,fBVP,t[0],ics,t[t.size-1],n-1, False, rho = rho)  
    sol1 = np.transpose(y)
    return sol1


def init():
    for line, pt in zip(lines, pts):
        line.set_data([], [])
        line.set_3d_properties([])

        pt.set_data([], [])
        pt.set_3d_properties([])
    return lines + pts

def animate(i,x_t):
    # we'll step two time-steps per frame.  This leads to nice results.
    i = (2 * i) % x_t.shape[1]
    
    for line, pt, xi in zip(lines, pts, x_t):
        x, y, z = xi[:i].T
        line.set_data(x, y)
        line.set_3d_properties(z)

        pt.set_data(x[-1:], y[-1:])
        pt.set_3d_properties(z[-1:])

    #ax.view_init(30, 0.3 * i)
    fig.canvas.draw()
    return lines + pts

def update(num, data,line):
    line.set_data(data[:2, :num])
    line.set_3d_properties(data[2, :num])
if __name__ == "__main__":
    prob = 1    #Can be one of 3

    folder = 'Prob1' #Folder where your files are stored
    speed = False
    n = 100000


    rho = 100
    sigma = 10.0
    beta = 8/3


    if prob == 1:
        tlen = 60
    elif prob == 2:
        tlen = 30
    elif prob == 3:
        tlen = 10
    
    pather = os.getcwd() + '\\Data'
    try:
        os.mkdir(pather)
    except:
        pass
 

    pather = os.getcwd() + '\\Data' + "\\" + folder
    try:
        os.mkdir(pather)
    except:
        pass

    t = np.linspace(0, tlen, n)

    if prob == 1: #This Problem graphs the cool 3d to show transient conditions don't matter
        Ntraj = 2
        ics = -15 + 30 * (2*np.random.random((Ntraj, 3))-1)

        tstart = time.time()
        p = mp.Pool(num_processes)
        if speed:
            solns = p.starmap(lor.solve, zip(ics, repeat(rho),repeat(t)))
        else:
            solns = p.starmap(solverk45, zip(ics, repeat(rho),repeat(t)))
        solns = np.array(solns)
        tend = time.time()
        print(f"multiprocessing: {tend-tstart} seconds")
    
        # Set up figure & 3D axis for animation
        fig = plt.figure()
        ax = fig.add_axes([0, 0, 1, 1], projection='3d')
        ax.axis('off')
        colors = plt.cm.jet(np.linspace(0, 1, Ntraj))

        # set up lines and points
        lines = sum([ax.plot([], [], [], '-', c=c)
                    for c in colors], [])
        pts = sum([ax.plot([], [], [], 'o', c=c)
                for c in colors], [])

        # prepare the axes limits
        print(solns.shape)
        print(min(solns[0,0,:]))
        ax.set_xlim((min(solns[0,:,0]), max(solns[0,:,0])))
        ax.set_ylim((min(solns[0,:,1]), max(solns[0,:,1])))
        ax.set_zlim((min(solns[0,:,2]), max(solns[0,:,2])))

        # set point-of-view: specified by (altitude degrees, azimuth degrees)
        ax.view_init(0, 0)
        anim = animation.FuncAnimation(fig, animate, init_func=init,fargs=[solns],
                                frames=500, interval=30, blit=True)

    # Save as mp4. This requires mplayer or ffmpeg to be installed
        plt.show()

        pat = pather + f"\\prob{int(prob)}.gif"
        anim.save(pat, writer='imagemagick', fps=30)

    elif prob == 2: #Illustrates Chaotic Nature from small change in conditions
        t = np.linspace(0, tlen, n) #Time for lorenz system
        ics = np.array([[1.000,0.001,0.1],[1.000,0.002,0.1]]) #Initial condition small differnece
        p = mp.Pool(num_processes)
        if speed: #scipy integrator
            solns = p.starmap(lor.solve, zip(ics, repeat(rho),repeat(t)))
        else: #Our integrator
            solns = p.starmap(solverk45, zip(ics, repeat(rho),repeat(t)))
        plt.title('Chaos')
        plt.plot(t,solns[0][:,0])
        plt.plot(t,solns[1][:,0])
        plt.xlabel('time')
        plt.ylabel('x')
        plt.show()

    elif prob == 3: #Transients part 2
        t = np.linspace(0, tlen, n)
 
        ics = np.array([100,100,100])
        if speed: #Uses scipy integrator
            soln = lor.solve(ics,rho,t)
        else: #Uses our integrator 
            solverk45(ics, rho,t)

        plt.plot(t, soln[:,0])
        plt.title('Transients disappear')
        plt.xlabel('time')
        plt.ylabel('x')
        plt.show()
        
