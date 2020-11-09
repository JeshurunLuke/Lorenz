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

prob = 2
folder = 'Prob1'
if prob == 4:
    rho = 28
else:
    rho = 39




sigma = 10.0
beta = 8/3
if prob == 1:
    tlen = 5
elif prob == 2:
    tlen = 30
elif prob == 3:
    tlen = 10
elif prob == 4:
    tlen = 30

n = 1600
def lorenz(q, t, sigma, rho, beta):
    x, y, z = q
    return [sigma*(y - x), x*(rho - z) - y, x*y - beta*z]
def solve(ic):
    t = np.linspace(0, tlen, n)
    sigma = 10.0
    beta = 8/3

    sol = odeint(lorenz, ic, t, args=(sigma, rho, beta), rtol=1e-10, atol=1e-12)
    return sol
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

    ax.view_init(30, 0.3 * i)
    fig.canvas.draw()
    return lines + pts

def update(num, data,line):
    line.set_data(data[:2, :num])
    line.set_3d_properties(data[2, :num])
if __name__ == "__main__":

    pather = os.getcwd() + '\\Data' + "\\" + folder
    try:
        os.mkdir(pather)
    except:
        pass
    if prob == 1: #This Problem graphs the cool 3d to show transient conditions don't matter
        Ntraj = 20
        np.random.seed(1)
        ics = -15 + 30 * np.random.random((Ntraj, 3))

        tstart = time.time()
        p = mp.Pool(5)
        solns = p.map(solve, ics)
        solns = np.array(solns)
        print(solns.shape)
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
        ax.set_xlim((-25, 25))
        ax.set_ylim((-35, 35))
        ax.set_zlim((5, 55))

        # set point-of-view: specified by (altitude degrees, azimuth degrees)
        ax.view_init(30, 0)
        anim = animation.FuncAnimation(fig, animate, init_func=init,fargs=[solns],
                                frames=500, interval=30, blit=True)

    # Save as mp4. This requires mplayer or ffmpeg to be installed
        plt.show()

        pat = pather + f"\\prob{int(prob)}.gif"
        anim.save(pat, writer='imagemagick', fps=30)

    elif prob == 2: #Illustrates Chaotic Nature from small change in conditions
        t = np.linspace(0, tlen, n)
        ics = np.array([[1.000,0,0],[1.0001,0,0]])
        p = mp.Pool(5)
        solns = p.map(solve, ics)
        plt.title('Chaos')
        plt.plot(t,solns[0][:,0])
        plt.plot(t,solns[1][:,0])
        plt.xlabel('time')
        plt.ylabel('x')
        plt.show()

    elif prob == 3: #Transients part 2
        t = np.linspace(0, tlen, n)
 
        ics = np.array([100,100,100])
        soln = solve(ics)
        plt.plot(t, soln[:,0])
        plt.title('Transients disappear')
        plt.xlabel('time')
        plt.ylabel('x')
        plt.show()
        
    elif prob == 4:
        Ntraj = 1
        t = np.linspace(0, tlen, n)
        ics = np.array([rho**(1/2),rho**(1/2),rho])
        soln = solve(ics)   
        fig = plt.figure()
        ax = p3.Axes3D(fig)

        col = ['black','green','cyan','blue','red','black','black','black','black']

        y = np.transpose(soln)
        N = y.shape[1]-1
        line, = ax.plot(y[0, :], y[1, :], y[2, :], '-', color=col[1])

        # Setting the axes properties
        ax.set_xlim3d([min(y[0])-1, max(y[0])+1])
        ax.set_zlim3d([min(y[2])-1, max(y[2])+1])
        ax.set_zlabel('Z')
        ax.legend()
        ani = animation.FuncAnimation(fig, update, N, fargs=(y, line), interval=1000/N, blit=False)
        plt.show()        
        pat = pather + f"\\prob{int(prob)}with{int(rho)}.gif"
        ani.save(pat, writer='imagemagick', fps=100)

