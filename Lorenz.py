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
from itertools import repeat
import os 
from argparse import RawTextHelpFormatter
import odestep as step
import Lorenz3 as lor

sigma = 10.0
beta = 8/3

def check(x,y,rho):

    n = len(x)       

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    norm = plt.Normalize(min(y[0,:]),max(y[0,:]))(y[0,:])
    print(f"It started: {rho}")
    for i in range(n-1):
        ax.plot(y[0,i:i+2], y[1,i:i+2], y[2,i:i+2], color=plt.cm.viridis(norm[i]))
    print(f"It ended: {rho}")
    ax.set_title(f'r = {rho}')
    ax.set_xlim3d([min(y[0])-1, max(y[0])+1])
    ax.set_xlabel('X')

    ax.set_ylim3d([min(y[1])-1, max(y[1])+1])
    ax.set_ylabel('Y')

    ax.set_zlim3d([min(y[2])-1, max(y[2])+1])
    ax.set_zlabel('Z')
    ax.legend()
    folder = 'Tests'
    
    pather = os.getcwd() + '\\Data' + "\\" + folder
    try:
        os.mkdir(pather)
    except:
        pass
    
    pat = pather + f"\\r{int(rho)}.png"
    print(pat)
    plt.savefig(pat)
    print('Yayay')
    plt.show()


def integ(rho, n,tlen):
    stepper = 'rk45'
    tlen = 100
    t = np.linspace(0, tlen, n)
    ics = np.array([50.0,50.0,50.0])
    fINT,fORD,fRHS,fBVP = lor.ode_init(stepper,False)
    x,y,it = fINT(fRHS,fORD,fBVP,t[0],ics,t[t.size-1],n-1, False, rho = rho) 
    check(x,y,rho)
def solve(ic,rho,t):
    sol = odeint(lor.lorenz, ic, t, args=(sigma, rho, beta), rtol=1e-6, atol=1e-6)
    return sol,rho
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument("stepper",type=str,default='euler',
                        help="stepping function:\n"
                             "   odeint: Production grade integrator\n"
                             "   rk4  : Runge-Kutta 2nd order\n"
                             "   rk45  : Runge-Kutta 4th order\n")

    parser.add_argument("namefolder",type=str,default=1,
                        help="paramaters tested")

    args   = parser.parse_args()
    stepper= args.stepper
    Name = args.namefolder



    '''
    3 varables that need attention
    '''
    n = 50000 #Needs to be atleast above 100,000 for decent recovery
    tlen = 100 #Need atleast above 100 for decent masking and a certain ratio between n/tlen has to be maintained for convergence dt<0
    #But if you try tlen = 10 with binary it breaks horribly for r>24.8 why?

    if stepper == 'odeint':
        speed = True
    else:
        speed = False





    t = np.linspace(0, tlen, n)
    ics = np.array([50.0,50.0,50.0])
    rho = [10, 22, 24.5, 100, 126.52, 150, 166.3, 212]

    



    if speed:
        print("Yo")
        num_processes = 7 
        p = mp.Pool(num_processes)
        mp_solutions = p.starmap(integ, zip(repeat(ics), rho, repeat(tlen)))
        for sol1 in mp_solutions:
            sol2,rho = sol1
            y = np.transpose(sol2)
            check(t,y,rho)
    else:
        num_processes = 7 
        p = mp.Pool(num_processes)
        mp_solutions = p.starmap(integ, zip(rho, repeat(n), repeat(tlen)))


    plt.show()