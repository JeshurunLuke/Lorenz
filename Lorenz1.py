import time
import multiprocessing as mp
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

sigma = 10.0
beta = 8/3
tlen = 5
rho = 39
n = 1000
def lorenz(q, t, sigma, rho, beta):
    x, y, z = q
    return [sigma*(y - x), x*(rho - z) - y, x*y - beta*z]
def solve(ic):
    t = np.linspace(0, tlen, n)
    sigma = 10.0
    beta = 8/3

    sol = odeint(lorenz, ic, t, args=(sigma, rho, beta), rtol=1e-10, atol=1e-12)
    return sol

if __name__ == "__main__":
    Ntraj = 10
    np.random.seed(1)
    ics = -15 + 30 * np.random.random((Ntraj, 3))

    tstart = time.time()
    p = mp.Pool(5)
    solns = p.map(solve, ics)
    tend = time.time()
    print(f"multiprocessing: {tend-tstart} seconds")
  
    fig = plt.figure()
    ax = fig.add_axes([0, 0, 1, 1], projection='3d')
    colors = plt.cm.jet(np.linspace(0, 1, Ntraj))
    print(colors)