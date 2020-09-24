import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D

rho = 10

sigma = 10
beta = 8.0 / 3.0

def f(state, t):
    x, y, z = state  # Unpack the state vector
    return sigma * (y - x), x * (rho - z) - y, x * y - beta * z  # Derivatives

state0 = [0, 0 , 0]
t = np.arange(0.0, 100.0, 0.01)

states = odeint(f, state0, t)

fig = plt.figure()
ax = fig.gca(projection="3d")
ax.plot(states[:, 0], states[:, 1], states[:, 2])

state0 = np.array([np.sqrt(beta*(rho-1)), np.sqrt(beta*(rho-1)), rho-1]) #+ np.array([0.01,0.01,0.01])
t = np.arange(0.0, 40.0, 0.01)
states = odeint(f, state0, t)

ax = fig.gca(projection="3d")
ax.plot(states[:, 0], states[:, 1], states[:, 2], color = 'green')
plt.draw()
plt.show()