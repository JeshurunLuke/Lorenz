#Setup: conda install pyaudio

'''
What does this Function do?: The Lorenz2.py was used to generate all the plots for 
the synchronization portion of the report.
Plots include:
    Optimization for n: steps required for 10^-6 (EXTENT OF CONVERGENCE)
    Optimizaiton for r and ic: (OVERALL PERFORMACE)
    Optimization for r: (RATE OF CONVERGENCE)

Arguments:
No command line arguments
Look at ___name___ = ____main____ for arguments



Output:
Problem Dependant
'''



import numpy as np
import time as time
import matplotlib.pyplot as plt

import multiprocessing as mp
from itertools import repeat
import Lorenz3 as lor
# Constants
sigma = 10.0
beta = 8/3
num_processes = 7  # Number of processors you want to use to generate data: has to be less than your # of logical processors





###
# Unperturbed function
# Performance syncriozation with signal = 0 to check for convergence
def Unperturbed(sol1,rho, ic,t, stepper,speed):
    signal = np.zeros(sol1[:,0].size)
    driving = signal + sol1[:,0]

    if speed == False: #Receiver differential equations with our integrator
        tstart = time.time()

        fINT,fORD,fRHS,fBVP = lor.ode_init(stepper, True)
        x,yrec,it = fINT(fRHS,fORD,fBVP,t[0],ic,t[t.size-1],signal.size-1,True, driving = driving, rho = rho)  

        tend = time.time()
        print(f"{rho}: It took approximatly {tend-tstart} seconds")
    else: #Receiver differential equation with scipy integrator
        tstart = time.time()
        sol2 = lor.solvrec(ic,rho, driving,t)
        yrec = np.transpose(sol2)
        tend = time.time()
        print(f"{rho}: It took approximatly {tend-tstart} seconds")

    y = np.transpose(sol1) 

    return y,yrec


##Looks at extent of integration used to find level of synchronization for different n
def Problem2(n,rho,speed,tlen):
    stepper = 'rk45'

    t = np.linspace(0, tlen, n) #t space for Chaotic system normally array from 0 to 100
    ics = np.ones(3)*50.0 #Initial Conditions

    if speed == True: #Uses scipy integrator to get trajectory: Driving
        ics = list(ics)
        sol1 = lor.solve(ics,rho,t)

    else: #Uses our integrators to solve for trajectories: Driving
        ics = np.array(ics)
        fINT,fORD,fRHS,fBVP = lor.ode_init(stepper,False)
        x,y,it = fINT(fRHS,fORD,fBVP,t[0],ics,t[t.size-1],n-1, False, rho = rho)  
        sol1 = np.transpose(y)
    #Preturbs initial conditions by a vector [10,10,10]
    ics += np.ones(3)*10
    y, yrec = Unperturbed(sol1,rho, ics,t, stepper,speed) #Solves receiver with a driving = 0

    #Gets the error
    error = np.abs(y[0]-yrec[0])
    
    #Gets the halfway point length 
    ss = y[0][n//2:-1].size
    #finds the minimium value as an approximation to syncrhoization time after the half way point
    m = np.min(np.abs(y[0][n//2:-1]-yrec[0][n//2:-1]))
    
    location = list(np.abs(y[0][n//2:-1]-yrec[0][n//2:-1])).index(m)
    location += ss
    
    #Averages from minimimum to end of integration time
    m = np.average(error[location:-1])

    return n,m, location, error

##Looks at Extent and Rate of integration with an average used to find level of synchronization for different r and initial conditions
def Problem3(rho,n,ic,speed,tlen):
    print(f'{rho} with {ic} intialization')
    stepper = 'rk45'

    t = np.linspace(0, tlen, n)
    if speed: #Uses scipy integrator to get trajectory: Driving
        ics = list(ic)
        sol1 = lor.solve(ics,rho,t)

    else: #Uses our integrators to solve for trajectories: Driving
        ics = np.array(ic)
        fINT,fORD,fRHS,fBVP = lor.ode_init(stepper,False)
        x,y,it = fINT(fRHS,fORD,fBVP,t[0],ics,t[-1],n-1, False, rho = rho)  
        sol1 = np.transpose(y)
    print(f'{rho} with {ic} Synchronizaiton')

    #Preturbs initial conditions by a vector [10,10,10]
    ics += np.ones(3)*10
    y, yrec = Unperturbed(sol1,rho, ics,t, stepper,speed)  #Solves receiver with a driving = 0
    
    #Gets the error and averages   
    error = np.abs(y[0]-yrec[0])
    m = np.average(error)

    return rho,m, error

#Gets lyaponov and lyaponov derivate evaluation
def lyap(x, y, yrec):
    e = np.zeros((3,y[0].size))
    e[0] = y[0]-yrec[0]
    e[1] = y[1]-yrec[1]
    e[2] = y[2]-yrec[2]
    

    Ldot = -(e[0]-0.5*e[1])**2-0.75*(e[1]**2)-4*beta*(e[2]**2) #Lyaponev Derivative

    L = 0.5*(1/sigma)*(e[0]**2) + e[1]**2 + e[2]**2 #Lyaponev Function

    return L,Ldot


##Looks at Rate of integration with an average used to find level of synchronization for different r 
def Problem4(rho,n,speed,tlen):
    stepper = 'rk45'

    t = np.linspace(0, tlen, n)
    ics = np.ones(3)*50.0
    if speed:#Uses scipy integrator to get trajectory: Driving
        ics = list(ics)

        sol1 = lor.solve(ics,rho,t)

    else: #Uses our integrators to solve for trajectories: Driving
        ics = np.array(ics)
        fINT,fORD,fRHS,fBVP = lor.ode_init(stepper,False)
        x,y,it = fINT(fRHS,fORD,fBVP,t[0],ics,t[t.size-1],n-1, False, rho = rho)  
        sol1 = np.transpose(y)
    ics += np.ones(3)*1000
    y, yrec = Unperturbed(sol1,rho, ics,t, stepper,speed) #Solves receiver with a driving = 0
    
    #Gets the Lyaponov and checks for when its below 0.5 as a way to quanity rate of rxn
    L,Ldot = lyap(t, y, yrec)
    print(f"Initial Lyaponov value: {L[0]}")
    for i in range(0,L.size):
        if L[i] < 0.5:
            rate = np.average(Ldot[0:i])
            time = t[i]
            break
    return rho,L,rate, time


'''
Default:
    ndependance = True
    icrdependance = False
    rdependance = False
'''
if __name__ == "__main__":
    ndependance = False ##Looks at extent of integration used to find level of synchronization for different n
    icrdependance = True ##Looks at Extent and Rate of integration with an average used to find level of synchronization for different r and initial conditions
    rdependance = False ##Looks at Rate of integration with an average used to find level of synchronization for different r and initial conditions

    tlen = 100
    #If speed = False --> rk45 integrator else odeint integrator
    speed = True


    #If speed = False --> rk45 integrator else odeint integrator

    if ndependance:

        rho = 60 #Sets r
        n = [100,1000,5000,10000,50000,100000,1000000,2500000,5000000]

        #Uses parallel run to run all the different n systems
        print("multiprocessing:", end='')
        tstart = time.time()
        p = mp.Pool(num_processes)
        mp_solutions = p.starmap(Problem2, zip(n, repeat(rho),repeat(speed),repeat(tlen))) 
        #serial_solutions = [Problem2(ic) for ic in n]
        tend = time.time()
        tmp = tend - tstart
        print(" %8.3f seconds" % tmp)


        narray = []
        minarray = []
        for sol1 in mp_solutions:
            n, minim, location, error = sol1
            narray.append(n)
            minarray.append(minim)
            print(n,location/n*100, minim)
            #Plots end of error if the requirement is met
            if minim < 5*10**(-6):
                t = np.linspace(0, tlen, n)

                plt.plot(t[location:-1],error[location:-1])
                plt.title(f'{n} error')
                plt.show()
            
        plt.plot(np.log(narray),minarray,color='r')
        stri = f'Parameters:\n Sigma = {sigma}\n beta = {round(beta,3)}\n r = {rho}'
        plt.text(13, 13,stri, style='italic', bbox = {'facecolor': 'white'})
        plt.title('Performance of system for various Steps')
        plt.xlabel('log(Step)')
        plt.ylabel('Minimal Error')
        plt.show()
    if icrdependance:

        ics = [[0.75,0.75,0.75],[5.0,5.0,5.0],[10.0,10.0,10.0],[100.0,100.0,100.0],[50.0,50.0,50.0]]
        n = 100000 #Sets number of steps
        #Iterates through each IC
        for ic in ics:
            print(f"IC in progress:{ic}")
            rho = np.linspace(25,70,15)

            ic2 = tuple(ic)
            
            #Uses parallel run to run all the different r systems
            print("multiprocessing:", end='')
            tstart = time.time()
            p = mp.Pool(num_processes)
            mp_solutions = p.starmap(Problem3, zip(rho, repeat(n), repeat(ic2),repeat(speed),repeat(tlen)))
            #serial_solutions = [Problem2(ic) for ic in n]
            tend = time.time()
            tmp = tend - tstart
            print(" %8.3f seconds" % tmp)



            rhoarray = []
            minarray = []
            for sol1 in mp_solutions:
                rho, minim, error = sol1
                rhoarray.append(rho)
                minarray.append(minim)

            plt.plot(rhoarray,minarray, label=f'IC = {ic}')
        plt.legend()
        stri = f'Parameters:\n Sigma = {sigma}\n beta = {round(beta,3)}'
        plt.text(25, 0.028,stri, style='italic', bbox = {'facecolor': 'white'})
        plt.title('Performance of system for various R')
        plt.xlabel('rho')
        plt.ylabel('Drive Function Average Error')
        plt.show()
    print("To")
    if rdependance: #Checks rate for different r 
        rho = np.arange(25,70,5) 
        n = 1000000
        
        #Uses parallel run to run all the different r systems
        print("multiprocessing:", end='')
        tstart = time.time()
        p = mp.Pool(num_processes)
        mp_solutions = p.starmap(Problem4, zip(rho, repeat(n),repeat(speed),repeat(tlen)))
        tend = time.time()
        tmp = tend - tstart
        print(" %8.3f seconds" % tmp)

        rate = []
        rhoarray = []
        times = []
        for sol1 in mp_solutions:
            rho1,L,percent,tc = sol1
            rate.append(percent)
            rhoarray.append(rho1)
            times.append(tc)

        
        plt.plot(rhoarray,np.abs(rate))
        stri = f'Parameters:\n Sigma = {sigma}\n beta = {round(beta,3)}'
        plt.text(190, 3,stri, style='italic', bbox = {'facecolor': 'white'})
        plt.title('Rate of Convergence of system for various R')
        plt.xlabel('rho')
        plt.ylabel('Average Lyaponov Derivative')
        plt.show()
        plt.plot(rhoarray,times)
        stri = f'Parameters:\n Sigma = {sigma}\n beta = {round(beta,3)}'
        plt.text(190, 3,stri, style='italic', bbox = {'facecolor': 'white'})
        plt.title('Rate of Convergence of system for various R')
        plt.xlabel('rho')
        plt.ylabel('Time of Convergence')
        plt.show()