#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Of bowls of soup and girlfriends. Works with python 2 and 3.
Dependencies: numpy, matplotlib
"""

import numpy as np
import pylab as pl

def SA(f, theta_0 = 0., a_0=10, n_steps=1000):
    """
    Stochastic approximation using the Robbins-Monro (1951) method
    """
    theta = theta_0
    a = float(a_0)
    theta_history = [theta_0]
    step_decrease_rate = 0.5

    n = 1
    i_steps = 1
    while i_steps < n_steps+1:
        theta = theta - (a/n**step_decrease_rate)*f(theta)

        theta_history.append(theta)
        i_steps += 1
        n += 1

    return theta_history

def SA2(f, theta_0 = 0., a_0=10, n_steps=1000):
    """
    Stochastic approximation using the Robbins-Monro (1951) method,
    Kesten-Delyon acceleration and Polyak-Ruppert averaging
    """
    theta = theta_0
    av_theta = theta_0
    y_prev = 0
    a = float(a_0)
    theta_history = [theta_0]
    av_theta_history = [theta_0]
    step_decrease_rate = 0.5

    n = 1
    i_steps = 1
    averaged_steps = 1
    while i_steps < n_steps+1:
        y = f(theta)
        theta = theta - (a/n**step_decrease_rate)*y
        theta_history.append(theta)

        # Polyak-Ruppert averaging
        # Start averaging after crossing the zero of the functions an arbitrary
        # number of times (20, in this case)
        if n == 20:
            av_theta = theta
        elif n > 20:
            av_theta = (float(averaged_steps)/(averaged_steps+1)*av_theta + 
                        1./(averaged_steps+1)*theta)
            averaged_steps += 1

        av_theta_history.append(av_theta)

        # Kesten-Delyon accel.; change of sign suggests step size is too large
        n += y*y_prev < 0
        y_prev = y

        i_steps += 1

    return theta_history, av_theta_history

def noiseless_function(x):
    """
    Target function. Root at 5
    """
    return 1/(1+np.exp(-x+5))-0.5

def perturb(x):
    """
    Noise in gf's decision. If you want to cap the salt so that it doesn't go
    negative, you should do it here.
    """
    gf_crazyness = 0.2
    try:
        # x is a vector
        return x+np.random.normal(scale=gf_crazyness, size=len(x))
    except:
        # x is a scalar
        return x+np.random.normal(scale=gf_crazyness) 

def noisy_function(x):
    return perturb(noiseless_function(x))

theta_0 = 100       # Initial estimation
n_steps = int(1e2)  # Realistic number of bowls of soup before she leaves me
a_0 = 100           # Initial step size

res = SA(noisy_function, theta_0 = theta_0, a_0 = a_0, n_steps = n_steps)
res2, av_res2 = SA2(noisy_function, theta_0 = theta_0, a_0 = a_0, 
                    n_steps = n_steps)

pl.subplot(1, 2, 1)
xs = np.linspace(0, 10, 1000)
pl.plot(xs, noiseless_function(xs))
pl.plot(xs, noisy_function(xs), '+')
pl.xlabel('Amount of salt (arbitrary units)')
pl.ylabel('(Abs): verbal level of dissatisfaction (a.u.)\n(Sign): Reason')
pl.legend(['Underlying noiseless function', 'Noisy samples'], loc='best')

print res[-1], res2[-1], av_res2[-1]

pl.subplot(1, 2, 2)
pl.plot(res, alpha=0.4)
pl.plot(res2, alpha=0.4)
pl.plot(av_res2, alpha=0.4)
pl.xlabel('Bowls of soup')
pl.ylabel('Amount of salt (a.u.)')

pl.legend(['Plain SA', 'K-D accel' , 'K-D accel + P-R average'], loc='best')
pl.show()
