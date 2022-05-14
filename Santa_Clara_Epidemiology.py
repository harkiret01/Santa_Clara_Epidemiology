# -*- coding: utf-8 -*-
"""
Created on Fri May 13 17:47:26 2022

@author: Harkiret Singh
"""

import matplotlib.pyplot as plt
from modsim import *

def make_system(alpha, beta, gamma, delta, epsilon, zeta, eta, theta, iota, kappa, labda, mu, nu, xi, omicron, phi):
    #s=30% of 1.94 million or 582000 people susceptible in Santa Clara County
    # I=2410 people infected
    init = State(S=385600, I=2410, D=0, A=0, R=0, T=0, H=0, E=0)
    init /= sum(init) #convert from number of people to fractions
    t0 = 0 #initial time
    t_end = 7 * 60 #number of days for 60 weeks (5 years)

    return System(init=init, t0=t0, t_end=t_end,alpha = alpha,beta = beta,gamma = gamma,delta = delta,epsilon = epsilon,
                  theta = theta,eta = eta,zeta = zeta,mu = mu,nu = nu,tau = tau,lamda = lamda,kappa = kappa,xi = xi,rho = rho,sigma = sigma)


contact_time = 15
recovery_time = 7
#values based on day 22 of lockdown at Italy from https://www.nature.com/articles/s41591-020-0883-7
#transmission rate
alpha = 0.36
beta = 0.005
gamma = 0.2
delta = 0.005
#rate of detection
epsilon = 0.5
theta = 0.6
#incubation rate
eta = 0.034
zeta = 0.034
#life threatening incubation rate
mu = 0.008
nu = 0.015
#mortality rate
tau = 0.03
#recovery rate
lamda = 0.08
kappa = 0.017
xi = 0.017
rho = 0.017
sigma = 0.017

system = make_system(alpha, beta, gamma, delta, epsilon, theta, eta, zeta, mu, nu, tau, lamda, kappa, xi, rho, sigma)


def update_func(state, T, system):
    constants  = ["alpha", "beta", "gamma", "delta", "epsilon", "theta", "eta", "zeta", "mu", "nu", "tau", "lamda", "kappa", "xi", "rho", "sigma"]
    (alpha, beta, gamma, delta, epsilon, theta, eta, zeta, mu, nu, tau, lamda, kappa, xi, rho, sigma) = [system[key] for key in constants]
    S,I,D,A,R,T,H,E = state
    
    infected = S*(alpha*I+beta*D+gamma*A*delta*R) - (epsilon + zeta + lamda )*I
    diagnosed = epsilon*I-(eta+rho)*D
    ailing = zeta*I-(theta+mu+kappa)*A
    recognized = eta*D+theta*A-(nu+tau)*R
    threatened = mu*A+nu*R-(sigma+tau)*T
    healed = lamda*I+rho*D+kappa*A+xi*R+sigma*T
    extinct = tau*T
    
    S -= infected
    I += infected
    D += diagnosed
    A += ailing 
    R += recognized 
    T += threatened 
    H += healed
    E += extinct
    return State(S=S,I=I,D=D,A=A,R=R,T=T,H=H,E=E)



init = State(S=385600, I=2410, D=0, A=0, R=0, T=0, H=0, E=0)
init /= sum(init)
state = update_func(init, 0, system)


def run_simulation(system, update_func):
    df = TimeFrame(columns=system.init.index)
    df.row[system.t0] = system.init

    state = system.init
    t0 = system.t0
    
    for t in linrange(system.t0, system.t_end):
        df.row[t+1] = update_func(df.row[t], t, system)
    
    return df



tc = 9      # time between contacts in days, it can be 9 days already 
tr = 20     # recovery time in days
#transmission rate * contacts per day
alpha = .004*6
beta = .02*6
gamma = .002*9
delta = .004*9
#rate of detection
epsilon = 0.02/7
theta = 0.05*0.07
#incubation rate
eta = 0.1/7 
zeta = 0.1/6
#life threatening incubation rate
mu=0.3*.02
nu=.02/7
# mortality rate
tau = .02
#recovery rate
lamda = 0.3*0.01
kappa = 0.15*0.04*0.06
xi = 0.2*0.04*0.06
rho = 0.3*0.0025
sigma = 0.04*0.002

system = make_system(alpha, beta, gamma, delta, epsilon, theta, eta, zeta, mu, nu, tau, lamda, kappa, xi, rho, sigma)
results = run_simulation(system, update_func)
results.head()

def plot_results(S, I, D, A, R, T, H,E):  
  
    
    plot(S, 'cyan', label='Susceptible')
    plot(I,'orange', label = "Infected")
    plot(D,'green', label = "Diagnosed")
    plot(A,'black', label = "Ailing")
    plot(R,'purple', label = "Recognized")
    plot(T,'orange', label = "Threatened")
    plot(H,'red', label = "Healed")
    plot(E,'blue', label = "Extinct")
    
    
    decorate(xlabel='Number of days', ylabel='Population %')                
plot_results(results.S, results.I, results.D, results.A, results.R, results.T, results.H,results.E)