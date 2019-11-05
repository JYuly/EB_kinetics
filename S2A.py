#standard math packages
import numpy as np
from scipy import linalg
import math

#packages for plotting
import matplotlib               
import matplotlib.pyplot as plt

from fncts import * #imports functions to use in this script


#define Marcus ET parameters
reorganization_E = 1

V0 = 0.1
reservoir_rate = 10**7


slope = 0.15

cofactor_distance = 10


boltzmann_k = 8.617333 * 10**(-5) #Boltzmann Constant in eV/K
beta_min = 5 #1/kT in 1/eV *293/373 for 100C (optimal growth condition for pyrococcus furiosus )
beta_max = 38.91*1.7
N = 90
beta_step = (beta_max - beta_min)/N

data  = []
data2 = []
data3 = []
for n in range(N):
    #initialize the cofactor data dictionary IMPORTANT: THIS DICTIONARY MUST BE INITIALIZED
    #LIKE THIS,to satisfy the way the rate matrix is represented internally
    cofactor_data = {"Red":{"Cofactor ID": 0, "Redox State": 2, "Oxidation Potential": 0.4, "Acceptors":{}},
                 "SQ":{"Cofactor ID": 0, "Redox State": 1, "Oxidation Potential": -0.4, "Acceptors":{}},
                "Ox":{"Cofactor ID": 0, "Redox State": 0, "Oxidation Potential": "N/A", "Acceptors":{}}}



    #Initiate cofactors
    Add_cofactor(cofactor_data, "H1", 0.4 - slope*1)
    Add_cofactor(cofactor_data, "H2", 0.4 - slope*2)
    Add_cofactor(cofactor_data, "L1", -0.4 + slope*1)
    Add_cofactor(cofactor_data, "L2", -0.4 + slope*2)

    #Connect the fully reduced form of the bifurcating cofactor to the proximal acceptors
    Add_connection(cofactor_data, "Red", "H1", cofactor_distance)
    Add_connection(cofactor_data, "Red", "L1", cofactor_distance) #this is a short-circuit

    Add_connection(cofactor_data, "Red", "L2", cofactor_distance*2) #this is a short-circuit
    Add_connection(cofactor_data, "Red", "H2", cofactor_distance*2)
    
    #Connect the semiquinone (SQ) form of the bifurcating cofactor to the proximal acceptors
    Add_connection(cofactor_data, "SQ", "L1", cofactor_distance)
    Add_connection(cofactor_data, "SQ", "H1", cofactor_distance) #this is a short-circuit

    Add_connection(cofactor_data, "SQ", "L2", cofactor_distance*2)
    Add_connection(cofactor_data, "SQ", "H2", cofactor_distance*2) #this is a short-circuit

    
    #Connect cofactors further down the branches
    Add_connection(cofactor_data, "H1", "H2", cofactor_distance)
    Add_connection(cofactor_data, "L1", "L2", cofactor_distance)

    #Add short-circuit rates directly between branches
    Add_connection(cofactor_data, "H1", "L1", cofactor_distance*2) #this is a short-circuit
    Add_connection(cofactor_data, "H2", "L1", cofactor_distance*3) #this is a short-circuit
    Add_connection(cofactor_data, "H1", "L2", cofactor_distance*3) #this is a short-circuit
    Add_connection(cofactor_data, "H2", "L2", cofactor_distance*4) #this is a short-circuit

    #Connect bifurcating and terminal cofactors to reservoirs
    Add_reservoir(cofactor_data, "Red", 2, "Two-electron reservoir", 0 ,reservoir_rate) 
    Add_reservoir(cofactor_data, "H2", 1, "High potential reservoir", 0 ,reservoir_rate)
    Add_reservoir(cofactor_data, "L2", 1, "Low potential reservoir", 0 ,reservoir_rate)  

    K = Construct_rate_matrix(cofactor_data, reorganization_E, V0, beta_min + n*beta_step)


    pop_init = np.zeros(len(K))
    pop_init[0]=1
    T = 20

    data.append(2*ETrate(0.15, cofactor_distance, reorganization_E, beta_min+beta_step*n, V0)* np.exp(-2*slope*(beta_min+beta_step*n))) 
    data2.append(Flux_reservoir(cofactor_data, "H2", Evolve(K,T,pop_init), beta_min + n*beta_step))
    data3.append(ETrate(-0.50, cofactor_distance*2, reorganization_E, beta_min+beta_step*n, V0)* np.exp(-2*slope*(beta_min+beta_step*n)))
x = np.linspace(beta_min, beta_max, N)

plt.rc('font', family='DejaVu Sans')
plt.rc('xtick', labelsize='large')
plt.rc('ytick', labelsize='large')
plt.rc('text', usetex=True)

fig = plt.figure(figsize=(4, 3))
ax = fig.add_subplot(1, 1, 1)
ax.ticklabel_format(style='sci',scilimits=(-3,4),axis='both')
plt.plot(x, data,  '#FF6400', linestyle='--')
plt.plot(x, data2,'#2D00C8')
plt.plot(x, data3,'#2E9100', linestyle='--')
plt.yscale('log')


ax.set_xlabel('$1/kT$ (eV$^{-1}$)',size='x-large')
ax.set_ylabel('Flux (Sec$^{-1}$)',size='x-large') 
plt.gcf().subplots_adjust(bottom=0.2)
plt.gcf().subplots_adjust(left=0.2)

plt.savefig("S2A.svg")



