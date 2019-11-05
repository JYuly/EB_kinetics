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
beta = 38.91 #1/kT in 1/eV
V0 = 0.1
reservoir_rate = 10**7 #This is the default rate constant for electron flow INTO the reservoirs. This can be changed directly using the last argument of the Add_reservoir() function


slope = 0.15 #this is the "slope" of the energy landscapes (i.e. the difference in reduction potentials of neighboring cofactors)

cofactor_distance = 10 #This is the distance between neighboring cofactors
    
N = 100 #The number of points to be plotted
res2emin = -0.1 #The range of energies of the 2-electron (D) reservoir to be plotted over
res2emax = 0.1
dx = (res2emax-res2emin)/N #energy step size
data  = [] #initiate arrays to store data
data2 = []
data3 = []
for n in range(N):
    #initialize the cofactor data dictionary IMPORTANT: THIS DICTIONARY MUST BE INITIALIZED
    #LIKE THIS,to satisfy the way the rate matrix is represented internally (see the companion fncts.py for more details)
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
    Add_connection(cofactor_data, "Red", "L1", 13)#this is a short-circuit (enhanced by shortening distance by 5 Angstrom to show robustness)

    Add_connection(cofactor_data, "Red", "L2", cofactor_distance*2)#this is a short-circuit
    Add_connection(cofactor_data, "Red", "H2", cofactor_distance*2)
    
    #Connect the semiquinone (SQ) form of the bifurcating cofactor to the proximal acceptors
    Add_connection(cofactor_data, "SQ", "L1", 13)
    Add_connection(cofactor_data, "SQ", "H1", cofactor_distance) #this is a short-circuit(enhanced by shortening distance by 5 Angstrom to show robustness)

    Add_connection(cofactor_data, "SQ", "L2", cofactor_distance*2)
    Add_connection(cofactor_data, "SQ", "H2", cofactor_distance*2) #this is a short-circuit

    
    #Connect cofactors further down the branches
    Add_connection(cofactor_data, "H1", "H2", cofactor_distance)
    Add_connection(cofactor_data, "L1", "L2", cofactor_distance)

    #Add short-circuit rates directly between branches
    Add_connection(cofactor_data, "H1", "L1", cofactor_distance*2)#this is a short-circuit
    Add_connection(cofactor_data, "H2", "L1", cofactor_distance*3)#this is a short-circuit
    Add_connection(cofactor_data, "H1", "L2", cofactor_distance*3)#this is a short-circuit
    Add_connection(cofactor_data, "H2", "L2", cofactor_distance*4)#this is a short-circuit

    #Connect bifurcating and terminal cofactors to reservoirs
    Add_reservoir(cofactor_data, "Red", 2, "Two-electron reservoir", res2emin + dx*n ,reservoir_rate) 
    Add_reservoir(cofactor_data, "H2", 1, "High potential reservoir", 0 ,reservoir_rate)
    Add_reservoir(cofactor_data, "L2", 1, "Low potential reservoir", 0 ,reservoir_rate)  

    K = Construct_rate_matrix(cofactor_data, reorganization_E, V0, beta)


    pop_init = np.zeros(len(K))
    pop_init[0]=1
    T = 20
    
    data.append(Flux_reservoir(cofactor_data, "Red", Evolve(K,T,pop_init), beta))
    data2.append(Flux_reservoir(cofactor_data, "H2", Evolve(K,T,pop_init), beta))
    data3.append(Flux_reservoir(cofactor_data, "L2", Evolve(K,T,pop_init), beta))
    
x = np.linspace(res2emin*1000, res2emax*1000, N) #*1000 to convert to meV

plt.rc('font', family='DejaVu Sans')
plt.rc('xtick', labelsize='large')
plt.rc('ytick', labelsize='large')
plt.rc('text', usetex=True)

fig = plt.figure(figsize=(4, 3))
ax = fig.add_subplot(1, 1, 1)
ax.ticklabel_format(style='sci',scilimits=(-3,4),axis='y')
plt.plot(x, data, '#91009B', linestyle='--')
plt.plot(x, data2,'#2D00C8')
plt.plot(x, data3,'#B40005')
ax.set_xlabel('$\Delta G_{bifurc}$ (meV)',size='x-large')
ax.set_ylabel('Flux (Sec$^{-1}$)',size='x-large') 
plt.gcf().subplots_adjust(bottom=0.2)
plt.gcf().subplots_adjust(left=0.2)

plt.savefig("S1E.svg")



