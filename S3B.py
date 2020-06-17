import numpy as np
from scipy import linalg
import math
#packages for plotting
import matplotlib               
import matplotlib.pyplot as plt
#packages for interactive widgets 
from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widget

from fncts import * #imports functions specific to this electron bifurcation kinetics model



reorganization_E = 0.7
beta = 38.91 #1/kT in 1/eV
V0 = 0.1
reservoir_rate =  10**6

#initialize the cofactor data dictionary IMPORTANT: THIS DICTIONARY MUST BE INITIALIZED
#EXACTLY LIKE THIS,to satisfy the way the rate matrix is represented internally
N = 100
res2emax = 0.1
res2emin = -0.1
dx = (res2emax-res2emin)/N
data = []
data2 = []
data3 = []
data4 = []


for n in range(N):
    cofactor_data = {"Red":{"Cofactor ID": 0, "Redox State": 2, "Oxidation Potential": 0.5, "Acceptors":{}},
                 "SQ":{"Cofactor ID": 0, "Redox State": 1, "Oxidation Potential": -0.35, "Acceptors":{}},
                "Ox":{"Cofactor ID": 0, "Redox State": 0, "Oxidation Potential": "N/A", "Acceptors":{}}}



    Add_cofactor(cofactor_data, "H1", 0.285)
    Add_cofactor(cofactor_data, "H2", 0.2)
    Add_cofactor(cofactor_data, "L1", -0.088)
    Add_cofactor(cofactor_data, "L2", 0.050)

    Add_connection(cofactor_data, "Red", "L1", 12.4)
    Add_connection(cofactor_data, "Red", "H1", 16.7) #Artificially increase distance since this is known to be slower
    Add_connection(cofactor_data, "H1", "H2", 10)
    #Add_connection(cofactor_data, "L1", "L2", 12.3) #Remove rate due to inhibitor/heme knockout
    Add_connection(cofactor_data, "L1", "H1", 7+12.4) #leak rate directly from L to H
    

    Add_connection(cofactor_data, "SQ", "L1", 12.4)
    Add_connection(cofactor_data, "SQ", "H1", 7)
    

    Add_reservoir(cofactor_data, "Red", 2, "Two-electron reservoir", res2emin + dx*n, 10**7) 
    Add_reservoir(cofactor_data, "H2", 1, "High potential reservoir", 0.0, 10**7)
    Add_reservoir(cofactor_data, "L2", 1, "Low potential reservoir", -0.040,10**7)  

    K = Construct_rate_matrix(cofactor_data, reorganization_E, V0, beta)


    pop_init = np.zeros(len(K))
    pop_init[0]=1
    T = 100
    data.append(Flux_reservoir(cofactor_data, "Red", Evolve(K,T,pop_init), beta))
    data2.append(Flux_reservoir(cofactor_data, "H2", Evolve(K,T,pop_init), beta))
    data3.append(Flux_reservoir(cofactor_data, "L2", Evolve(K,T,pop_init), beta))
        
x = np.linspace((res2emin)*1000, (res2emax)*1000, N) #*1000 to convert to meV

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
ax.set_xlabel('$\Delta G_{Q_o}$ (meV)',size='x-large')
ax.set_ylabel('Flux (Sec$^{-1}$)',size='x-large') 
plt.gcf().subplots_adjust(bottom=0.2)
plt.gcf().subplots_adjust(left=0.2)

plt.savefig("S3B.svg")

