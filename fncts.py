#Import necessary packages

import numpy as np
from scipy import linalg
import math
#packages for plotting
import matplotlib               
import matplotlib.pyplot as plt



# cofactor_data dictionary functions
#=========================================


#  The cofactor_data dictionary and associated functions contains provides a kind of "user interface" to work with the model parameters before the rate matrix is constructed. These functions provide a clean way to interact with this dictionary.
# For a clean view of the current contents of this dictionary, use Print_cofactor_data


#Add_cofactor adds a cofactor to the cofactor_data dictionary, along with its reduction
def Add_cofactor(cofactor_data, cofactor_name, reduction_potential):
    cofactor_ID = len(cofactor_data) - 2 #subtract 2 because Ox and ASQ are included in cofactor_data, despite being the same cofactor
    cofactor_data[cofactor_name] = {"Cofactor ID": cofactor_ID, "Redox State": 1, "Oxidation Potential": reduction_potential,"Acceptors":{}}
    return cofactor_data

# Add_connection stores the tunneling distance between a donor and acceptor cofactor. This distance is used later on to calculate tunneling rate constants, along with the difference in reduction potentials of the cofactors.
def Add_connection(cofactor_data, donor_cofactor, acceptor_cofactor, distance):
    cofactors = (donor_cofactor, acceptor_cofactor)
    if all(x in cofactor_data for x in cofactors):
        cofactor_data[donor_cofactor]["Acceptors"][acceptor_cofactor] = distance
        return cofactor_data
    raise ValueError("Cofactors not present")

# Add_reservoir() adds information about an electron reservoir to the cofactor_data dictionary: which cofactor it exchanges electrons with, how many electrons are exchanged at a time, the delta_G of the exchange, and the rate constant for electron flow INTO the reservoir.
# (WARNING: this function accepts free energy as an argument, not the reservoir reduction potential! The difference in
# reduction potential must be converted to free energy to be fed into this function)
def Add_reservoir(cofactor_data, donor_cofactor, num_electrons, reservoir_name, delta_G, in_rate):
    if donor_cofactor in cofactor_data:
        res = {"Reservoir name": reservoir_name, "Delta G":delta_G, "Num Electrons": num_electrons, "In Rate": in_rate}
        cofactor_data[donor_cofactor]["Accepting Reservoir"] = res
        return cofactor_data
    raise ValueError("Cofactor not present")

#Print_cofactor_data is a function that allows the user to see the contents of the kinetics model in a human-readable way
def Print_cofactor_data(cofactor_data):
    for cofactor_name,cofactor_info in cofactor_data.items():
        print(cofactor_name)
        print("---------")
        for lineitem,lineitemdata in cofactor_info.items():
            if (lineitem=="Acceptors"):
                print("Acceptors: \t Name \t Distance (Ang)")
                for label, data in lineitemdata.items():
                    print("\t \t", label, "\t", data)
            else:
                if (lineitem=="Accepting Reservoir"):
                    print("Accepting Reservior: \t Name \t Delta G \t Num Electrons")
                    print("\t \t \t",lineitemdata["Reservoir name"], "\t", lineitemdata["Delta G"], "\t \t", lineitemdata["Num Electrons"])
                else: 
                        print(lineitem,lineitemdata)
        print()
    return


# kinetics model functions
#=========================================

#These functions construct and manipulate the rate matrix and associated dynamics of the kinetic model. The logic of the functions works with the occupation number representation of the redox state of the bifurcase [n_B,n_L1,n_L2,...,n_H1, n_H2, ..], but this is converted to an index which is used to track the state of the system as a numpy array, so that the numpy matrix functions can be used. The kinetic model is completely scalable in the lengths of each branch

#(NOTE: there is currently no way to allow cofactors other than the electron bifurcating site to hold more than one electron. This may be implemented in a future version.)



#Construct_rate_matrix takes the cofactor_data dictionary into which the cofactor parameters are entered by the user, and contructs the rate matrix using another function Connectwdetailed balance. 

def Construct_rate_matrix(cofactor_data, reorg_E, V0, beta):
    num_cofactors = len(cofactor_data) - 2 #subtract 2 because Ox and ASQ are included in cofactor_data, despite being the same cofactor
    K = np.zeros((3*2**(num_cofactors-1), 3*2**(num_cofactors-1))) #initialize rate matrix
    

    for donor_name, donor in cofactor_data.items():
        
        donor_cofactor_ID = donor["Cofactor ID"]
        donor_redox_state = donor["Redox State"]
        donor_redox_potential = donor["Oxidation Potential"]
        
        for acceptor_name, distance in donor["Acceptors"].items():
            
            acceptor_cofactor_ID = cofactor_data[acceptor_name]["Cofactor ID"]
            acceptor_redox_state = cofactor_data[acceptor_name]["Redox State"]-1
            acceptor_redox_potential = cofactor_data[acceptor_name]["Oxidation Potential"]
            k = ETrate(-(acceptor_redox_potential-donor_redox_potential), distance, reorg_E, beta, V0)
            K = Connectwdetailedbalance(K,donor_cofactor_ID, donor_redox_state, acceptor_cofactor_ID, acceptor_redox_state, k, -(acceptor_redox_potential-donor_redox_potential), beta)
            
        if "Accepting Reservoir" in donor:
            donor_final_redox_state = donor_redox_state - donor["Accepting Reservoir"]["Num Electrons"]
            deltaG = donor["Accepting Reservoir"]["Delta G"]
            in_rate = donor["Accepting Reservoir"]["In Rate"]
            K = Reservoirwdetailedbalance(K, donor_cofactor_ID, donor_redox_state, donor_final_redox_state, in_rate, deltaG, beta)
    
    return K


#Evolve returns population vector after time t given an initial vector pop_init   
def Evolve(K,t,pop_init):
    return linalg.expm(K *t)@pop_init



#Pop finds the population of a cofactor in a given redox state. This requires "tracing out" the redox state of all other cofactors.
def Pop(P, cofactor_data, cofactor, redox_state):
    cofactor_ID = cofactor_data[cofactor]["Cofactor ID"]
    num_cofactors = int(np.log2(len(P)/3)+1)
    Pop = 0
    for i in range(len(P)):
        if (state(i, num_cofactors)[cofactor_ID] == redox_state):
           Pop = Pop +  P[i]
    return Pop




#This function plots the population of a given cofactor in the redox state found in the cofactor_data dictionary
    
def Popplot(K,cofactor_data, T,pop_init, cofactor):
    cofactor_ID = cofactor_data[cofactor]["Cofactor ID"]
    redox_state = cofactor_data[cofactor]["Redox State"]
    t = np.linspace(0,T,250)
    prob1 = [Pop(Evolve(K,t,pop_init),cofactor_ID,redox_state) for t in t]
    plt.plot(t, prob1)
    return;

#This function plots the net flux into a reservoir versus time
def Plotflux(K,cofactor_data, reservoir_rate, T,pop_init, cofactor,beta):
    t = np.linspace(0,T,250)
    flux = [Flux_reservoir(cofactor_data, reservoir_rate, cofactor, Evolve(K,t,pop_init),beta) for t in t]
    plt.plot(t, flux)
    return;

#index() returns the index (internal labeling scheme) corresponding to state = [nb, nl1, nh1, nl2, nh2,...] which is in the occupation number representation.
#state[0] is special because there are three possibilities for nb, no electrons, one electron, and two electrons
def index(state):
    index = state[0]
    for i in range(len(state)-1):
        index = index + state[i+1]*3*2**(i)
    return index


#state() returns the state = [nb, nl1, nh1, nl2, nh2,...] corresponding to the index index
#state[0] is special because there are three possibilities for nb, no electrons, one electron, and two electrons
def state(index, num_cofactors):
    state = [int(index%3)]
    num = (index - index%3)/3
    for i in range(num_cofactors-1):
        state.append(int(num%2))
        num=(num - num%2)/2     
    return state



#Connect adds rate constant k between electron donor (cofactor_i_ID) and acceptor (cofactor_f_ID), which are INITIALLY in redox_state_i and redox_state_f, respectively
#This function uses the internal labelling of the states which uses one index which maps to the occupation number representation [n_B,n_L1,n_L2,...,n_H1, n_H2, ..] and back using the index() and state() functions.
def Connect(K,cofactor_i_ID, redox_state_i, cofactor_f_ID, redox_state_f, k):
    num_cofactors = int(np.log2(len(K)/3)+1)
    
    for i in range(3*2**(num_cofactors-1)): #loop through all the states, to look for initial (donor) state
    
        if (state(i, num_cofactors)[cofactor_i_ID] == redox_state_i)and(state(i, num_cofactors)[cofactor_f_ID] == redox_state_f): #checking that the initial (donor) state has the donor and acceptor in the correct oxidation states
            
            for j in range(3*2**(num_cofactors-1)): #looping to find the final (acceptor) state    
                
                if (state(j, num_cofactors)[cofactor_i_ID] == redox_state_i -1)and(state(j, num_cofactors)[cofactor_f_ID] == redox_state_f + 1): #check that the final (acceptor) state has the correct oxidation states
                    I = np.delete(state(i, num_cofactors),[cofactor_i_ID,cofactor_f_ID]) #delete the information about the donor and acceptor oxidation states (already checked)
                    J = np.delete(state(j, num_cofactors),[cofactor_i_ID,cofactor_f_ID])
                    
                    if np.array_equal(I,J): #check that the remaining cofactors do not change state (to conserve electrons).
                        K[j][i]= K[j][i] + k #add population of final state
                        K[i][i]= K[i][i]-k #remove population at initial state
    return K



#Connectwdetailedbalance connects cofactor_i and cofactor_f with forward rate constant kf and reverse rate constant satisfying detailed balance
def Connectwdetailedbalance(K,cofactor_i, redox_state_i, cofactor_f, redox_state_f, kf, deltaG, beta):
    K = Connect(K,cofactor_i, redox_state_i, cofactor_f, redox_state_f, kf)
    K = Connect(K,cofactor_f, redox_state_f+1, cofactor_i, redox_state_i-1, kf*np.exp(beta*deltaG))
    return K



#ETrate calculates the nonadiabatic (Marcus) ET rate
def ETrate(deltaG, R, reorgE, beta, V0):
    hbar = 6.5821 * 10**(-16) #hbar in eV sec
    decay = 1 #distance decay of ET rate (in per Angstrom)
    return (2*math.pi/hbar)*(V0**2)*np.exp(-decay*R)*(1/(math.sqrt(4*math.pi*(1/beta)*reorgE)))*np.exp(-beta*(deltaG + reorgE)**2/(4*reorgE))



#Reservoir adds rate constants between redox_state_i and redox_state_f of a cofactor (specified by cofactor_ID). (WARNING: the redox state of the reservoir is never modeled.)
def Reservoir(K, cofactor_ID, redox_state_i, redox_state_f, k):
    num_cofactors = int(np.log2(len(K)/3)+1)
    for i in range(3*2**(num_cofactors-1)): #loop through all the states, to look for initial (donor) state
        if (state(i, num_cofactors)[cofactor_ID] == redox_state_i): #checking that the initial (donor) state has the donor correct oxidation states
            for j in range(3*2**(num_cofactors-1)): #looping to find the final (acceptor) state
                if (state(j, num_cofactors)[cofactor_ID] == redox_state_f): #check that the final (acceptor) state has the correct oxidation states
                    #delete the information about the donor oxidation states (already checked)
                    I = np.delete(state(i, num_cofactors),[cofactor_ID])
                    J = np.delete(state(j, num_cofactors),[cofactor_ID])
                    if np.array_equal(I,J): #check that the remaining cofactors do not change state (to conserve electrons).
                        K[i][i]= K[i][i]-k #remove population of initial state
                        K[j][i]= K[j][i]+k #add population of final state
    return K




#Reservoirwdetailedbalance calls Reservoir twice to satisfy detailed balance.
def Reservoirwdetailedbalance(K, cofactor, redox_state_i, redox_state_f, kf, deltaG, beta):
    K = Reservoir(K, cofactor, redox_state_i, redox_state_f, kf)
    K = Reservoir(K, cofactor, redox_state_f, redox_state_i, kf*np.exp(beta*deltaG))
    return K


#Flux_reservoir calculates the instantaneous net flux into the reservoir connected to cofactor.    
def Flux_reservoir(cofactor_data, cofactor, P, beta):
    cofactor_ID = cofactor_data[cofactor]["Cofactor ID"]
    in_rate =  cofactor_data[cofactor]["Accepting Reservoir"]["In Rate"]
    deltaG = cofactor_data[cofactor]["Accepting Reservoir"]["Delta G"]
    reverse_rate = in_rate*np.exp(beta*deltaG)
    donor_redox_state = cofactor_data[cofactor]["Redox State"]
    electrons_transferred = cofactor_data[cofactor]["Accepting Reservoir"]["Num Electrons"]
    final_redox_state = donor_redox_state - electrons_transferred

    
    netflux = (Pop(P, cofactor_data, cofactor, donor_redox_state)*in_rate - Pop(P,cofactor_data, cofactor, final_redox_state)*reverse_rate)*electrons_transferred
    return netflux

#Rate_cofactor() calculates the instantaneous forward rate from cofactor_i to cofactor_f
def Rate_cofactor(cofactor_data, K, cofactor_i, cofactor_f, P):
    
    num_cofactors = int(np.log2(len(P)/3)+1)
    cofactor_i_ID = cofactor_data[cofactor_i]["Cofactor ID"]
    redox_state_i = cofactor_data[cofactor_i]["Redox State"]
    cofactor_f_ID = cofactor_data[cofactor_f]["Cofactor ID"]
    redox_state_f = cofactor_data[cofactor_f]["Redox State"] - 1
    flux = 0
    for i in range(3*2**(num_cofactors-1)): #loop through all the states, to look for initial (donor) state
    
        if (state(i, num_cofactors)[cofactor_i_ID] == redox_state_i)and(state(i, num_cofactors)[cofactor_f_ID] == redox_state_f): #checking that the initial (donor) state has the donor and acceptor in the correct oxidation states
            
            for j in range(3*2**(num_cofactors-1)): #looping to find the final (acceptor) state    
                
                if (state(j, num_cofactors)[cofactor_i_ID] == redox_state_i -1)and(state(j, num_cofactors)[cofactor_f_ID] == redox_state_f + 1): #check that the final (acceptor) state has the correct oxidation states
                    I = np.delete(state(i, num_cofactors),[cofactor_i_ID,cofactor_f_ID]) #delete the information about the donor and acceptor oxidation states (already checked)
                    J = np.delete(state(j, num_cofactors),[cofactor_i_ID,cofactor_f_ID])
                    
                    if np.array_equal(I,J): #check that the remaining cofactors do not change state (to conserve electrons).
                       flux = flux + K[j][i]*P[i] #add contribution to flux
    return flux

#Flux_cofactor() calculates the instantaneous NET flux from cofactor_i to cofactor_f, by calling Rate_cofactor() twice, once for the forward rate, and once for the reverse rate
def Flux_cofactor(cofactor_data, K, cofactor_i, cofactor_f, P):
    return Rate_cofactor(cofactor_data, K, cofactor_i, cofactor_f, P) - Rate_cofactor(cofactor_data, K, cofactor_f, cofactor_i, P)

