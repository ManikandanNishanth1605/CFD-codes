"""
Created by Manikandan Nishanth
ME21B105
"""
####--------IMPORTS--------####
import numpy as np # Numpy module
import matplotlib.pyplot as plt # Matplotlib for plotting

####--------INPUTS---------####

L = 1     # Length of the rod
N_x = 11   # Number of nodes on rod

t_g = [0.1,0.5,1,2,10]      # Time stamps when rod is studied
del_t = [0.1,0.01,0.001]    # Values of time steps

alpha = 1  # Thermal diffusivity

####-----OPERATION--------####

def operation(L,N_x,del_t,t_g,alpha):
    """For the given inputs, mesh is created on rod and temperature distribution is generated.
    Plotting is also done."""
    for delt in del_t:
        # Looping through every timesteps
        print(f"For t={delt}:") # Printing values of time step
        del_x = L/(N_x-1) # Distance between nodes
        nu = alpha*delt/(del_x**2) # Value of constants
        print(f"For del_t={delt}:, nu={nu}") # Printing values of time step and constant
        T = [0]*N_x # Setting up temperatures of the nodes
        time = (np.append(np.arange(0,10,delt),10)).tolist() # All the times of the simulation
        fig,ax = plt.subplots() # Setting up plots
        plt.rcParams['text.usetex'] = True
        x = (np.arange(0,N_x)*(del_x)).tolist() # Nodes of the rod
        for t in time:
            # Looping through the times in simulation
            T1 = [0]*N_x # Temperatures at next time step
            # Boundary conditions are imposed
            T1[0] = 1
            T1[N_x-1] = 0
            for i in range(1,N_x - 1):
                # For every nodes
                T1[i] = nu*(T[i+1] + T[i-1]) + (1 - 2*nu)*T[i] # Evaluating temperature at each nodes from current time step
            T = T1 # 
            if t in t_g:
                # If the current time is present in the time we need to study 
                ax.plot(x,T1,label=f"t={t}s") # Plotting temperature at that time step
                print(f"at t={t}s, distribution is {T1}") # Printing temperatures values
        plt.grid(visible=True) # Enabling grids
        plt.xticks(x) # Xticks
        plt.yticks(np.linspace(0,1,11)) # Yticks
        plt.xlim(0,x[-1]) # Limits of X
        plt.ylim(0,1) # Limits of Y
        plt.title(f"Temperature distribution for $\\alpha={alpha}$,$\Delta t$={delt}s,$\\Delta x = {del_x}$") # Title of the plot
        plt.xlabel("Distance $(x)$ in metres") # Xlabel
        plt.ylabel("Non dimensional Temperature ($T$)") # Ylabel
        plt.legend() # Enabling legend
        fig.savefig(f"{alpha}_{del_x}m_{delt}s.png",dpi=300) # Saving the image
        plt.show() # To display the plot
    return None

####-----MAIN FUNCTION--------####
if __name__=="__main__":
    operation(L,N_x,del_t,t_g,alpha) # Running the function