"""
Created By Manikandan Nishanth M
ME21B105
"""
############################################----IMPORT----#######################################
from time import * # To calculate the execution time
import numpy as np # Array operations
import matplotlib.pyplot as plt # Plotting

############################################----SOLVER----#######################################
class Laplace:
    """Package to solve Laplace Equation over recangular domain
    
    Parameters
    ----------
    l_x
        Length along x-axis
    l_y
        Length along y-axis
    n_x
        Number of nodes along x-axis
    n_y
        Number of nodes along y-axis
    T_0
        Initial temperature
    """
    def __init__(self,l_x,l_y,n_x,n_y,T_0) -> None:
        self.Lx = l_x # Length of plate along x-axis
        self.Ly = l_y # Length of plate along y-axis
        self.Nx = n_x # Number of grid points along x-axis
        self.Ny = n_y # Number of grid points along y-axis
        self.T0 = T_0 # Initial temperature of plate (Not required actually)

    def generateMesh(self):
        """Create mesh for FDM for given mesh
        information"""
        # self.T_k = self.T0*np.zeros((self.Nx,self.Ny)) # initially at zero temperature
        self.T_k = self.T0*np.ones((self.Ny,self.Nx)) # initially at temperature T0 
        # Use any one initial condition of two, comment out the other
        self.del_x = self.Lx/(self.Nx-1) # Node spacing
        self.del_y = self.Ly/(self.Ny-1) # Node spacing
        self.beta = self.del_x/self.del_y # Aspect ratio
        return None

    def applyBoundaryConditions(self,Tu,Tl,Tr,Tb):
        """Imposes boundary conditions on the rectangular mesh
        
        Parameters
        ----------
        Tu
            Temperature of Upper face
        Tl
            Temperature of Left face
        Tr
            Temperature of Right face
        Tb
            Temperature of Bottom face
        """
        self.T_k[:,0] *= 0
        self.T_k[:,0] += Tl
        self.T_k[:,-1] *= 0
        self.T_k[:,-1] += Tr
        self.T_k[0,:] *= 0
        self.T_k[0,:] += Tu
        self.T_k[-1,:] *= 0
        self.T_k[-1,:] += Tb

    def error(self,T1,T2):
        """Calculates error between temperatures of two consequent iterations
        
        Parameter
        ---------
        T1
            Temperature of every nodes in mesh
        T2
            Temperature of every node in mesh of another iteration
        
        Returns
        -------
        Error between temperature of two meshes"""
        return np.sum(np.abs(T1-T2))
    
    def plotPlate(self,plate,ax,iterations,scheme):
        """Plots contour, but in a better way
        
        Parameters
        ----------
        plate
            The value of the temperatures at every node in meshgrid
        ax
            axis object which is returned from matplotlib.pyplot.subplots
        iterations
            The number of iterations needed to create this distribution
        scheme
            Scheme executed"""
        pl = ax.contourf(plate,levels=range(0,101,2),cmap=plt.cm.jet,vmax=100,vmin=0)
        plt.colorbar(pl,ax=ax)
        ax.set_title(f"Temperature Distribution with {scheme},converged in {iterations} iters")
        ax1 = plt.gca()
        ax1.set_ylim(ax.get_ylim()[::-1])

    def pointGaussSiedel(self):
        """Executes point Gauss Siedel method to solve for mesh"""
        Tk = self.T_k.copy() # Keeps on updating for every iteration
        Tk1 = self.T_k.copy() # Holds new values of temperature
        Tk0 = self.T_k.copy() # Holds previous values of temperature
        # print(Tk.shape)
        iterations = 0 # iteration counter
        fig,ax = plt.subplots() # Figure
        start = time()
        while(True):
            # Until convergence, calculation goes on
            Tk0 = Tk.copy() # Holding old value for error calculation
            for i in range(self.Ny-2,0,-1):
                for j in range(1,self.Nx-1,1):
                    # Moving bottom left corner to right top corner
                    Tk1[i,j] = (Tk1[i,j-1] + Tk[i,j+1] + (self.beta**2)*(Tk1[i-1,j] + Tk[i+1,j]))/(2*(1+(self.beta**2))) # Equation (4)
            Tk = Tk1.copy() # Updating the temperature
            iterations += 1 # Updating the counter
            if self.error(Tk0,Tk) < 0.01:
                # Convergence check
                break
        end = time()
        print(end-start)
        # Plotting the result
        self.plotPlate(Tk,ax,iterations,"PGSM")
        fig.savefig(f"PGSM_{iterations}.svg")
        fig.savefig(f"PGSM_{iterations}.png")
        plt.show()
        return Tk,iterations,self.error(Tk0,Tk),(end-start)
    
    def pSOR(self,w):
        """Executes Point Successive Over-relaxation method"""
        Tk = self.T_k.copy() # Keeps on updating for every iteration
        Tk1 = self.T_k.copy() # Holds new values of temperature
        Tk0 = self.T_k.copy() # Holds previous values of temperature
        # print(Tk.shape)
        iterations = 0 # iteration counter
        fig,ax = plt.subplots() # Figure
        start = time()
        while(True):
            # Until convergence, calculation goes on
            Tk0 = Tk.copy() # Holding old value for error calculation
            for i in range(self.Ny-2,0,-1):
                for j in range(1,self.Nx-1,1):
                    # Moving bottom left corner to right top corner
                    Tk1[i,j] = (1-w)*Tk[i,j] + w*(Tk1[i,j-1] + Tk[i,j+1] + (self.beta**2)*(Tk1[i-1,j] + Tk[i+1,j]))/(2*(1+(self.beta**2))) # Equation (5)
            Tk = Tk1.copy() # Updating the temperature
            iterations += 1 # Updating the counter
            if self.error(Tk0,Tk) < 0.01:
                # Convergence check
                break
        end = time()
        print(end-start)
        # Plotting the result
        self.plotPlate(Tk,ax,iterations,"PSOR")
        fig.savefig(f"PSOR_{iterations}.svg")
        fig.savefig(f"PSOR_{iterations}.png")
        plt.show()
        return Tk,iterations,self.error(Tk0,Tk),(end-start)

    def lineGaussSiedel(self):
        """Executes Line Gauss Siedel method"""
        Tk = self.T_k.copy()
        Tk1 = self.T_k.copy()
        Tk0 = self.T_k.copy()
        iterations = 0
        fig,ax = plt.subplots()
        start = time()
        while(True):
            # Iterate until error converge
            # Alternate approach - collect the tri-diagonal elements and solve the problem
            Tk0 = Tk.copy() # Take a copy of previous iteration
            for i in range(self.Ny-2,0,-1):
                # Solve for each row
                A = np.zeros((self.Nx-2,self.Nx-2)) # Tri-diagonal matrix
                F1 = np.zeros((self.Nx-2,1)) # Vector of known values
                print(f"Iter:{iterations}, row:{i}")
                for j in range(1,self.Nx-1,1):
                    # Assembling the matrix for this row by assembling the matrix
                    if j==1:
                        # For the left-most node 
                        a = np.array([2*(1+self.beta**2),-1])
                        f1 = (self.beta**2)*(Tk1[i-1,j] + Tk[i+1,j]) + Tk1[i,j-1]
                        # Assemble
                        A[j-1,j-1:j+1] += a
                        F1[j-1,0] += f1
                    elif j==self.Nx-2:
                        # For right-most node
                        a = np.array([-1,2*(1+self.beta**2)])
                        f1 = (self.beta**2)*(Tk1[i-1,j] + Tk[i+1,j]) + Tk1[i,j+1]
                        # Assemble
                        A[j-1,j-2:j] += a
                        F1[j-1,0] += f1
                    else:
                        # For other nodes in the row
                        a = np.array([-1,2*(1+self.beta**2),-1])
                        f1 = (self.beta**2)*(Tk1[i-1,j] + Tk[i+1,j])
                        # Assemble
                        A[j-1,j-2:j+1] += a
                        F1[j-1,0] += f1
                # As the loop is over, add to the row
                Tk1[i,1:self.Nx-1] = np.transpose(np.matmul(np.linalg.inv(A),F1))
            # When the whole mesh is computed for one iteration, check for error
            Tk = Tk1.copy()
            iterations += 1 # Update the counter
            if self.error(Tk0,Tk) < 0.01:
                break
        end = time()
        print(end-start)
        # Plot the error
        self.plotPlate(Tk,ax,iterations,"LGSM")
        fig.savefig(f"LGSM_{iterations}.svg")
        fig.savefig(f"LGSM_{iterations}.png")
        plt.show()
        return Tk,iterations,self.error(Tk0,Tk),(end-start)

    def lSOR(self,w):
        """Executes Line Successive Over-relaxation method"""
        Tk = self.T_k.copy()
        Tk1 = self.T_k.copy()
        Tk0 = self.T_k.copy()
        iterations = 0
        fig,ax = plt.subplots()
        start = time()
        while(True):
            # Iterate until error converge
            # Alternate approach - collect the tri-diagonal elements and solve the problem
            Tk0 = Tk.copy() # Take a copy of previous iteration
            for i in range(self.Ny-2,0,-1):
                # Solve for each row
                A = np.zeros((self.Nx-2,self.Nx-2)) # Tri-diagonal matrix
                F1 = np.zeros((self.Nx-2,1)) # Vector of known values
                print(f"Iter:{iterations}, row:{i}")
                for j in range(1,self.Nx-1,1):
                    # Assembling the matrix for this row by assembling the matrix
                    if j==1:
                        # For the left-most node 
                        a = np.array([2*(1+self.beta**2),-w])
                        f1 = w*(self.beta**2)*(Tk1[i-1,j] + Tk[i+1,j]) + w*Tk1[i,j-1] + 2*(1-w)*(1+self.beta**2)*Tk[i,j]
                        # Assemble
                        A[j-1,j-1:j+1] += a
                        F1[j-1,0] += f1
                    elif j==self.Nx-2:
                        # For right-most node
                        a = np.array([-w,2*(1+self.beta**2)])
                        f1 = w*(self.beta**2)*(Tk1[i-1,j] + Tk[i+1,j]) + w*Tk1[i,j+1] + 2*(1-w)*(1+self.beta**2)*Tk[i,j]
                        # Assemble
                        A[j-1,j-2:j] += a
                        F1[j-1,0] += f1
                    else:
                        # For other nodes in the row
                        a = np.array([-w,2*(1+self.beta**2),-w])
                        f1 = w*(self.beta**2)*(Tk1[i-1,j] + Tk[i+1,j]) + 2*(1-w)*(1+self.beta**2)*Tk[i,j]
                        # Assemble
                        A[j-1,j-2:j+1] += a
                        F1[j-1,0] += f1
                # As the loop is over, add to the row
                Tk1[i,1:self.Nx-1] = np.transpose(np.matmul(np.linalg.inv(A),F1))
            # When the whole mesh is computed for one iteration, check for error
            Tk = Tk1.copy()
            iterations += 1 # Update the counter
            print(self.error(Tk0,Tk))
            if self.error(Tk0,Tk) < 0.01:
                break
        end = time()
        print(end-start)
        # Plot the error
        self.plotPlate(Tk,ax,iterations,"LSOR")
        fig.savefig(f"LSOR_{iterations}.svg")
        fig.savefig(f"LSOR_{iterations}.png")
        plt.show()
        return Tk,iterations,self.error(Tk0,Tk),(end-start)

    def adi(self):
        """Executes Alternate Direction Implicit method"""
        Tk = self.T_k.copy()
        Tk1 = self.T_k.copy()
        Tk0 = self.T_k.copy()
        iterations = 0
        fig,ax = plt.subplots()
        start = time()
        while(True):
            # Iterate until error converge
            # Alternate approach - collect the tri-diagonal elements and solve the problem
            Tk0 = Tk.copy() # Take a copy of previous iteration
            for i in range(self.Ny-2,0,-1):
                # Solve for each row
                A = np.zeros((self.Nx-2,self.Nx-2)) # Tri-diagonal matrix
                F1 = np.zeros((self.Nx-2,1)) # Vector of known values
                print(f"Iter:{iterations}, row:{i}")
                for j in range(1,self.Nx-1,1):
                    # Assembling the matrix for this row by assembling the matrix
                    if j==1:
                        # For the left-most node 
                        a = np.array([2*(1+self.beta**2),-1])
                        f1 = (self.beta**2)*(Tk1[i-1,j] + Tk[i+1,j]) + Tk1[i,j-1]
                        # Assemble
                        A[j-1,j-1:j+1] += a
                        F1[j-1,0] += f1
                    elif j==self.Nx-2:
                        # For right-most node
                        a = np.array([-1,2*(1+self.beta**2)])
                        f1 = (self.beta**2)*(Tk1[i-1,j] + Tk[i+1,j]) + Tk1[i,j+1]
                        # Assemble
                        A[j-1,j-2:j] += a
                        F1[j-1,0] += f1
                    else:
                        # For other nodes in the row
                        a = np.array([-1,2*(1+self.beta**2),-1])
                        f1 = (self.beta**2)*(Tk1[i-1,j] + Tk[i+1,j])
                        # Assemble
                        A[j-1,j-2:j+1] += a
                        F1[j-1,0] += f1
                # As the loop is over, add to the row
                Tk1[i,1:self.Nx-1] = np.transpose(np.matmul(np.linalg.inv(A),F1))
            Tk = Tk1.copy()
            for j in range(self.Nx-2,0,-1):
                # Solve for each column
                A = np.zeros((self.Ny-2,self.Ny-2)) # Tri-diagonal matrix
                F1 = np.zeros((self.Ny-2,1)) # Vector of known values
                print(f"Iter:{iterations}, column:{i}")
                for i in range(1,self.Ny-1,1):
                    # Assembling the matrix for this column by assembling the matrix
                    b = self.beta**2
                    if i==1:
                        # For the top-most node 
                        a = np.array([2*(1+self.beta**2),-b])
                        f1 = (Tk1[i,j-1] + Tk[i,j+1]) + b*Tk1[i-1,j]
                        # Assemble
                        A[i-1,i-1:i+1] += a
                        F1[i-1,0] += f1
                    elif i==self.Ny-2:
                        # For bottom-most node
                        a = np.array([-b,2*(1+self.beta**2)])
                        f1 = (Tk1[i,j-1] + Tk[i,j+1]) + b*Tk1[i+1,j]
                        # Assemble
                        A[i-1,i-2:i] += a
                        F1[i-1,0] += f1
                    else:
                        # For other nodes in the column
                        a = np.array([-b,2*(1+self.beta**2),-b])
                        f1 = (Tk1[i,j-1] + Tk[i,j+1])
                        # Assemble
                        A[i-1,i-2:i+1] += a
                        F1[i-1,0] += f1
                # As the loop is over, add to the column
                Tk1[1:self.Ny-1,j] = np.transpose(np.matmul(np.linalg.inv(A),F1))
            Tk = Tk1.copy()
            iterations += 1
            if self.error(Tk0,Tk) < 0.01:
                break
        end = time()
        print(end-start)
        # Plot the error
        self.plotPlate(Tk,ax,iterations,"ADI")
        fig.savefig(f"ADI_{iterations}.svg")
        fig.savefig(f"ADI_{iterations}.png")
        plt.show()
        return Tk,iterations,self.error(Tk0,Tk),(end-start)

    def adiSOR(self,w):
        """Executes ADI with Successive Over-relaxation"""
        Tk = self.T_k.copy()
        Tk1 = self.T_k.copy()
        Tk0 = self.T_k.copy()
        iterations = 0
        fig,ax = plt.subplots()
        start = time()
        while(True):
            # Iterate until error converge
            # Alternate approach - collect the tri-diagonal elements and solve the problem
            Tk0 = Tk.copy() # Take a copy of previous iteration
            for i in range(self.Ny-2,0,-1):
                # Solve for each row
                A = np.zeros((self.Nx-2,self.Nx-2)) # Tri-diagonal matrix
                F1 = np.zeros((self.Nx-2,1)) # Vector of known values
                print(f"Iter:{iterations}, row:{i}")
                for j in range(1,self.Nx-1,1):
                    # Assembling the matrix for this row by assembling the matrix
                    if j==1:
                        # For the left-most node 
                        a = np.array([2*(1+self.beta**2),-w])
                        f1 = w*(self.beta**2)*(Tk1[i-1,j] + Tk[i+1,j]) + w*Tk1[i,j-1] + 2*(1-w)*(1+self.beta**2)*Tk[i,j]
                        # Assemble
                        A[j-1,j-1:j+1] += a
                        F1[j-1,0] += f1
                    elif j==self.Nx-2:
                        # For right-most node
                        a = np.array([-w,2*(1+self.beta**2)])
                        f1 = w*(self.beta**2)*(Tk1[i-1,j] + Tk[i+1,j]) + w*Tk1[i,j+1] + 2*(1-w)*(1+self.beta**2)*Tk[i,j]
                        # Assemble
                        A[j-1,j-2:j] += a
                        F1[j-1,0] += f1
                    else:
                        # For other nodes in the row
                        a = np.array([-w,2*(1+self.beta**2),-w])
                        f1 = w*(self.beta**2)*(Tk1[i-1,j] + Tk[i+1,j]) + 2*(1-w)*(1+self.beta**2)*Tk[i,j]
                        # Assemble
                        A[j-1,j-2:j+1] += a
                        F1[j-1,0] += f1
                # As the loop is over, add to the row
                Tk1[i,1:self.Nx-1] = np.transpose(np.matmul(np.linalg.inv(A),F1))
            Tk = Tk1.copy()
            for j in range(1,self.Nx-1,1):
                # Solve for each column
                A = np.zeros((self.Ny-2,self.Ny-2)) # Tri-diagonal matrix
                F1 = np.zeros((self.Ny-2,1)) # Vector of known values
                print(f"Iter:{iterations}, column:{j}")
                for i in range(1,self.Ny-1,1):
                    # Assembling the matrix for this column by assembling the matrix
                    b = self.beta**2
                    if i==1:
                        # For the top-most node 
                        a = np.array([2*(1+self.beta**2),-w*b])
                        f1 = w*(Tk1[i,j-1] + Tk[i,j+1]) + w*b*Tk1[i-1,j] + 2*(1-w)*(1+self.beta**2)*Tk[i,j]
                        # Assemble
                        A[i-1,i-1:i+1] += a
                        F1[i-1,0] += f1
                    elif i==self.Ny-2:
                        # For bottom-most node
                        a = np.array([-w*b,2*(1+self.beta**2)])
                        f1 = w*(Tk1[i,j-1] + Tk[i,j+1]) + w*b*Tk1[i+1,j] + 2*(1-w)*(1+self.beta**2)*Tk[i,j]
                        # Assemble
                        A[i-1,i-2:i] += a
                        F1[i-1,0] += f1
                    else:
                        # For other nodes in the column
                        a = np.array([-w*b,2*(1+self.beta**2),-w*b])
                        f1 = w*(Tk1[i,j-1] + Tk[i,j+1]) + 2*(1-w)*(1+self.beta**2)*Tk[i,j]
                        # Assemble
                        A[i-1,i-2:i+1] += a
                        F1[i-1,0] += f1
                # As the loop is over, add to the column
                Tk1[1:self.Ny-1,j] = np.transpose(np.matmul(np.linalg.inv(A),F1))
            Tk = Tk1.copy()
            iterations += 1
            if self.error(Tk0,Tk) < 0.01:
                break
        end = time()
        print(end-start)
        # Plot the error
        self.plotPlate(Tk,ax,iterations,"ADISOR")
        fig.savefig(f"ADISOR_{iterations}.svg")
        fig.savefig(f"ADISOR_{iterations}.png")
        plt.show()
        return Tk,iterations,self.error(Tk0,Tk),(end-start)

##########################################----END OF PROGRAM----##################################
