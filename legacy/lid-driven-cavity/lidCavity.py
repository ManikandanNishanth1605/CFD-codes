# IMPORTS
import numpy as np                 # Numpy for array operations
import matplotlib.pyplot as plt    # for plotting results

# GIVEN PARAMETERS
Nx = 129           # Number of nodes in the x-direction
L = 1             # Domain size
Wall_Velocity = 1 # Top wall velocity
rho = 1
Re = 100          # Reynolds number
nu = Wall_Velocity*L/Re  # Kinematic viscosity
dt = 0.001        # Time step
maxIt = 50000     # Maximum iterations
maxe = 1e-7       # Maximum error for convergence

# SETUP 1D GRID
Ny = Nx           # Assuming the number of nodes are same along y-direction
h = L / (Nx - 1)  # Size of the grid
x = np.linspace(0, L, Nx) # X-coords
y = np.linspace(0, L, Ny) # Y-coords
# Arrays for access of nodes along iterations
im = np.arange(0, Nx - 2) 
i = np.arange(1, Nx - 1)
ip = np.arange(2, Nx)
jm = np.arange(0, Ny - 2)
j = np.arange(1, Ny - 1)
jp = np.arange(2, Ny)

# PREALLOCATE MATRICES
Vo = np.zeros((Nx, Ny))  # Vorticity
St = np.zeros((Nx, Ny))  # Stream function
Vop = np.zeros_like(Vo)  # Vorticity for previous iteration
u = np.zeros_like(Vo)    # Velocity in x-direction
v = np.zeros_like(Vo)    # Velocity in y-direction
p = np.zeros_like(Vo)    # Pressure

# SOLVE LOOP SIMILAR TO GAUSS-SEIDEL METHOD
for iter in range(1, maxIt + 1):
    # 1. ASSIGN BOUNDARY CONDITIONS
    Vo[:, -1] = -2 * St[:, -2] / (h**2) - 2 * Wall_Velocity / h  # Top
    Vo[:, 0] = -2 * St[:, 1] / (h**2)                            # Bottom
    Vo[0, :] = -2 * St[1, :] / (h**2)                            # Left
    Vo[-1, :] = -2 * St[-2, :] / (h**2)                          # Right

    # 2. SOLVE VORTICITY TRANSPORT EQUATION FOR VORTICITY
    Vop = Vo.copy() # Storing the old value
    Vo[i[:, None], j] = Vop[i[:, None], j] + (
        - (St[i[:, None], jp] - St[i[:, None], jm]) / (2 * h) * 
          (Vop[ip[:, None], j] - Vop[im[:, None], j]) / (2 * h) +
        (St[ip[:, None], j] - St[im[:, None], j]) / (2 * h) * 
          (Vop[i[:, None], jp] - Vop[i[:, None], jm]) / (2 * h) +
        nu * (Vop[ip[:, None], j] + Vop[im[:, None], j] - 4 * Vop[i[:, None], j] + 
                    Vop[i[:, None], jp] + Vop[i[:, None], jm]) / (h**2)
    ) * dt

    # 3. SOLVE POISSON EQUATION FOR STREAM FUNCTION USING POINT GAUSS SIEDEL
    St[i[:, None], j] = (
        Vo[i[:, None], j] * h**2 +
        St[ip[:, None], j] + St[i[:, None], jp] +
        St[i[:, None], jm] + St[im[:, None], j]
    ) / 4

    # CHECK FOR CONVERGENCE
    if iter > 10:
        error = np.max(np.abs(Vo - Vop))
        if error < maxe:
            print(f'Converged in {iter} iterations with error: {error:.2e}')
            break
else:
    print('Reached maximum iterations without convergence.')

# FIND VELOCITY AND PRESSURE FROM STREAM FUNCTION
u[1:Nx - 1, -1] = Wall_Velocity
u[i[:, None], j] = (St[i[:, None], jp] - St[i[:, None], jm]) / (2 * h)
v[i[:, None], j] = -(St[ip[:, None], j] - St[im[:, None], j]) / (2 * h)
# PRESSURE USING GAUSS ITERATIONS
for iters in range(100):
	pold = p.copy()
	rhs = 2*rho*((St[ip[:, None], j] - 2 * St[i[:, None], j] + St[im[:, None], j]) / (h**2)
			* (St[i[:, None], jp] - 2 * St[i[:, None], j] + St[i[:, None], jm]) / (h**2)
			- (1 / (4*(h**2))) * (St[ip[:, None], jp] - St[ip[:, None], jm] - St[im[:, None], jp] + St[im[:, None], jm]))
	p[i[:, None], j] = (
			rhs * h**2 +
			p[ip[:, None], j] + p[i[:, None], jp] +
			p[i[:, None], jm] + p[im[:, None], j]
	) / 4
    
	if iters > 10:
		if np.max(np.abs(p - pold)) < 1e-4:
			print(f"Pressure converged in {iters}")
			break
print(p.shape)

# PLOTS
cm = plt.cm.jet(np.linspace(0, 1, 100)[::-1])  # Colormap

# Streamlines
N = 1000
np.random.seed(42)
xstart = np.random.rand(N) * max(x)
ystart = np.random.rand(N) * max(y)
X, Y = np.meshgrid(x, y)
plt.figure(1)
plt.streamplot(X, Y, u.T, v.T, start_points=np.array([xstart, ystart]).T, color='k')
plt.title('Stream Function')
plt.xlabel('X-Location')
plt.ylabel('Y-Location')
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
plt.savefig(f"stream_{Nx}_{Re}.svg")
plt.figure(1)
arr = np.array([-1e-10,-1e-7,-1e-5,-1e-4,-0.01,-0.03,-0.05,-0.07,-0.09,-0.1])
arr = arr[::-1]
arr = np.append(arr,np.array([1e-8,1e-7,1e-6,1e-5]))
plt.contour(X, Y, St.T, arr)
plt.title('Streamlines')
plt.xlabel('X')
plt.ylabel('Y')
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
plt.savefig(f"streamcontour_{Nx}_{Re}.svg")

# Vorticity Contour
plt.figure(2)
plt.contourf(x, y, Vo.T, 23, cmap=plt.cm.jet)
plt.colorbar(label='Vorticity', orientation='horizontal')
plt.title('Vorticity Contour')
plt.xlabel('X-Location')
plt.ylabel('Y-Location')
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
plt.savefig(f"vorticity_{Nx}_{Re}.svg")

# Centerline X-Direction Velocity
plt.figure(3)
plt.plot(u[Nx // 2, :], y)
plt.title('Centerline X-Direction Velocity')
plt.xlabel('u/U')
plt.ylabel('Y/L')
plt.grid()
plt.savefig(f"u_{Nx}_{Re}.svg")

# Centerline Y-Direction Velocity
plt.figure(4)
plt.plot(x, v[:, Ny // 2])
plt.title('Centerline Y-Direction Velocity')
plt.xlabel('X/L')
plt.ylabel('v/U')
plt.grid()
plt.savefig(f"v_{Nx}_{Re}.svg")

# Pressure distribution contour
plt.figure(5)
plt.contourf(x, y, p.T, 23, cmap=plt.cm.jet)
plt.colorbar(label='Pressure', orientation='horizontal')
plt.title('Pressure Contour')
plt.xlabel('X-Location')
plt.ylabel('Y-Location')
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
plt.savefig(f"pressure_{Nx}_{Re}.svg")

# pointsx = np.array([1,8,9,10,14,23,37,59,65,80,95,110,123,124,125,126,129])
# pointsy = np.array([1,9,10,11,13,21,30,31,65,104,111,117,122,123,124,125,129])

# vel_u = u[Nx//2, pointsx-1]
# vel_v = v[pointsy-1, Ny//2]

# np.savetxt("u_pts",vel_u,delimiter=',')
# np.savetxt("v_pts",vel_v,delimiter=',')

plt.show()