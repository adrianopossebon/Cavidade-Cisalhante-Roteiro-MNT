# In this file we have the functions used to calculate the pressure.


# Libraries.
from numba import njit


# Function to calculate the pressure.
@njit
def calculate_pressure(p, iterations, u_star, v_star, Nx, Ny, dx, dy, dt, omega, tol):
    
    
    # Initial value of the error. 
    error_max = 10.0
    
    # Number of iterations in SOR. 
    iterations = 0

    while error_max >= tol:

        R_max = 0.0 

        for i in range(0,Nx):
            for j in range(0,Ny):

                # Left bottom.
                if i == 0 and j == 0:
                    a_ii = -1.0/(dx*dx) - 1.0/(dy*dy) 
                    b = (u_star[i+1,j] - u_star[i,j])/(dx*dt) + (v_star[i,j+1] - 
                       v_star[i,j])/(dy*dt)
                    R = (p[i+1,j] - 1.0*p[i,j])/(dx*dx) 
                    R = R + (p[i,j+1] - 1.0*p[i,j])/(dy*dy)

                # Left top.
                elif i == 0 and j == Ny-1:
                    a_ii = -1.0/(dx*dx) - 1.0/(dy*dy) 
                    b = (u_star[i+1,j] - u_star[i,j])/dx + (v_star[i,j+1] - 
                       v_star[i,j])/dy
                    b = b/dt
                    R = (p[i+1,j] - 1.0*p[i,j])/(dx*dx)
                    R = R + (-1.0*p[i,j] + p[i,j-1])/(dy*dy)
                
                # Right bottom.
                elif i == Nx-1 and j == 0:
                    a_ii = -1.0/(dx*dx) - 1.0/(dy*dy) 
                    b = (u_star[i+1,j] - u_star[i,j])/dx + (v_star[i,j+1] - 
                       v_star[i,j])/dy
                    b = b/dt
                    R = (- 1.0*p[i,j] + p[i-1,j])/(dx*dx)
                    R = R + (p[i,j+1] - 1.0*p[i,j])/(dy*dy)
                
                # Right top.
                elif i == Nx-1 and j == Ny-1:
                    a_ii = -1.0/(dx*dx) - 1.0/(dy*dy) 
                    b = (u_star[i+1,j] - u_star[i,j])/dx + (v_star[i,j+1] - 
                       v_star[i,j])/dy
                    b = b/dt
                    R = (- 1.0*p[i,j] + p[i-1,j])/(dx*dx)
                    R = R + (- 1.0*p[i,j] + p[i,j-1])/(dy*dy)

                # Left wall. 
                elif i == 0:
                    a_ii = -1.0/(dx*dx) - 2.0/(dy*dy) 
                    b = (u_star[i+1,j] - u_star[i,j])/dx + (v_star[i,j+1] - 
                       v_star[i,j])/dy
                    b = b/dt
                    R = (p[i+1,j] - 1.0*p[i,j])/(dx*dx)
                    R = R + (p[i,j+1] - 2.0*p[i,j] + p[i,j-1])/(dy*dy)

                # Right wall. 
                elif i == Nx-1:
                    a_ii = -1.0/(dx*dx) - 2.0/(dy*dy) 
                    b = (u_star[i+1,j] - u_star[i,j])/dx + (v_star[i,j+1] - 
                       v_star[i,j])/dy
                    b = b/dt
                    R = (- 1.0*p[i,j] + p[i-1,j])/(dx*dx)
                    R = R + (p[i,j+1] - 2.0*p[i,j] + p[i,j-1])/(dy*dy)

                # Bottom wall. 
                elif j == 0:
                    a_ii = -2.0/(dx*dx) - 1.0/(dy*dy) 
                    b = (u_star[i+1,j] - u_star[i,j])/dx + (v_star[i,j+1] - 
                       v_star[i,j])/dy
                    b = b/dt
                    R = (p[i+1,j] - 2.0*p[i,j] + p[i-1,j])/(dx*dx)
                    R = R + (p[i,j+1] - 1.0*p[i,j])/(dy*dy)

                # Top wall. 
                elif j == Ny-1:
                    a_ii = -2.0/(dx*dx) - 1.0/(dy*dy) 
                    b = (u_star[i+1,j] - u_star[i,j])/dx + (v_star[i,j+1] - 
                       v_star[i,j])/dy
                    b = b/dt
                    R = (p[i+1,j] - 2.0*p[i,j] + p[i-1,j])/(dx*dx)
                    R = R + (- 1.0*p[i,j] + p[i,j-1])/(dy*dy)

                # Internal points. 
                else:
                    a_ii = -2.0/(dx*dx) - 2.0/(dy*dy) 
                    b = (u_star[i+1,j] - u_star[i,j])/dx + (v_star[i,j+1] - 
                       v_star[i,j])/dy
                    b = b/dt
                    R = (p[i+1,j] - 2.0*p[i,j] + p[i-1,j])/(dx*dx)
                    R = R + (p[i,j+1] - 2.0*p[i,j] + p[i,j-1])/(dy*dy)
                  
                R = (1.0/a_ii)*(b - R)
                p[i,j] = p[i,j] + omega*R

                if abs(R) > R_max: 
                    R_max = abs(R)

        iterations = iterations + 1
        error_max = R_max

    # Updating the Ghost Points for Pressure.
    for i in range(0, Nx):
        p[i,-1] = p[i,0]
        p[i,Ny] = p[i,Ny-1]
    
    for j in range(0, Ny):
        p[-1,j] = p[0,j]
        p[Nx,j] = p[Nx-1,j]

    p[-1,-1] = p[0,0]
    p[-1,Ny] = p[0,Ny-1]
    p[Nx,-1] = p[Nx-1,0]
    p[Nx,Ny] = p[Nx-1,Ny-1]

    return p, iterations
