# In this file we have the functions used to calculate the velocities, 
# both star and new. 


# Libraries.
from numba import njit


# Function to calculate u_star. 
# njit below is a decorator. It changes a 
# function. njit is used to make the code faster. 
# Here we are using the numba library. 
@njit    
def calculate_u_star(u_star, u, v, u_top, Re, Nx, Ny, dx, dy, dt):

    # Calculation of u_star for the internal points.
    for i in range(1, Nx):
        for j in range(0, Ny):
            # Interpolation of v.
            v_interp = 0.25*(v[i,j+1] + v[i,j] + v[i-1,j+1] + v[i-1,j])
            # Advection term. 
            adv = u[i,j]*(u[i+1,j] - u[i-1,j])/(2.0*dx)  
            adv += v_interp*(u[i,j+1] - u[i,j-1])/(2.0*dy)
            # Viscous term. 
            visc = (u[i+1,j] - 2.0*u[i,j] + u[i-1,j])/(Re*dx*dx)
            visc += (u[i,j+1] - 2.0*u[i,j] + u[i,j-1])/(Re*dy*dy)
            # u_star. 
            u_star[i,j] = u[i,j] + dt*(-adv + visc)
      
    # Boundary update.  
    for i in range(1,Nx):
        u_star[i,-1] = -u_star[i,0]
        u_star[i,Ny] = 2.0*u_top[i] - u_star[i,Ny-1]

    return u_star


# Function to calculate v_star.
@njit
def calculate_v_star(v_star, u, v, Re, Nx, Ny, dx, dy, dt):

    # Calculation of v_star for the internal points.
    for i in range(0, Nx):
        for j in range(1, Ny):
            # Interpolation of u.
            u_interp = 0.25*(u[i+1,j] + u[i,j] + u[i+1,j-1] + u[i,j-1])
            # Advection term. 
            adv = u_interp*(v[i+1,j] - v[i-1,j])/(2.0*dx) 
            adv += v[i,j]*(v[i,j+1] - v[i,j-1])/(2.0*dy)
            # Viscous term.
            visc = (v[i+1,j] - 2.0*v[i,j] + v[i-1,j])/(Re*dx*dx)
            visc += (v[i,j+1] - 2.0*v[i,j] + v[i,j-1])/(Re*dy*dy)
            # v_star.
            v_star[i,j] = v[i,j] + dt*(-adv + visc)
      
    #Updating the ghost points.     
    for j in range(1, Ny):
        v_star[-1,j] = -v_star[0,j]
        v_star[Nx,j] = -v_star[Nx-1,j]

    return v_star
    
    
# Function to calculate the new velocity u.
@njit
def calculate_u_new(u, u_star, p, Nx, Ny, dx, dt):

    for i in range(1,Nx):
        for j in range(-1,Ny+1):
            u[i,j] = u_star[i,j] - dt*((p[i,j] - p[i-1,j])/dx)

    return u


# Function to calculate the new velocity v.
@njit
def calculate_v_new(v, v_star, p, Nx, Ny, dy, dt):

    for i in range(-1,Nx+1):
        for j in range(1,Ny):
            v[i,j] = v_star[i,j] - dt*((p[i,j] - p[i,j-1])/dy)

    return v
