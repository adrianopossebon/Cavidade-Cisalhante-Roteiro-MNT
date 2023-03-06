# In this file we have auxiliar functions that are used to export 
# the results, calculate the stream function etc. 


# Libraries.
import numpy as np
from numba import njit
import pandas as pd


# Prints on screen. 
def print_on_screen(time, time_final, pressure_iterations, div):
    """
    This function prints some informations on the
    screen during the simulation. 
    """
    
    print('\n\n')
    print(f' Time: {time: 6.4f} .')
    print(f' Final time: {time_final: 6.4f} .')
    print(f' Completed: {100.0*time/time_final: 4.1f} % .')
    print(f' Number of iterations in pressure: {pressure_iterations} .')
    print(f' Maximum value of divergence: {np.max(div): .4e} .')

    return None
    

# Checks the input parameters.
def check_parameters(Re, dx, dy, dt):
    """
    This function monitor the value of dt, dx and dy. 
    
    If the conditions for stability are not satisfied, 
    the simulation stops. 
    """

    if ((dt >= dx) or (dt >= 0.25*Re*dx*dx) or (dt >= 0.25*Re*dy*dy) 
        or (dx >= 1.0/np.sqrt(Re)) or (dy >= 1.0/np.sqrt(Re))):
        
        print('--------')
        print('--------')
        print('ATTENTION!! WRONG PARAMETERS! THE SIMULATION STOPS HERE!')
        print('--------')
        print('--------')

        if (dt >= dx):
            print('\n dt >= dx \n')
        if (dt >= 0.25*Re*dx*dx):
            print('\n dt >= 0.25*Re*dx*dx \n')
        if (dt >= 0.25*Re*dy*dy):
            print('\n dt >= 0.25*Re*dy*dy \n')
        if (dx >= 1.0/np.sqrt(Re)):
            print('\n dx >= 1.0/np.sqrt(Re)) \n')        
        if (dy >= 1.0/np.sqrt(Re)):
            print('\n dy >= 1.0/np.sqrt(Re)) \n')    
        
        quit()
        
    return None
    

# Calculates the divergence of the velocity.
@njit
def check_for_incompressibility(div, u, v, Nx, Ny, dx, dy):
    """
    This function calculates the divergence of the velocity
    in every point of the mesh. 
    """
    
    for i in range(0, Nx):
        for j in range(0, Ny):
            div[i, j] = ((u[i+1, j] - u[i, j])/dx
                + (v[i, j+1] - v[i, j])/dy)

    return div


# Vorticity.
@njit
def calculate_vorticity(vort_plot, u, v, Nx, Ny, dx, dy):
    """
    Calculates the vorticity. 
    """
    
    for i in range(0, Nx+1):
        for j in range(0, Ny+1):
            
            vort_plot[i, j] = ((v[i, j]-v[i-1, j])/dx 
                - (u[i, j]-u[i, j-1])/dy)

    return vort_plot


@njit
def calculate_stream_function(
        psi, 
        u, v,
        Nx, Ny,
        dx, dy, 
        omega, 
        tol, iterations
        ):
    """
    This function calculates the streamfunction. 
    It is not necessary to call it at every time step. 
    It can be called only at the end of the simulation. 
    Here, the SOR method is used.
    """
    
    dx2 = dx*dx
    dy2 = dy*dy
    
    error_max = 10.0

    iterations = 0

    # Loop of SOR. 
    while error_max >= tol:

        # Error.
        R_max = 0.0

        for i in range(1, Nx):
            for j in range(1, Ny):

                c = -2.0/(dx2) - 2.0/(dy2)

                Ax = ((psi[i+1, j] - 2.0*psi[i, j] + psi[i-1, j])/(dx2)
                    + (psi[i, j+1] - 2.0*psi[i, j] + psi[i, j-1])/(dy2))

                b = -(v[i, j] - v[i-1, j])/dx + (u[i, j] - u[i, j-1])/dy

                R = (b - Ax)/c

                # psi update.
                psi[i, j] = psi[i, j] + omega*R

                # Error update.
                if abs(R) > R_max:
                    R_max = abs(R)
                    
        iterations += 1

        error_max = R_max

    return psi, iterations    


# Exports the results. 
def export_results(u, v, pressure):
    """
    This function exports the results in csv format. 
    """

    # Velocity u.
    df = pd.DataFrame(u)
    df.to_csv('output_u.csv')
    
    # Velocity v.
    df = pd.DataFrame(v)
    df.to_csv('output_v.csv')

    # Pressure. 
    df = pd.DataFrame(pressure)
    df.to_csv('output_p.csv')
    
    return None
