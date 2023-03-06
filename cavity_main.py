# Driven Cavity Problem solved by Finite Differences and
# coupled with an Explicit First Order Projection Method.
# Implemented by Adriano Possebon Rosa.
# Undergraduate Course: Numerical Methods in Thermofluids. 
# University of Brasília.
# 2023/03/02.


# We start the code by importing the libraries.
# Numpy is used to implement the arrays.
import numpy as np
# Matplotlib is used to make the plots, later.
import matplotlib.pyplot as plt
# Timeit library is used only to calculate the machine time
# of the simulation.
import timeit

# We will also use functions that are in other files. 
# We put functions in separated files just to keep it 
# more organized. 
from cavity_velocity import *
from cavity_pressure import *
from cavity_general_functions import *
from cavity_plot import *


# start is the machine time at the beginning of the simulation.
start = timeit.time.time()


# We need to specify some parameters of the simulation. 
# These are the inputs. 

# Reynolds number.
Re = 1000.0    

# Dimensions of the cavity.
Lx = 1.0    # x direction.
Ly = 1.0    # y direction.

# Number of spaces in the x and y directions. 
Nx = 100
Ny = 100

# Time-step.
dt = 0.005

# Final time.
time_final = 40.0

# Value of omega used in SOR for pressure.
omega_p = 1.95

# Value of the tolerance when solving the pressure with SOR.
tol_p = 1.e-7

# Velocity on the top wall.
u_top = np.ones((Nx+1), float)    # Velocity is uniform.


# Ok. We have given all of the inputs. 
# Now we will define other parameters, based on the inputs.
# And then we will create the vectors and arrays. 

# dx and dy.
dx = Lx/Nx
dy = Ly/Ny

# Before we start the simulation, let's check if the 
# parameters are acceptable, based on the stability conditions. 
check_parameters(Re, dx, dy, dt)
  
# We are going to create the arrays.
# The most important arrays are the arrays
# of u, v, u_star, v_star and pressure. 

# Velocity in the x direction.
u = np.zeros((Nx+1,Ny+2), float)
u_star = np.copy(u)

# Boundary condition for u in the top.
u[:,Ny] = 2.0*u_top[:]

# Velocity in the y direction. 
v = np.zeros((Nx+2,Ny+1), float)
v_star = np.copy(v)
v_new = np.copy(v)

# Pressure array.
pressure = np.zeros((Nx+2,Ny+2),float)

# It is important to know the number of iterations necessary
# in the pressure calculation. 
pressure_iterations = 0

# Divergence of Velocity at each point.
# We will calculate the divergence at every point.
# In this way we will know if the divergence of the
# velocity is indeed small. 
div = np.zeros((Nx,Ny),float)


# Here is the main part of the code: the time loop.
# We make the time evolution using a while loop. 
# We need to create the time variable. 
time = 0.0

# It is also important to know the iteration number, 
# in time. For this purpose, we create the integer k. 
k = 0

# Time loop starts here. 
while time < time_final - 1.e-5*dt:

    # u_star.
    u_star = calculate_u_star(u_star, u, v, u_top, Re, Nx, Ny, dx, dy, dt)
    # v_star.
    v_star = calculate_v_star(v_star, u, v, Re, Nx, Ny, dx, dy, dt)

    # pressure.
    # The pressure calculation is the most expensive part 
    # of the simulation. 
    pressure, pressure_iterations = calculate_pressure(
        pressure, 
        pressure_iterations, 
        u_star, v_star, 
        Nx, Ny, 
        dx, dy, dt, 
        omega_p, tol_p)

    # New velocity u.
    u = calculate_u_new(u, u_star, pressure, Nx, Ny, dx, dt)
    # New velocity v. 
    v = calculate_v_new(v, v_star, pressure, Nx, Ny, dy, dt)

    # Calculating the divergence of the velocity.
    # This is optional. 
    div = check_for_incompressibility(div, u, v, Nx, Ny, dx, dy)

    # Time update. 
    time = time + dt  
    # Time iterations update. 
    k = k + 1 
    
    # Now, we will print some information on the screen.
    # This is important so that we can monitor the progress 
    # of the simulation.
    print_on_screen(time, time_final, pressure_iterations, div)


# The time loop ends here. 

# Here we take the machine time and print the execution total time 
# in the screen. 
end = timeit.time.time()
print(f'\n\n Execution time: {end - start:.4} s. \n\n')

# Export results. 
export_results(u, v, pressure)

# Plot the results. 
plot_results(u, v, pressure, Re, Lx, Ly, Nx, Ny, dx, dy)

