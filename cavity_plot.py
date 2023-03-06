# In this file we have functions that are used to plot
# the results.

# To run this file independently, you will need to enter 
# the required parameters at the end of the file.

# Attention: if the plots of vorticity and streamfunction are not
# good, you have to change the intervals vector for the contourplot. 


# Libraries.
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from cavity_general_functions import *


# Plot the results. 
def plot_results(u, v, pressure, Re, Lx, Ly, Nx, Ny, dx, dy):
    """
    Plot the results of the simulation. 
    """

    # To prepare for plotting, we must first 
    # interpolate the velocity and pressure at 
    # the mesh points. This is necessary because 
    # our calculations utilize a staggered grid. 
    # However, for the purposes of visualization, 
    # we require the variables to be located at 
    # the same positions.
    u_plot = np.zeros((Nx+1, Ny+1), float)
    v_plot = np.zeros((Nx+1, Ny+1), float)
    pressure_plot = np.zeros((Nx+1, Ny+1), float)

    for i in range(0, Nx+1):
        for j in range(0, Ny+1):
        
            u_plot[i, j] = 0.5*(u[i, j]+u[i, j-1])
            
            v_plot[i, j] = 0.5*(v[i, j]+v[i-1, j])
            
            pressure_plot[i, j] = 0.25*(pressure[i, j] 
                + pressure[i-1, j] + pressure[i, j-1] 
                + pressure[i-1, j-1])

    # The vectors x and y below are used to plot. 
    x = np.linspace(0.0, Lx, Nx+1)    # Vetor x.
    y = np.linspace(0.0, Ly, Ny+1)    # Vetor y.

    # This function plots the streamlines and the pressure. 
    plot_streamlines_and_pressure(
        u_plot, v_plot, 
        pressure_plot, Re, Nx, Ny,
        Lx, Ly, 
        dx, dy, 
        x, y, 
        )
        
    # This function plots the vorticity.     
    plot_vorticity(
        u, v, 
        Re, 
        Nx, Ny,
        Lx, Ly, 
        dx, dy, 
        x, y, 
        )

    # This function plots the stream function.
    plot_stream_function(
        u, v, 
        Re, 
        Nx, Ny,
        Lx, Ly, 
        dx, dy, 
        x, y, 
        )

    # The following function is designed to compare 
    # our results with those obtained by Ghia, Ghia 
    # and Shin (1982) and Marchi, Suero and Araki (2009)
    # for a square cavity of size 1 and Reynolds numbers 
    # (Re) of 1, 10, 100, 400, or 1000. Specifically, we 
    # compare the u velocity at x = 0.5 for each set of results.
    if (Re == 1.0 or Re == 10.0 or Re == 100.0 or 
        Re == 400.0 or Re == 1000.0) and (Lx == 1.0 and
        Ly == 1.0):

        plot_u_versus_y(
            u_plot, v_plot, 
            Re, Nx, Ny,
            Lx, Ly, 
            dx, dy, 
            x, y, 
            )
                
    return None
    

# Plot of streamlines and pressure. 
def plot_streamlines_and_pressure(
        u_plot, v_plot, 
        pressure_plot,
        Re, Nx, Ny,
        Lx, Ly, 
        dx, dy, 
        x, y, 
        ):
    """
    Plot the streamlines generated by python 
    and the pressure contour. 
    """
    
    fig, ax = plt.subplots()

    # Pressure.
    cs = ax.contourf(
        x, y, np.transpose(pressure_plot), 
        cmap='coolwarm'
        )

    fig.colorbar(cs)

    # Streamlines. 
    cs = ax.streamplot(
        x, y, 
        np.transpose(u_plot), np.transpose(v_plot),
        density=2, minlength=0.2, color='k'
        )

    ax.set_xlabel(r'$x$', fontsize=24, usetex=True)
    ax.set_ylabel(r'$y$', fontsize=24, usetex=True)

    # Grid.
    # ax.grid(color='gray')

    plt.xticks(size=16, usetex=True)
    plt.yticks(size=16, usetex=True)

    minor_locator = AutoMinorLocator(5)
    ax.xaxis.set_minor_locator(minor_locator)
    minor_locator = AutoMinorLocator(5)
    ax.yaxis.set_minor_locator(minor_locator)

    # Scaled cavity, for rectangular geometries. 
    plt.axis('scaled')

    ax.set_xlim([0.0, Lx])
    ax.set_ylim([0.0, Ly])

    # Save the figure. 
    plt.savefig(
        'fig1_streamlines_and_pressure',
        format='pdf', 
        dpi=100, 
        bbox_inches='tight'
        )
                
    plt.close()


# Plot of vorticity. 
def plot_vorticity(
        u, v, 
        Re, 
        Nx, Ny,
        Lx, Ly, 
        dx, dy, 
        x, y, 
        ):
    """
    Plot of the vorticity vector. 
    """
    
    # Vorticity array. 
    vort_plot = np.zeros((Nx+1, Ny+1), float)
    
    # Function to calculate the vorticity.
    calculate_vorticity(
        vort_plot,
        u, v,
        Nx, Ny,
        dx, dy
        )

    # Creating the interval vector for plotting. 
    # This vector improves visualization and should 
    # be adjusted according to the result.
    intervals = np.hstack((
        np.linspace(np.min(vort_plot), -0.01, 10),
        np.linspace(0.01, np.max(vort_plot), 10)
        ))

    # For Re = 1000 and Lx = Ly = 1, the intervals below is good. 
    if Re == 1000.0:
        intervals = [
            np.min(vort_plot), -10.0, -5.0,  
            -3.0, -2.08, -1.5, 0.0, 2.0, 5.0, 
            np.max(vort_plot)
            ] 
 
    fig, ax = plt.subplots()

    cs = ax.contour(
        x, y, np.transpose(vort_plot), 
        intervals, 
        colors='black'
        )

    ax.clabel(cs, inline=True)

    # Colors. 
    '''
    cs = ax.contourf(
        x, y, np.transpose(vort_plot), 
        intervals, 
        cmap='coolwarm'
        )
    '''
    
    #fig.colorbar(cs)

    ax.set_xlabel(r'$x$', fontsize=24, usetex=True)
    ax.set_ylabel(r'$y$', fontsize=24, usetex=True)

    # Grid.
    # ax.grid(color='gray')

    plt.xticks(size=16, usetex=True)
    plt.yticks(size=16, usetex=True)

    minor_locator = AutoMinorLocator(5)
    ax.xaxis.set_minor_locator(minor_locator)
    ax.yaxis.set_minor_locator(minor_locator)

    # Scaled cavity, for rectangular geometries.     
    plt.axis('scaled')

    # Save the figure. 
    plt.savefig(
        'fig2_vorticity', 
        format='pdf',
        dpi=100,
        bbox_inches='tight', 
        )
                
    plt.close()                

    return


# Plot of streamfunction.
def plot_stream_function(
        u, v, 
        Re, 
        Nx, Ny,
        Lx, Ly, 
        dx, dy, 
        x, y, 
        ):
    """
    Plot of the streamfunction.
    """

    # Streamfunction array. 
    psi = np.zeros((Nx+1, Ny+1), float)

    # Calculating the streamfunction. 
    iterations_psi = 0
    omega = 1.95
    tol = 1.e-8
    psi, iterations_psi = calculate_stream_function(
        psi,
        u, v, 
        Nx, Ny,
        dx, dy, 
        omega, tol,
        iterations_psi
        )
        
    # Print on screen the maximum value of psi. 
    print(f'\n \n  Max abs value of Psi: {np.amax(abs(psi))} \n \n')        

    fig, ax = plt.subplots()

    # Creating the interval vector for plotting psi. 
    # This vector improves visualization and should 
    # be adjusted according to the result.
    # For Re = 1000 and Lx = Ly = 1, the intervals below is good. 
    intervals = np.hstack((
        np.linspace(np.min(psi), -1.1e-2, 10),
        np.linspace(-4.e-3, -1.e-5, 2),
        np.linspace(1.e-6, 1.e-5, 2),
        np.linspace(1.1e-4, max(1.e-3, np.max(psi)), 6)
        ))
    
    # The code below is only to avoid error when psi is zero
    # or when the intervals vector is not good. 
    if np.max(np.abs(psi)) <= 1.e-6:
        pass
    else:
        try:   
            cs = ax.contour(x, y, np.transpose(psi), intervals, colors='black')
            cs = ax.contourf(x, y, np.transpose(psi), intervals, cmap='rainbow')    
        except:
            cs = ax.contour(x, y, np.transpose(psi), colors='black')
            cs = ax.contourf(x, y, np.transpose(psi), cmap='coolwarm')

    # Line labels.
    #ax.clabel(cs, inline=True)

    # Colorbar. 
    #fig.colorbar(cs)

    ax.set_xlabel(r'$x$', fontsize=24, usetex=True)
    ax.set_ylabel(r'$y$', fontsize=24, usetex=True)

    # Grid. 
    # ax.grid(color='gray')

    plt.xticks(size=16, usetex=True)
    plt.yticks(size=16, usetex=True)

    minor_locator = AutoMinorLocator(5)
    ax.xaxis.set_minor_locator(minor_locator)
    ax.yaxis.set_minor_locator(minor_locator)

    # Scaled cavity, for rectangular geometries. 
    plt.axis('scaled')

    # Save the figure. 
    plt.savefig(
        'fig3_streamfunction',
        format='pdf',
        dpi=100, 
        bbox_inches='tight', 
        )
    
    plt.close()
    
    return


# Plot u versus y and compares with other works.
def plot_u_versus_y(
        u_plot, v_plot, 
        Re, Nx, Ny,
        Lx, Ly, 
        dx, dy, 
        x, y, 
        ):
    """
    Plots u as a function of y for x = 0.5 and
    compares with the results of Ghia, Ghia and 
    Shin (1982) and Marchi, Suero and Araki (2009). 
    """

    fig, ax = plt.subplots()

    if Nx % 2 == 0:
        ax.plot(
            u_plot[int(Nx/2), :], y, '-k',
            label='This work')
    else:
        aux = int(Nx/2)
        ax.plot(
            0.5*(u_plot[aux, :] + u_plot[aux + 1, :]), y, '-k',
            label='This work')

    # Matrix with the results from Ghia. 
    # The first column shows the values of y, 
    # the second shows the values of u for Re = 100, 
    # the third for Re = 400, and the fourth for Re = 1000.
    results_ghia = np.array([
        1.0000,	 1.0000, 1.00000, 1.0000,	     
        0.9766,	0.84123, 0.75837, 0.65928,
        0.9688,	0.78871, 0.68439, 0.57492,
        0.9609,	0.73722, 0.61756, 0.51117,
        0.9531,	0.68717, 0.55892, 0.46604,
        0.8516,	0.23151, 0.29093, 0.33304,
        0.7344,	0.00332, 0.16256, 0.18719,
        0.6172,	-0.1364, 0.02135, 0.05702,
        0.5000,	-0.2058, -0.1148, -0.0608,	
        0.4531,	-0.2109, -0.1712, -0.1065,
        0.2813,	-0.1566, -0.3273, -0.2781,
        0.1719,	-0.1015, -0.2430, -0.3829,
        0.1016,	-0.0643, -0.1461, -0.2973,
        0.0703,	-0.0478, -0.1034, -0.2222,
        0.0625,	-0.0419, -0.0927, -0.2020,
        0.0547,	-0.0372, -0.0819, -0.1811,
        0.0000,	 0.0000,  0.0000,  0.0000
        ])
    
    y_ghia = results_ghia[0::4]    
    u_ghia = ((Re==100)*results_ghia[1::4]
        + (Re==400)*results_ghia[2::4]        
        + (Re==1000)*results_ghia[3::4])

    if Re == 100 or Re == 400 or Re == 1000:  
        ax.plot(
            u_ghia, y_ghia, 'sk', 
            fillstyle='none', ms=7.0,
            mew=1.5,
            label='Ghia, Ghia and Shin (1982)'
            )         

    # Now we get the results of Marchi (2009). 
    y_marchi = np.linspace(0.0625, 0.9375, 15)
    
    if Re == 0.01:
        u_marchi = np.array([
            -3.85275436e-2, -6.9584425e-2, -9.6906717e-2,
            -1.22595555e-1, -1.47461728e-1, -1.71067124e-1, 
            -1.91535923e-1, -2.05191715e-1, -2.06089397e-1, 
            -1.85581148e-1, -1.322092275e-1, -3.2443684e-2, 
             1.27054983e-1,  3.55228331e-1, 6.51176326e-1  
            ])        
    if Re == 10:
        u_marchi = np.array([
            -3.85425800e-2, -6.96238561e-2, -9.6983962e-2, 
            -1.22721979e-1, -1.47636199e-1, -1.71260757e-1, 
            -1.91677043e-1, -2.05164738e-1, -2.05770198e-1, 
            -1.84928116e-1, -1.313892353e-1, -3.1879308e-2, 
            1.26912095e-1, 3.54430364e-1, 6.50529292e-1
            ])

    if Re == 100:
        u_marchi = np.array([
            -4.1974991e-2, -7.7125399e-2, -1.09816214e-1, 
            -1.41930064e-1, -1.72712391e-1, -1.98470859e-1, 
            -2.12962392e-1, -2.091491418e-1, -1.82080595e-1, 
            -1.31256301e-1, -6.0245594e-2, 2.7874448e-2,
            1.40425325e-1, 3.1055709e-1, 5.97466694e-1
            ])
    
    if Re == 400:
        u_marchi = np.array([
            -9.259926e-2, -1.78748051e-1, -2.6391720e-1, 
            -3.2122908e-1, -3.2025109e-1, -2.6630635e-1,
            -1.9073056e-1, -1.15053628e-1, -4.2568947e-2, 
            3.024302e-2, 1.0545601e-1, 1.8130685e-1, 
            2.5220384e-1, 3.1682969e-1, 4.69580199e-1
            ])        
    
    if Re == 1000:
        u_marchi = np.array([
            -2.02330048e-1, -3.478451e-1, -3.844094e-1,
            -3.189461e-1, -2.456937e-1, -1.837321e-1,
            -1.2341046e-1, -6.205613e-2, 5.6180e-4, 
            6.5248742e-2, 1.3357257e-1, 2.0791461e-1, 
            2.884424e-1, 3.625454e-1, 4.229321e-1
            ])
        
    if (Re == 0.01 or Re == 10.0 or Re == 100.0 or 
        Re == 400.0 or Re == 1000.0): 
        ax.plot(
            u_marchi, y_marchi, 'ok', 
            fillstyle='none', ms=7.0,
            mew=1.5,
            label='Marchi, Suero and Araki (2009)'
            )
    
    ax.legend()                
        
    ax.set_xlabel(r'$u$', fontsize=24, usetex=True)
    ax.set_ylabel(r'$y$', fontsize=24, usetex=True)

    # Grid.
    ax.grid(color='gray')
    ax.grid(which='minor', color='gray', alpha=0.1)

    plt.xticks(size=16, usetex=True)
    plt.yticks(size=16, usetex=True)

    minor_locator = AutoMinorLocator(5)
    ax.xaxis.set_minor_locator(minor_locator)
    minor_locator = AutoMinorLocator(5)
    ax.yaxis.set_minor_locator(minor_locator)

    ax.set_xlim([np.min(u_plot) - 0.05, np.max(u_plot) + 0.05])
    ax.set_ylim([-.05, Ly + 0.05])

    # Save the figure. 
    plt.savefig(
        'fig4_u_versus_y',
        format='pdf', 
        dpi=100, 
        bbox_inches='tight'
        )
                
    plt.close()

    return
    
    
# This is the main function. 
# We can run this file independently.  
def plot_main():

    # Reynolds number.
    Re = 1000.0    

    # Dimensions of the cavity.
    Lx = 1.0    # x direction.
    Ly = 1.0    # y direction.

    # Number of spaces in the x and y directions. 
    Nx = 100
    Ny = 100
    
    u = np.array(
        pd.DataFrame.to_numpy(pd.read_csv('output_u.csv'))[:, 1:]
        )
        
    v = np.array(
        pd.DataFrame.to_numpy(pd.read_csv('output_v.csv'))[:, 1:]
        )
        
    pressure = np.array(
        pd.DataFrame.to_numpy(pd.read_csv('output_p.csv'))[:, 1:]
        )
        
    plot_results(u, v, pressure, Re, Lx, Ly, Nx, Ny, Lx/Nx, Ly/Ny)
    

# The plots can be automatically generated from 
# the file "cavity_main.py", or they can be 
# generated from this file. 
# That is why we have the following "if" statement.
if __name__ == '__main__':
    plot_main()
