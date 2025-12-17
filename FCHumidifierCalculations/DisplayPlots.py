import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import numpy as np

def Display3DPlot(x, xName, y, yName, z,zName):
    # Let's assume X, Y, and Z are already your lists of data
    X = np.asarray(x)
    Y = np.asarray(y)
    Z = np.asarray(z)

    # Create a meshgrid for X and Y
    X, Y = np.meshgrid(X, Y)

    # If your Z is 1D (100^2 values), reshape it into a 2D array (100x100)
    #Z = Z.reshape((100, 100))

    # Create a figure and 3D axis
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot the surface, with custom colormap for hillshading effect
    surf = ax.plot_surface(X, Y, Z, cmap='RdYlBu', edgecolor='none')

    # Customize the color bar
    fig.colorbar(surf, shrink=0.5, aspect=5)

    # Optionally, adjust the view angle for better visualization
    ax.view_init(azim=45, elev=30)

    ax.set_xlabel(xName)  # Label for X axis
    ax.set_ylabel(yName)  # Label for Y axis
    ax.set_zlabel(zName)  # Label for Z axis

    # Show the plot
    plt.show()


def Display2DPlot(x, xName, y, yName):
    fig, ax1 = plt.subplots()
    ax1.plot(x, y)

    ax1.set(xlabel=xName, ylabel=yName,
   )
    ax1.grid()
    plt.show()