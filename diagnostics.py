"""
Diagnostic functions for pyCCD scripts
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse

def plotSources(data, sources, circle=False):
    """
    Plot sources detected by SEP on top of an image
    
    Parameters
    ----------
    data : array-like
        Image data for the CCD frame
    sources : astropy Table object
        Source catalog outputted by SEP for the frame
    circle : bool, optional
        Toggle to switch to circular indicator
        Default = False [elliptical apertures used]
    
    Returns
    -------
    None
    """
    fig, ax = plt.subplots()
    m, s = np.mean(data), np.std(data)
    im = ax.imshow(data, interpolation='nearest', cmap='gray',
                   vmin=m-s, vmax=m+s, origin='lower')

    for i in range(len(sources)):
        if circle:
            c = Circle(xy=(sources['x'][i], sources['y'][i]),
                       radius=3)
            c.set_facecolor('none')
            c.set_edgecolor('red')
            ax.add_artist(c)
        else:
            e = Ellipse(xy=(sources['x'][i], sources['y'][i]),
                        width=6*sources['a'][i],
                        height=6*sources['b'][i],
                        angle=sources['theta'][i] * 180. / np.pi)
            e.set_facecolor('none')
            e.set_edgecolor('red')
            ax.add_artist(e)
        
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()
    plt.close(fig)

def plotXY(data, x, y):
    """
    Plot xy markers on top of an image
    
    Parameters
    ----------
    data : array-like
        Image data for the CCD
    x, y : array-like
        xy coords for the markers to be placed
    
    Returns
    -------
    None
    """
    fig, ax = plt.subplots()
    m, s = np.mean(data), np.std(data)
    im = ax.imshow(data, interpolation='nearest', cmap='gray',
                   vmin=m-s, vmax=m+s, origin='lower')
    
    for (i, j) in zip(x, y):
        c = Circle(xy=(i, j), radius=3)
        c.set_facecolor('none')
        c.set_edgecolor('red')
        ax.add_artist(c)
    
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()
    plt.close(fig)
