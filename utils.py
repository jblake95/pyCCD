"""
General functions for pyCCD scripts
"""

import numpy as np

def pruneNansFromTable(table):
    """
    Remove rows from a table if they contain nans
    
    Parameters
    ----------
    table : astropy Table object
        Table to be pruned of nans
    
    Returns
    -------
    table : astropy Table object
        Table pruned of nans
    """
    nan_in_row = np.zeros(len(table), dtype=np.bool)
    for col in table.colnames:
        nan_in_row |= np.isnan(table[col])
        
    return table[~nan_in_row]

def getTrailLength(exptime, platescale, rate=15.034):
    """
    Calculate a rough estimate of expected trail length based on 
    an informed rate and exposure time
    """
    return (exptime * rate) / platescale
