"""
Functions for performing source extraction on CCD images
-Background subtraction
-Threshold extraction of sources
"""

import sep

def subtractBackground(data, mask=None, box_width=32, box_height=32, 
                       filter_width=3, filter_height=3):
    """
    Determine the spatially varying sky background using SEP
    and subtract from the image
    
    Parameters
    ----------
    data : array-like
        CCD data from which to subtract the background
    mask : array-like
        Bad pixel mask for the CCD frame
    box_width : int, optional
        Width of background boxes in pixels
        Default = 32
    box_height : int, optional
        Height of background boxes in pixels
        Default = 32
    filter_width : int, optional
        Width of filter in boxes
        Default = 3
    filter_height : int, optional
        Height of filter in boxes
        Default = 3
    
    Returns
    -------
    data_sub : array-like
        Data array with background signal subtracted
    bkg_rms : float
        Global rms of the spatially varying background, for use as a 
        backup threshold in the extraction procedure
    """
    
    # FITS files can be backwards byte order - SEP needs this fixed
    try:
        bkg = sep.Background(data, 
                             mask=mask,
                             bw=box_width, 
                             bh=box_height,
                             fw=filter_width, 
                             fh=filter_height)
    except:
        data = data.byteswap(True).newbyteorder()
        bkg = sep.Background(data, 
                             mask=mask,
                             bw=box_width, 
                             bh=box_height,
                             fw=filter_width, 
                             fh=filter_height)
    
    data_sub = data - bkg
    
    # calculate background rms as a backup for extraction threshold
    bkg_rms = bkg.globalrms
    
    return data_sub, bkg_rms
