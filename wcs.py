"""
Functions for obtaining World Coordinate System (WCS) solutions of 
CCD frames using the Astrometry.net package
"""

import os
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits

try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

def solveField(filename, file_prefix, bintable=True,
               input_dir='', output_dir=None, 
               ra=None, dec=None, radius=2.,
               scale_low=0.3, scale_high=0.4, nx=8176, ny=6132):
    """
    Solve field using Astrometry.net for images or source tables
    
    Parameters
    ----------
    filepath : str
        Path to FITS file containing either image data or a table with
        pixel centroids (x,y) and flux for detected sources
    file_prefix : str
        Name of output World Coordinates System file containing the
        astrometric solution and Correspondence file
    bintable : bool, optional
        Toggle to change between FITS image and FITS bintable format
        Default = True (bintable)
    out_dir : str, optional
        Path to directory in which to place output files
        Default = None, will use input directory
    ra : float, optional
        Right ascension of image center in degrees or hh:mm:ss
        Default = None
    dec : float, optional
        Declination of image center in degrees or [+-]dd:mm:ss
        Default = None
    radius : float, optional
        Only search in indices within 'radius' of the field
        center given by ('ra', 'dec')
        Default = 2.
    scale_low : float, optional
        Lower bound of image scale estimate in arcsecperpixel
        Default = 0.3 [INT WFC]
    scale_high : float, optional
        Upper bound of image scale estimate in arcsecperpixel
        Default = 0.4 [INT WFC]
    nx : int, optional
        Width of image in pixels
        Default = 8176 [INT frame]
    ny : int, optional
        Height of image in pixels
        Default = 6132 [INT frame]
    
    Returns
    -------
    None
    """
    
    astrom_loc='/usr/bin/'              # location of astrometry.Net
    cmd = (
        '{}solve-field '                # call solve-field function 
        '{} '                           # file containing field to solve
        '--no-verify '                  # ignore existing wcs info
        '--no-fits2fits '               # don't sanitize fits files
        '--no-plots '                   # don't create plots of result
        '--crpix-center '               # set wcs reference to center
        '--new-fits none '              # no new fits file
        '--wcs {} '                     # name of wcs output file
        '--solved none '                # no solved file
        '--match none '                 # no match file
        '--rdls none '                  # no rdls file
        '--corr {} '                    # name of corr output file
        '--axy none '                   # no axy file
        '--index-xyls none '            # no index-xyls file
        '--overwrite '                  # overwrite existing outputs
        '--scale-low {} '               # lower bound of scale estimate
        '--scale-high {} '              # upper bound of scale estimate
        '--scale-units arcsecperpix '   # units of scale estimates
        ).format(astrom_loc,
                 filepath,
                 file_prefix+'.wcs',
                 file_prefix+'.corr',
                 scale_low,
                 scale_high)
    
    if out_dir is not None:
        cmd = (
            '{} '                       # previous command
            '--dir {} '                 # name of output directory
            ).format(cmd, 
                     out_dir)
    
    if bintable:
        cmd = (
            '{} '                       # previous command
            '--x-column x_det '         # name of column containing x
            '--y-column y_det '         # name of column containing y
            '--sort-column flux '       # column to sort by
            '--width {} --height {} '   # width/height of field (pixels)
            ).format(cmd, 
                     nx, 
                     ny)
    
    if ra and dec is not None:
        cmd = (
            '{} '                       # previous command
            '--ra {} '                  # right ascension of center
            '--dec {} '                 # declination of center
            '--radius {}'               # radius within which to look
            ).format(cmd, 
                     str(ra), 
                     str(dec), 
                     str(radius))
    print(cmd)
    os.system(cmd)
    
    return None

def convertToDetector(x, y, hdu_hdr):
    """
    Convert a list of xy pixel coordinates to detector coordinates for
    a CCD mosaic (e.g. INT WFC)
    
    Parameters
    ----------
    x, y : array-like
        Lists of xy positions in the CCD image
    hdu_hdr : str
        FITS header for HDU, containing transformation coefficients
    
    Returns
    -------
    x_det, y_det : array-like
        Detector coords corresponding to input xy coords
    """
    w = WCS(hdu_hdr)
    xy_coords = np.column_stack([x, y])
    
    # FITS convention, so use Fortran-like 1-based origin
    detector = w.all_pix2world(xy_coords, 1)
    x_det, y_det = detector[:, 0], detector[:, 1]
    
    return x_det, y_det

def convertToWCS(x, y, wcs_hdr):
    """
    Convert a list of xy pixel coordinates to (ra,dec) coordinates
    
    Parameters
    ----------
    x, y : array-like
        Lists of xy positions in the WCS-solved CCD image
    wcs_hdr : str
        FITS header containing WCS solution
    
    Returns
    -------
    ra, dec : array-like
        WCS coords corresponding to input xy coords
    """
    w = WCS(wcs_hdr)
    xy_coords = np.column_stack([x, y])
    
    # FITS convention, so use Fortran-like 1-based origin
    world = w.all_pix2world(xy_coords, 1)
    ra, dec = world[:, 0], world[:, 1]
    
    return ra, dec

def convertToPixels(ra, dec, wcs_hdr):
    """
    Convert a list of (ra,dec) coordinates to xy pixel coordinates
    
    Parameters
    ----------
    ra, dec : array-like
        Lists of world coords to be converted to pixel coords
    wcs_hdr : str
        FITS header containing WCS solution
    
    Returns
    -------
    x, y : array-like
        Pixel coords corresponding to input world coords
    """
    w = WCS(wcs_hdr)
    world_coords = np.column_stack([ra, dec])
    
    # FITS convention, so use Fortran-like 1-based origin
    pixel = w.all_world2pix(world_coords, 1)
    x, y = pixel[:, 0], pixel[:, 1]
    
    return x, y 
