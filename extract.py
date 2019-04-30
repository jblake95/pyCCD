"""
Functions for performing source extraction on CCD images
-Background subtraction
-Threshold extraction of sources
"""

from utils import pruneNansFromTable
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

def sourceExtract(data, thresh=3, bkg=False, bkg_rms=None, 
                  err=None, mask=None, min_area=5, 
                  deblend_cont=0.05, segment=False, extras=False):
    """
    Extract all sources above a certain threshold in the given image
    
    Parameters
    ----------
    data : array-like
        CCD image frame from which to extract sources
    thresh : float, optional
        Number of sigma a detection must be above the background to 
        be flagged as a source - if err not given, bkg_rms is needed
        Default = 3
    bkg : bool, optional
        Toggle to model spatially varying background and subtract from
        data - by default assumes this has been done separately
        Default = False
    bkg_rms : float, optional
        Estimation of the global background noise - used to determine
        threshold if err is None - can calculate global background
        rms when subtracting background model
        Default = None
    err : array-like, optional
        Error array for the CCD frame - supersedes bkg_rms when
        determining the threshold
        Default = None
    min_area : int, optional
        Minimum number of pixels to be flagged as a source
        Default = 5
    deblend_cont : float, optional
        Minimum contrast ratio used by SEP for deblending
        Default = 0.05
    segment : bool, optional
        Toggle to generate a segmentation map for the given image
        Default = False
    extras : bool, optional
        Toggle to calculate ellipticity, FWHM, Kron radius and 
        flux radius
        Default = False
    
    Returns
    -------
    sources : astropy Table object
        Table containing quantities determined by sep for each source
        detected in the given image
    segmentation_map : array-like, optional
        Array of integers with same shape as data - pixels not 
        belonging to any object have value 0, whilst all pixels 
        belonging to ith object have value (e.g. sources[i]) have
        value i+1 - only returned if seg_map is True
    
    Raises
    ------
    None
    """
    
    # subtract spatially varying background model if requested
    if bkg:
        data, bkg_rms = subtractBackground(data)
    
    # determine threshold for extraction
    if err is None:
        thresh *= bkg_rms
    
    # extract sources
    if not segment:
        sources = sep.extract(data, 
                              thresh, 
                              err=err,
                              mask=mask,
                              deblend_cont=deblend_cont)
    else:
        sources, seg_map = sep.extract(data, 
                                       thresh,
                                       err=err,
                                       mask=mask,
                                       deblend_cont=deblend_cont,
                                       segmentation_map=seg_map)
    
    sources = Table(sources)
    
    # remove nans from table
    sources = pruneNansFromTable(sources)
    
    if extras:
		# calculate ellipticity parameter
		sources['ellipticity'] = 1.0 - (sources['b'] / sources['a'])
		
		# calculate full width half maxima
		sources['fwhm'] = calculateFWHM(sources['a'], sources['b'])
		
		# compute kron radii
		try:
			sources['kronr'], krflag = sep.kron_radius(data, 
													   sources['x'], 
													   sources['y'], 
													   sources['a'], 
													   sources['b'], 
													   sources['theta'], 
													   6.0)
			sources['flag'] |= krflag
		except Exception as e:
			print(e)
			pass
		
		# compute flux radii
		try:
			sources['fluxr'], frflag = sep.flux_radius(data,
													   sources['x'],
													   sources['y'],
													   6.0*sources['a'],
													   0.5,
													   subpix=5)
			sources['flag'] |= frflag
		except Exception as e:
			print(e)
			pass
    
    if segment:
        return sources, seg_map
    else:
        return sources

def calculateFWHM(a, b):
    """
    Calculate the FWHM of sources detected by SEP
    
    Parameters
    ----------
    a, b : array-like
        Elliptical parameters supplied by SEP
    
    Returns
    -------
    fwhm : array-like
        Full width half maxima of sources
    """
    return 2 * np.sqrt(np.log(2) * (a**2 + b**2))
