"""
Align and subtract two overlapping INT WFC FITS CCD images
***TODO: generalise to include sidereally tracked images***
"""

from extract import (
    subtractBackground,
    sourceExtract,
    )
from wcs import (
    convertToDetector,
    convertToWCS,
    )
from utils import getTrailLength
from diagnostics import plotSources
import argparse as ap
import numpy as np
import cv2
import sep
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky

try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

def argParse():
    """
    Argument parser settings
    
    Parameters
    ----------
    None
    
    Returns
    -------
    args : array-like
        Array of command line arguments
    """
    parser = ap.ArgumentParser()
    
    parser.add_argument('img_1',
                        help='path to CCD frame 1',
                        type=str)
    
    parser.add_argument('img_2',
                        help='path to CCD frame 2',
                        type=str)
    
    parser.add_argument('wcs_1',
                        help='path to WCS solution for CCD frame 1',
                        type=str)
    
    parser.add_argument('wcs_2',
                        help='path to WCS solution for CCD frame 2',
                        type=str)
    
    parser.add_argument('bp_mask',
                        help='path to bad pixel mask for the frames',
                        type=str)
    
    parser.add_argument('--diagnostics',
                        help='include sanity checks?',
                        action='store_true')
    
    return parser.parse_args()

if __name__ == "__main__":
	
	args = argParse()
	
	## load the images, WCS headers and bad pixel mask
	print('Loading image 1...')
	try:
		with fits.open(args.img_1) as f1:
			success = np.arange(0, len(f1))
			while True:
				check = input('Please specify relevant HDU: ')
				if int(check) in success:
					print('Proceeding with HDU {}...'.format(check))
					break
				else:
					print('Invalid selection...\n'
						  'Options: {}--{}'.format(str(min(success)),
												   str(max(success))))
			
			img_1 = f1[int(check)].data.astype(np.float64)
			hdr_1 = f1[int(check)].header
			primhdr_1 = f1[0].header
		
	except FileNotFoundError:
		print('Image 1 not found...')
		quit()
	
	print('Loading image 2...')
	try:
		with fits.open(args.img_2) as f2:
			success = np.arange(0, len(f2))
			while True:
				check = input('Please specify relevant HDU: ')
				if int(check) in success:
					print('Proceeding with HDU {}...'.format(check))
					break
				else:
					print('Invalid selection...\n'
						  'Options: {}--{}'.format(str(min(success)),
												   str(max(success))))
		
			img_2 = f2[int(check)].data.astype(np.float64)
			hdr_2 = f2[int(check)].header
			primhdr_2 = f2[0].header
		
	except FileNotFoundError:
		print('Image 2 not found...')
		quit()
	
	print('Loading wcs information...')
	try:
		with fits.open(args.wcs_1) as w1:
			wcs_1 = w1[0].header
	except FileNotFoundError:
		print('No WCS information for image 1...')
		quit()
	try:
		with fits.open(args.wcs_2) as w2:
			wcs_2 = w2[0].header
	except FileNotFoundError:
		print('No WCS information for image 2...')
		quit()
	
	print('Loading bad pixel mask...')
	try:
		with fits.open(args.bp_mask) as bpm:
			success = np.arange(0, len(bpm))
			while True:
				check = input('Please specify relevant HDU: ')
				if int(check) in success:
					print('Proceeding with HDU {}...'.format(check))
					break
				else:
					print('Invalid selection...\n'
						  'Options: {}--{}'.format(str(min(success)),
												   str(max(success))))
		
			mask = bpm[int(check)].data.astype(np.bool)
		
	except FileNotFoundError:
		print('Bad pixel mask not found...')
		quit()
	
	## feature extraction
	img_1, bkg_rms = subtractBackground(img_1, 
	                                    mask=mask)
	
	stars_1 = sourceExtract(img_1,
	                        thresh=5,
	                        mask=mask,
	                        bkg_rms=bkg_rms,
	                        deblend_cont=1.)
	
	# change header keywords as appropriate
	l_trail = getTrailLength(primhdr_1['EXPTIME'],
	                         primhdr_1['SECPPIX']*primhdr_1['CCDXBIN'])
	dl_max = 0.2*l_trail
	
	# remove bad sources, keep stars
	remove_idx = []
	for s, star in enumerate(stars_1):
		l_sep = np.sqrt((star['xmax'] - star['xmin'])**2 + 
		                (star['ymax'] - star['ymin'])**2)
		if abs(l_sep - l_trail) > dl_max:
			remove_idx.append(s)
	
	stars_1.remove_rows(remove_idx)
	
	if args.diagnostics:
		plotSources(img_1, stars_1, circle=True)
	
	img_2, bkg_rms = subtractBackground(img_2, 
	                                    mask=mask)
	
	stars_2 = sourceExtract(img_2,
	                        thresh=2,
	                        mask=mask,
	                        bkg_rms=bkg_rms,
	                        deblend_cont=1.)
	
	# change header keywords as appropriate
	l_trail = getTrailLength(primhdr_2['EXPTIME'],
	                         primhdr_2['SECPPIX']*primhdr_2['CCDXBIN'])
	dl_max = 0.2*l_trail
	
	# remove bad sources, keep stars
	remove_idx = []
	for s, star in enumerate(stars_2):
		l_sep = np.sqrt((star['xmax'] - star['xmin'])**2 + 
		                (star['ymax'] - star['ymin'])**2)
		if abs(l_sep - l_trail) > dl_max:
			remove_idx.append(s)
	
	stars_2.remove_rows(remove_idx)
	
	if args.diagnostics:
		plotSources(img_2, stars_2, circle=True)
	
	## match features between images
	xdet_1, ydet_1 = convertToDetector(stars_1['x'],
									   stars_1['y'],
									   hdr_1)
	xdet_2, ydet_2 = convertToDetector(stars_2['x'],
									   stars_2['y'],
									   hdr_2)
	ra_1, dec_1 = convertToWCS(xdet_1,
							   ydet_1,
							   wcs_1)
	ra_2, dec_2 = convertToWCS(xdet_2,
							   ydet_2,
							   wcs_2)
	stars_1['ra'] = ra_1
	stars_1['dec'] = dec_1
	stars_2['ra'] = ra_2
	stars_2['dec'] = dec_2
	
	coords_1 = SkyCoord(ra=stars_1['ra'], 
					    dec=stars_1['dec'],
					    unit=(u.deg, u.deg))
	coords_2 = SkyCoord(ra=stars_2['ra'], 
					    dec=stars_2['dec'],
					    unit=(u.deg, u.deg))
	
	match_idx, match_sep, _ = match_coordinates_sky(coords_1, coords_2)
	
	print(sorted(match_sep))
