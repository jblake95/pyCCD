"""
Align and subtract two overlapping INT WFC FITS CCD images
***TODO: generalise to include sidereally tracked images***
"""

from wcs import (
    convertToDetector,
    convertToWCS,
    convertToPixels,
    )
from diagnostics import plotXY
import argparse as ap
import numpy as np
import cv2
import random
import warnings
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from astropy.utils.exceptions import AstropyWarning
from astropy.coordinates import SkyCoord, match_coordinates_sky

try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

# disable astropy warnings - INT WCS is deprecated
warnings.simplefilter('ignore', category=AstropyWarning)

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
	
	####################################################################
	######### load the images, WCS headers and bad pixel mask ##########
	####################################################################
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
	
	####################################################################
	################# feature extraction and matching ##################
	####################################################################
	# inject random sampling points into image 1
	x_1 = []
	y_1 = []
	while len(x_1) < 1000:
		x_1.append(random.uniform(0, hdr_1['NAXIS1']))
		y_1.append(random.uniform(0, hdr_1['NAXIS2']))
	
	if args.diagnostics:
		plotXY(img_1, x_1, y_1)
	
	# convert xy to radec
	xdet_1, ydet_1 = convertToDetector(x_1, 
	                                   y_1, 
	                                   hdr_1)
	ra_1, dec_1 = convertToWCS(xdet_1, 
	                           ydet_1, 
	                           wcs_1)
	
	# match with image 2
	xdet_2, ydet_2 = convertToPixels(ra_1,
	                                 dec_1,
	                                 wcs_2)
	x_2, y_2 = convertToPixels(xdet_2,
	                           ydet_2,
	                           hdr_2)
	matches = Table([x_1, y_1, 
	                 xdet_1, ydet_1, 
	                 ra_1, dec_1, 
	                 xdet_2, ydet_2, 
	                 x_2, y_2],
	                 names=['x_1', 'y_1',
	                        'xdet_1', 'ydet_1',
	                        'ra_1', 'dec_1',
	                        'xdet_2', 'ydet_2',
	                        'x_2', 'y_2'])
	
	# filter out points that don't overlap between the images
	remove_idx = []
	for m, match in enumerate(matches):
		if (match['x_2'] < 0 or
		    match['x_2'] > hdr_2['NAXIS1'] or
		    match['y_2'] < 0 or
		    match['y_2'] > hdr_2['NAXIS2']):
			   remove_idx.append(m)
	matches.remove_rows(remove_idx)
	
	if args.diagnostics:
		plotXY(img_2, matches['x_2'], matches['y_2'])
	
	####################################################################
	#################### alignment and subtraction #####################
	####################################################################
	
