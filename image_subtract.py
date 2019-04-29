"""
Align and subtract two overlapping FITS CCD images
"""

import argparse as ap
import numpy as np
import cv2
import sep
from astropy.io import fits

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
    
    return parser.parse_args()

if __name__ == "__main__":
	
	args = argParse()
	
	# load the images as numpy arrays
	try:
		with fits.open(args.img_1) as f1:
			if len(f1) > 1:
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
	except FileNotFoundError:
		print('Image 1 not found...')
		quit()
	
	try:
		with fits.open(args.img_2) as f2:
			if len(f2) > 1:
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
	except FileNotFoundError:
		print('Image 2 not found...')
		quit()
