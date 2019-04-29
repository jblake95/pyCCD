"""
Align and subtract two overlapping CCD images
"""

import argparse as ap
import cv2
import sep

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
	
	
