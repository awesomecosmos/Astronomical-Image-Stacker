# -*- coding: utf-8 -*-

# This file contains all necessary packages and functions required for the 
# Sidereal Image Stacking of the images. This file is required to run the file 
# sidereal_img_stacker.py.
# Import it as: 
# from sis_funcs import *

###############################################################################
#-------------------SECTION ONE: IMPORTING PACKAGES---------------------------#
###############################################################################

import astroalign as aa
import astropy.units as u
from astropy.io import fits
from astropy.nddata import CCDData

# path-type packages
import os
import glob
from pathlib import Path

###############################################################################
#--------------------SECTION ONE: HELPER FUNCTIONS----------------------------#
###############################################################################

def registered_img_writer(img,registered_img,SIS_stack_path):
    """
    This function writes the original FITS header to a FITS file of the 
    registered image.
    
    Parameters
    ----------
    img : str
        Filename of image from which we want to extract header.
    
    registered_img : array
        Registered image as an array.
    
    SIS_stack_path : WindowsPath object
        Path to directory where sidereally-stacked image is to be saved.
    
    Returns
    -------
    Nothing.
    """
    img_hdul = fits.open(img)
    img_hdr1 = img_hdul[0].header
                
    run_filename = img_hdr1['RUN'].strip(' ')
    target_name = img_hdr1['FIELD'].strip(' ')
    exptime = img_hdr1['EXPTIME']
    filter_colour = img_hdr1['COLOUR'].strip(' ')
    obs_set = img_hdr1['SET'].strip(' ')
    chip_num = img_hdr1['CHIP']
                
    filename_to_write = "stacked-{}-{}-{}-{}-{}-{}.fit".format(run_filename,target_name,exptime,
                                                                   filter_colour,obs_set,
                                                                        chip_num)
        
    registered_img_ccd = CCDData(registered_img,unit='adu')
    registered_img_ccd.header = img_hdr1
    registered_img_ccd.write(SIS_stack_path/filename_to_write,overwrite=True)
