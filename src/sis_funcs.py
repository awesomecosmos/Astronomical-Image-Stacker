# -*- coding: utf-8 -*-

# This file contains all necessary packages and functions required for the 
# Sidereal Image Stacking of the images. This file is required to run the file 
# sidereal_img_stacker.py.
# Import it as: 
# from sis_funcs import *

###############################################################################
#-------------------SECTION ONE: IMPORTING PACKAGES---------------------------#
###############################################################################

# path-type packages
import os
import glob
from pathlib import Path

# basic Python packages
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Astropy packages
import astropy
import astroalign as aa
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from astropy.stats import mad_std
from astropy.nddata import CCDData
from astropy.visualization import hist
from astropy.utils.data import get_pkg_data_filename

# ccdproc packages
import ccdproc as ccdp
from ccdproc import Combiner
from ccdproc import ImageFileCollection
from ccdproc.utils.sample_directory import sample_directory_with_files

# Misc. packages
import seaborn as sns
from datetime import datetime

# warnings
import warnings
warnings.filterwarnings('ignore')

# user-defined packages
from convenience_functions import show_image

###############################################################################
#--------------------SECTION ONE: HELPER FUNCTIONS----------------------------#
###############################################################################

def registered_img_writer(img,registered_img,SIS_stack_path,obsdate):
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
    
    obsdate : str
        Date of observation for display/saving purposes. 
        Ideal format should be "YYYY-MM-DD".
    
    Returns
    -------
    registered_img_ccd : astropy.nddata.ccddata.CCDData
        CCDData object of registered image.
    """
    img_hdul = fits.open(img)
    img_hdr1 = img_hdul[0].header
                
    run_filename = img_hdr1['RUN'].strip(' ')
    target_name = img_hdr1['FIELD'].strip(' ')
    exptime = img_hdr1['EXPTIME']
    filter_colour = img_hdr1['COLOUR'].strip(' ')
    obs_set = img_hdr1['SET'].strip(' ')
    chip_num = img_hdr1['CHIP']
                
    filename_to_write = "{}-aligned-{}-{}-{}-{}-{}-{}.fit".format(obsdate,run_filename,target_name,
                                                               exptime,filter_colour,
                                                               obs_set,chip_num)
        
    registered_img_ccd = CCDData(registered_img,unit='adu')
    registered_img_ccd.header = img_hdr1
    registered_img_ccd.write(SIS_stack_path/filename_to_write,overwrite=True)
    
    return registered_img_ccd

###############################################################################
#--------------------SECTION TWO: PLOTTING FUNCTIONS--------------------------#
###############################################################################

def img_series_stats(image_ccd_lst,plots_path,obsdate):
    """
    This function plots count statistics for a series of images.
    
    Parameters
    ----------
    image_ccd_lst : list
        List of CCDData objects.
    
    plots_path : WindowsPath object
        Path to directory where plots are to be saved.
    
    obsdate : str
        Date of observation for display/saving purposes. 
        Ideal format should be "YYYY-MM-DD".
    
    Returns
    -------
    Nothing.
    """
    median_count = []
    mean_count = []
    
    source_hdu = CCDData(image_ccd_lst[0],unit='adu')
    source_image_data = source_hdu.data.astype(float) 
    source_image_hdr = source_hdu.header
    target_name = source_image_hdr['FIELD'].strip(' ')
    exptime = source_image_hdr['EXPTIME']
    chip_num = source_image_hdr['CHIP']
        
    for a_file in image_ccd_lst:
        hdu = CCDData(a_file,unit='adu')
        image_data = hdu.data.astype(float) 
        image_hdr = hdu.header
        
        median_count.append(np.median(a_file))
        mean_count.append(np.mean(a_file))
    
    min_count_for_median = np.min(median_count)
    min_count_for_mean = np.min(mean_count)
    max_count_for_median = np.max(median_count)
    max_count_for_mean = np.max(mean_count)
    
    plt.figure()
    plt.plot(mean_count, label='mean',color="palevioletred")
    plt.axhline(y=min_count_for_mean,linestyle='-',linewidth=0.5,color='blue',label='min mean {:.2f}'.format(min_count_for_mean),alpha=1)
    plt.axhline(y=max_count_for_mean,linestyle='-',linewidth=0.5,color='blue',label='max mean {:.2f}'.format(max_count_for_mean),alpha=1)
    plt.xlabel('Image number')
    plt.ylabel('Count (ADU)')
    plt.title('Mean pixel value for aligned images')
    plt.legend()
    plt.grid()
    plt.savefig(plots_path/"{}-{}-{}-aligned_stats_mean.jpg".format(obsdate,
                                                                    target_name,
                                                                    exptime,chip_num),
                                                                    dpi=900)
    plt.show()

    plt.figure()
    plt.plot(median_count, label='median',color="darkviolet")
    plt.axhline(y=min_count_for_median,linestyle='-',linewidth=0.5,color='red',label='min median {:.2f}'.format(min_count_for_median),alpha=1)
    plt.axhline(y=max_count_for_median,linestyle='-',linewidth=0.5,color='red',label='max median {:.2f}'.format(max_count_for_median),alpha=1)               
    plt.xlabel('Image number')
    plt.ylabel('Count (ADU)')
    plt.title('Median pixel value for aligned images')
    plt.legend()
    plt.grid()
    plt.savefig(plots_path/"{}-{}-{}-aligned_stats_median.jpg".format(obsdate,
                                                                      target_name,
                                                                      exptime,chip_num),
                                                                      dpi=900)
    plt.show()

def aligned_comparison_stats(unaligned_image_ccd_lst,aligned_image_ccd_lst,
                             plots_path,obsdate):
    """
    This function plots count statistics for a series of images to compare
    non-aligned and aligned images.
    
    Parameters
    ----------
    unaligned_image_ccd_lst : list
        List of CCDData objects of non-aligned images.
    
    aligned_image_ccd_lst : list
        List of CCDData objects of aligned images.
    
    plots_path : WindowsPath object
        Path to directory where plots are to be saved.
    
    obsdate : str
        Date of observation for display/saving purposes. 
        Ideal format should be "YYYY-MM-DD".
    
    Returns
    -------
    Nothing.
    """
    aligned_mean_count = []
    aligned_median_count = []
    unaligned_mean_count = []
    unaligned_median_count = []
    
    # getting information from source image for saving purposes
    source_hdu = CCDData(aligned_image_ccd_lst[0],unit='adu')
    source_image_data = source_hdu.data.astype(float) 
    source_image_hdr = source_hdu.header
    target_name = source_image_hdr['FIELD'].strip(' ')
    exptime = source_image_hdr['EXPTIME']
    chip_num = source_image_hdr['CHIP']
    
    for a_file in unaligned_image_ccd_lst[1:]:
        hdu = CCDData(a_file,unit='adu')
        image_data = hdu.data.astype(float) 
        image_hdr = hdu.header
        
        unaligned_mean_count.append(np.mean(a_file))
        unaligned_median_count.append(np.median(a_file))
    
    min_count_for_unaligned_mean = np.min(unaligned_mean_count)
    max_count_for_unaligned_mean = np.max(unaligned_mean_count)
    min_count_for_unaligned_median = np.min(unaligned_median_count)
    max_count_for_unaligned_median = np.max(unaligned_median_count)
        
    for a_file in aligned_image_ccd_lst:
        hdu = CCDData(a_file,unit='adu')
        image_data = hdu.data.astype(float) 
        image_hdr = hdu.header
        
        aligned_mean_count.append(np.mean(a_file))
        aligned_median_count.append(np.median(a_file))
    
    min_count_for_aligned_mean = np.min(aligned_mean_count)
    max_count_for_aligned_mean = np.max(aligned_mean_count)
    min_count_for_aligned_median = np.min(aligned_median_count)
    max_count_for_aligned_median = np.max(aligned_median_count)
    
    # plotting
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15,8),sharex=True, sharey=False)
    fig.suptitle('Comparison Statistics for Non-Aligned Images vs Aligned Images for {}'.format(obsdate))
    
    # unaligned mean
    ax1.plot(unaligned_mean_count, label='mean',color="deeppink")
    ax1.axhline(y=min_count_for_unaligned_mean,linestyle='-',linewidth=1,
                color='black',label='min mean {:.2f}'.format(min_count_for_unaligned_mean),alpha=1)
    ax1.axhline(y=max_count_for_unaligned_mean,linestyle='-',linewidth=1,
                color='black',label='max mean {:.2f}'.format(max_count_for_unaligned_mean),alpha=1)
    ax1.set_xlabel('Image number')
    ax1.set_ylabel('Count (ADU)')
    ax1.set_title('Mean pixel value for unaligned images')
    ax1.legend(loc="best")
    ax1.grid(b=True, which='both', axis='both')
    # ax1.sharex(ax3)
    
    # unaligned median
    ax2.plot(unaligned_median_count, label='median',color="deeppink")
    ax2.axhline(y=min_count_for_unaligned_median,linestyle='-',linewidth=1,
                color='black',label='min median {:.2f}'.format(min_count_for_unaligned_median),alpha=1)
    ax2.axhline(y=max_count_for_unaligned_median,linestyle='-',linewidth=1,
                color='black',label='max median {:.2f}'.format(max_count_for_unaligned_median),alpha=1)
    ax2.set_xlabel('Image number')
    ax2.set_ylabel('Count (ADU)')
    ax2.set_title('Median pixel value for unaligned images')
    ax2.legend(loc="best")
    ax2.grid(b=True, which='both', axis='both')
    # ax2.sharex(ax4)
    
    # aligned mean
    ax3.plot(aligned_mean_count, label='mean',color="darkviolet")
    ax3.axhline(y=min_count_for_aligned_mean,linestyle='-',linewidth=1,
                color='black',label='min mean {:.2f}'.format(min_count_for_aligned_mean),alpha=1)
    ax3.axhline(y=max_count_for_aligned_mean,linestyle='-',linewidth=1,
                color='black',label='max mean {:.2f}'.format(max_count_for_aligned_mean),alpha=1)
    ax3.set_xlabel('Image number')
    ax3.set_ylabel('Count (ADU)')
    ax3.set_title('Mean pixel value for aligned images')
    ax3.legend(loc="best")
    ax3.grid(b=True, which='both', axis='both')
    
    
    # aligned median
    ax4.plot(aligned_median_count, label='median',color="darkviolet")
    ax4.axhline(y=min_count_for_aligned_median,linestyle='-',linewidth=1,
                color='black',label='min median {:.2f}'.format(min_count_for_aligned_median),alpha=1)
    ax4.axhline(y=max_count_for_aligned_median,linestyle='-',linewidth=1,
                color='black',label='max median {:.2f}'.format(max_count_for_aligned_median),alpha=1)
    ax4.set_xlabel('Image number')
    ax4.set_ylabel('Count (ADU)')
    ax4.set_title('Median pixel value for unaligned images')
    ax4.legend(loc="best")
    ax4.grid(b=True, which='both', axis='both')
    
    
    for ax in fig.get_axes():
        ax.label_outer()
    
    fig.savefig(plots_path/"{}-{}-{}-comparison_stats.jpg".format(obsdate,
                                                                  target_name,
                                                                  exptime,chip_num),
                                                                    dpi=1000)
    fig.show()

def image_comparison(unaligned_image_ccd_lst,aligned_image_ccd_lst,stacked_img_ccd,outputs_path,obsdate):
    """
    This function produces plots of unaligned images vs aligned images, and of
    the source image vs the stacked image.
    
    Parameters
    ----------
    unaligned_image_ccd_lst : list
        List of CCDData objects of non-aligned images.
    
    aligned_image_ccd_lst : list
        List of CCDData objects of aligned images.
    
    stacked_img_ccd : astropy.nddata.ccddata.CCDData
        CCDData object of stacked image.
    
    outputs_path : WindowsPath object
        Path to directory where plots are to be saved.
    
    obsdate : str
        Date of observation for display/saving purposes. 
        Ideal format should be "YYYY-MM-DD".
    
    Returns
    -------
    Nothing.
    """
    source_hdu = CCDData(unaligned_image_ccd_lst[0],unit='adu')
    source_image_hdr = source_hdu.header
    run_filename = source_image_hdr['RUN'].strip(' ')
    target_name = source_image_hdr['FIELD'].strip(' ')
    exptime = source_image_hdr['EXPTIME']
    chip_num = source_image_hdr['CHIP']
    
    # compare unaligned vs aligned images
    for i, unaligned_img in enumerate(unaligned_image_ccd_lst[1:]):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10), tight_layout=True)
        
        # source_hdu = CCDData(unaligned_image_ccd_lst[0],unit='adu')
        image_hdr = unaligned_img.header
        run_filename = image_hdr['RUN'].strip(' ')
        target_name = image_hdr['FIELD'].strip(' ')
        exptime = image_hdr['EXPTIME']
        chip_num = image_hdr['CHIP']
        
        show_image(unaligned_img, cmap='gray', ax=ax1, fig=fig, percl=90)
        ax1.set_title('Unaligned Image for {}-{}-{}-{}s ({})'.format(run_filename,target_name,chip_num,exptime,obsdate))
        
        show_image(aligned_image_ccd_lst[i], cmap='gray', ax=ax2, fig=fig, percl=90)
        ax2.set_title('Aligned Image for {}-{}-{}-{}s ({})'.format(run_filename,target_name,chip_num,exptime,obsdate))
        
        plt.savefig(outputs_path/"unaligned_vs_aligned_{}-{}-{}-{}.jpg".format(run_filename,target_name,chip_num,exptime),dpi=900)
        plt.show()
    
    # compare source image to stacked image
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10), tight_layout=True)
    
    show_image(unaligned_image_ccd_lst[0], cmap='gray', ax=ax1, fig=fig, percl=90)
    ax1.set_title('Unaligned Source Image for {}-{}-{}s ({})'.format(target_name,chip_num,exptime,obsdate))
    
    show_image(stacked_img_ccd, cmap='gray', ax=ax2, fig=fig, percl=90)
    ax2.set_title('Aligned Stacked Image for {}-{}-{}s ({})'.format(target_name,chip_num,exptime,obsdate))
    
    plt.savefig(outputs_path/"source_vs_stacked_{}-{}-{}.jpg".format(target_name,chip_num,exptime),dpi=900)
    plt.show()
    
###############################################################################
#-------------------------------END OF CODE-----------------------------------#
###############################################################################