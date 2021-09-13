# -*- coding: utf-8 -*-

# SIS = Sidereal Image Stacker
# NIS = Non-sidereal Image Stacker

###############################################################################
#-------------------SECTION ZERO: IMPORTING PACKAGES--------------------------#
###############################################################################

# initialising timer so we can count elapsed time
from pytictoc import TicToc
t = TicToc() # create TicToc instance
t.tic() # Start timer

import os

# initialising starting directory
drp_code_path = "C:/Users/ave41/OneDrive - University of Canterbury/Master's 2021/ASTR480 Research/ASTR480 Code/01 Data Reduction Pipeline/DataReductionPipeline/src"
os.chdir(drp_code_path) #from now on, we are in this directory

# importing functions
from drp_funcs import *
from asp_funcs import *

# initialising starting directory
sis_code_path = "C:/Users/ave41/OneDrive - University of Canterbury/Master's 2021/ASTR480 Research/ASTR480 Code/02 Data Analysis/Image-Stacker/src"
os.chdir(sis_code_path) #from now on, we are in this directory

# importing functions
from sis_funcs import *

#%%
###############################################################################
#----------------------SECTION ONE: INITIALISATION----------------------------#
###############################################################################


# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv CHANGES vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# defining observation date for display/saving purposes
obsdate = "2021-01-21"

# "\\spcsfs\ave41\astro\ave41\C2021_A7"

# data_path = "//spcsfs/ave41/astro/ave41/C2021_A7"
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ CHANGES ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#%%
# defining paths
# test_data_path = "//spcsfs/ave41/astro/ave41/SIS_TestData_v1"
data_path = "//spcsfs/ave41/astro/ave41/ObsData-{}/ALERT/Reduced ALERT/WCS Calibrated".format(obsdate)

# defining paths
SIS_path = Path(data_path)
SIS_stack_path = path_checker(SIS_path,'Sidereally Stacked Images')
outputs_path = path_checker(SIS_stack_path,'Outputs')

# reading in data files
lst_of_files = [n for n in os.listdir(SIS_path) if (n.endswith('fit'))]
data_lst = []
for i in lst_of_files:
    data_lst.append(os.path.join(SIS_path,i))

# making list of CCDData objects of these non-aligned images
data_lst_images = []
for i in data_lst:
    img = fits.open(i)
    img_hdr = img[0].header
    data_array = CCDData(img[0].data,unit='adu')
    data_array.header = img_hdr
    data_lst_images.append(data_array)

#%%
###############################################################################
#---------------------SECTION TWO: ALIGNING IMAGES----------------------------#
###############################################################################

# aligning/registering images
source_img = data_lst[0]
source_img_hdu = fits.open(source_img)
source_img_data = source_img_hdu[0].data
  
registered_image_ccd_lst = []
  
for i in range(len(data_lst[1:])):
    target_img = data_lst[i+1]
    target_img_hdu = fits.open(target_img)
    target_img_data = target_img_hdu[0].data
    
    # checking for endian-ness of data
    # my machine is a little-endian compiler so need to convert from big-endian
    if not source_img_data.dtype.byteorder == '<':
        source_img_data2 = source_img_data.byteswap().newbyteorder()
        target_img_data2 = target_img_data.byteswap().newbyteorder()
        
        # the first argument of aa.register matches to the second argument
        # so the second argument is the 'reference' image to which the first image
        # is being matched 
        registered_image, footprint = aa.register(target_img_data2, source_img_data2)
        # registered_image, footprint = aa.register(source_img_data2,target_img_data2)
        registered_image_ccd = registered_img_writer(target_img,registered_image,SIS_stack_path,obsdate)
        registered_image_ccd_lst.append(registered_image_ccd)

    else:
        registered_image, footprint = aa.register(target_img_data, source_img_data)
        # registered_image, footprint = aa.register(source_img_data,target_img_data)
        registered_image_ccd = registered_img_writer(target_img,registered_image,SIS_stack_path,obsdate)
        registered_image_ccd_lst.append(registered_image_ccd)

registered_image_ccd_data = []
for i in registered_image_ccd_lst:
    registered_image_ccd_data.append(i.data)

#%%
# plotting stats for comapring non-aligned vs aligned images
aligned_comparison_stats(data_lst_images,registered_image_ccd_lst,outputs_path,obsdate)

#%%
###############################################################################
#---------------------SECTION THREE: STACKING IMAGES--------------------------#
###############################################################################

# defining source image and extracting information from source image
source_img_hdr1 = source_img_hdu[0].header
target_name = source_img_hdr1['FIELD'].strip(' ')
exptime = source_img_hdr1['EXPTIME']
filter_colour = source_img_hdr1['COLOUR'].strip(' ')
obs_set = source_img_hdr1['SET'].strip(' ')
chip_num = source_img_hdr1['CHIP']

filename_to_write = "{}-stacked-{}-{}-{}-{}-{}.fit".format(obsdate,target_name,exptime,
                                                        filter_colour,obs_set,
                                                        chip_num)
plot_filename_to_write = "{}-stacked-stats-violin-{}-{}-{}-{}-{}.jpg".format(obsdate,target_name,exptime,
                                                        filter_colour,obs_set,
                                                        chip_num)
# stacking aligned images together
stacked_img = sum(registered_image_ccd_data) 
stacked_img_ccd = CCDData(stacked_img,unit='adu')
stacked_img_ccd.header = source_img_hdr1
stacked_img_ccd.write(SIS_stack_path/filename_to_write,overwrite=True)

#%%
# violin plot of stacked image counts
plt.figure()
sns.set_theme(style="whitegrid")
ax = sns.violinplot(x=stacked_img.flatten(),color="darkviolet")
ax.set_title('Distribution of counts for stacked image')
ax.set_xlabel('Counts')
plt.savefig(outputs_path/plot_filename_to_write,dpi=900)
plt.show()

#%%
# comparing aligned vs unaligned, and source vs stacked images
image_comparison(data_lst_images,registered_image_ccd_lst,stacked_img_ccd,outputs_path,obsdate)

#%%
#================================ don't touch ================================#

###############################################################################
#-------------------------------END OF CODE-----------------------------------#
###############################################################################

t.toc() # Print elapsed time

###############################################################################
#-------------------------------END OF CODE-----------------------------------#
###############################################################################