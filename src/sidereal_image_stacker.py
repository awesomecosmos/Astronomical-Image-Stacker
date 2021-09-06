# -*- coding: utf-8 -*-

# SIS = Sidereal Image Stacker
# NIS = Non-sidereal Image Stacker

###############################################################################
#-------------------SECTION ZERO: IMPORTING PACKAGES--------------------------#
###############################################################################

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
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv CHANGES vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# defining paths
# test_data_path = "//spcsfs/ave41/astro/ave41/SIS_TestData_v1"
data_path = "//spcsfs/ave41/astro/ave41/ObsData-2021-02-18/ALERT/Reduced ALERT/WCS Calibrated"

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ CHANGES ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#%%
# defining paths
SIS_path = Path(data_path)
SIS_stack_path = path_checker(SIS_path,'Sidereally Stacked Images')
outputs_path = path_checker(SIS_stack_path,'Outputs')

# reading in data files
lst_of_files = [n for n in os.listdir(SIS_path) if (n.endswith('fit'))]
data_lst = []
for i in lst_of_files:
    data_lst.append(os.path.join(SIS_path,i))

#%%
# aligning images
source_img = data_lst[0]
source_img_hdu = fits.open(source_img)
source_img_data = source_img_hdu[0].data
  
registered_image_ccd_lst = []
  
for i in range(len(data_lst)):
    target_img = data_lst[i]
    target_img_hdu = fits.open(target_img)
    target_img_data = target_img_hdu[0].data
    
    # checking for endian-ness of data
    # my machine is a little-endian compiler so need to convert from big-endian
    if not source_img_data.dtype.byteorder == '<':
        source_img_data2 = source_img_data.byteswap().newbyteorder()
        target_img_data2 = target_img_data.byteswap().newbyteorder()
    
        registered_image, footprint = aa.register(source_img_data2, target_img_data2)
        registered_image_ccd = registered_img_writer(target_img,registered_image,SIS_stack_path)
        registered_image_ccd_lst.append(registered_image_ccd)

    else:
        registered_image, footprint = aa.register(source_img_data, target_img_data)
        registered_image_ccd = registered_img_writer(target_img,registered_image,SIS_stack_path)
        registered_image_ccd_lst.append(registered_image_ccd)

#%%
registered_image_ccd_data = []
for i in registered_image_ccd_lst:
    registered_image_ccd_data.append(i.data)

# plotting stats for each aligned image
aligned_img_stats(registered_image_ccd_data,outputs_path)
#%%
source_img_hdr1 = source_img_hdu[0].header

target_name = source_img_hdr1['FIELD'].strip(' ')
exptime = source_img_hdr1['EXPTIME']
filter_colour = source_img_hdr1['COLOUR'].strip(' ')
obs_set = source_img_hdr1['SET'].strip(' ')
chip_num = source_img_hdr1['CHIP']

filename_to_write = "stacked-{}-{}-{}-{}-{}.fit".format(target_name,exptime,
                                                        filter_colour,obs_set,
                                                        chip_num)

stacked_img = sum(registered_image_ccd_data) / len(registered_image_ccd_data)
stacked_img_ccd = CCDData(stacked_img,unit='adu')
stacked_img_ccd.header = source_img_hdr1
       
stacked_img_ccd.write(SIS_stack_path/filename_to_write,overwrite=True)

#%%
plt.figure()
sns.set_theme(style="whitegrid")
ax = sns.violinplot(x=stacked_img.flatten(),color="darkviolet")
ax.set_title('Distribution of counts for stacked image')
ax.set_xlabel('Counts')
plt.savefig(outputs_path/"stacked_stats_violin.jpg",dpi=900)
plt.show()



























