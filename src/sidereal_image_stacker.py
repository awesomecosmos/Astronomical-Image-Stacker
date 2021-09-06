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
sis_code_path = "C:/Users/ave41/OneDrive - University of Canterbury/Master's 2021/ASTR480 Research/ASTR480 Code/03 Photometrical Analysis"
os.chdir(sis_code_path) #from now on, we are in this directory

# importing functions
from sis_funcs import *

#%%
# defining paths/reading in data files
test_data_path = "//spcsfs/ave41/astro/ave41/SIS_TestData_v1"
data_path = "//spcsfs/ave41/astro/ave41/ObsData-2021-02-18/ALERT/Reduced ALERT/WCS Calibrated"

SIS_path = Path(data_path)
SIS_stack_path = path_checker(SIS_path,'Sidereally Stacked Images')
outputs_path = path_checker(SIS_stack_path,'Outputs')

#%%
lst_of_files = [n for n in os.listdir(SIS_path) if (n.endswith('fit'))]

data_lst = []
for i in lst_of_files:
    data_lst.append(os.path.join(SIS_path,i))

#%%
source_img = data_lst[0]
source_img_hdu = fits.open(source_img)
source_img_data = source_img_hdu[0].data
    
for i in range(len(data_lst)):
    target_img = data_lst[i]
    target_img_hdu = fits.open(target_img)
    target_img_data = target_img_hdu[0].data
    
    if not source_img_data.dtype.byteorder == '<':
        source_img_data2 = source_img_data.byteswap().newbyteorder()
        target_img_data2 = target_img_data.byteswap().newbyteorder()
    
        registered_image, footprint = aa.register(source_img_data2, target_img_data2)
        registered_img_writer(target_img,registered_image,SIS_stack_path)

    else:
        registered_image, footprint = aa.register(source_img_data, target_img_data)
        registered_img_writer(target_img,registered_image,SIS_stack_path)





























