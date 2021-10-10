# -*- coding: utf-8 -*-
"""
Created on Mon Sep 27 13:32:33 2021

@author: ave41
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from photutils import DAOStarFinder, aperture_photometry
from astropy.stats import sigma_clipped_stats

from astropy.visualization import SqrtStretch, simple_norm
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils import CircularAperture, CircularAnnulus


# Read the g-band image
data, hdr = fits.getdata('2021-02-18-stacked-C2021_A7-300-R-a-3.fit', header=True)
#%%
z_min = np.percentile(data,2)
z_max = np.percentile(data,98)
plt.imshow(data,vmin=z_min,vmax=z_max)
plt.colorbar()
plt.xlim(1120,1220)
plt.ylim(2300,2400)
#%%

# Estimate the sky background level and its standard deviation
mean, median, std = sigma_clipped_stats(data, sigma=3.0)    
print('Mean, median, std =',(mean, median, std))
print()

# Start up the DAOStarFinder object and detect stars
#     sources is an astropy table - see http://docs.astropy.org/en/stable/table/index.html
#     We try to find all stars with widths of about 5 pixels and peak values more than 
#     5 x the standard deviation of the background

daofind = DAOStarFinder(fwhm=5.0, threshold=5.*std)    
sources = daofind(data - median)

# Print a table of what we have found
for col in sources.colnames:
    sources[col].info.format = '%.8g'  # for consistent table output
print(sources)

#%%
positions = (sources['xcentroid'], sources['ycentroid'])

the_lst = []
x_vals = sources['xcentroid']
y_vals = sources['ycentroid']
for index in range(len(x_vals)):
    to_return = (x_vals[index],y_vals[index])
    the_lst.append(to_return)

# Display the image with a "square root stretch" - this makes fainter stars show up better
z_min = np.percentile(data,2)
z_max = np.percentile(data,98)
norm = ImageNormalize(stretch=SqrtStretch())
plt.figure()#figsize=(10,10)
plt.imshow(data, cmap='Greys', origin='lower',vmin=z_min,vmax=z_max,norm=norm)

# Plot some circles (apertures) on top of the image at the locations
# of the detected stars. Note that 'xcentroid' and 'ycentroid' are columns
# in our sources table.
positions = (sources['xcentroid'], sources['ycentroid'])
apertures = CircularAperture(the_lst, r=20.0)
apertures.plot(color='blue', lw=1.5, alpha=0.5)

#%%
# Display the image with a "square root stretch" - this makes fainter stars show up better
z_min = np.percentile(data,2)
z_max = np.percentile(data,98)
norm = ImageNormalize(stretch=SqrtStretch())
plt.figure()#figsize=(10,10)
plt.imshow(data, cmap='Greys', origin='lower',vmin=z_min,vmax=z_max,norm=norm)

# Plot some circles (apertures) on top of the image at the locations
# of the detected stars. Note that 'xcentroid' and 'ycentroid' are columns
# in our sources table.
positions = (sources['xcentroid'], sources['ycentroid'])
apertures = CircularAperture(the_lst, r=20.0)
apertures.plot(color='blue', lw=1.5, alpha=0.5)
positions = the_lst

#%%
# Set up a set of circular apertures (one for each position) with a radius of 5 pixels and annuli with
# inner and outer radii of 10 and 15 pixels.

apertures = CircularAperture(positions, r=5)
annulus_apertures = CircularAnnulus(positions, r_in=10, r_out=15)

# Display the image and overlay the apertures and anulli - zoom in to check that they are OK 
plt.figure()
norm = simple_norm(data, 'sqrt', percent=96) # another way to "stretch" the image display
plt.imshow(data, norm=norm,origin='lower')
plt.colorbar()
apertures.plot(color='deeppink', lw=2)
annulus_apertures.plot(color='red', lw=2)
# plt.xlim(1500, 1600)
# plt.ylim(2000, 2100)
plt.savefig("apertures.jpeg",dpi=900)
plt.show()

plt.figure()
norm = simple_norm(data, 'sqrt', percent=96) # another way to "stretch" the image display
plt.imshow(data, norm=norm,origin='lower')
plt.colorbar()
apertures.plot(color='deeppink', lw=2)
annulus_apertures.plot(color='red', lw=2)
plt.xlim(1120,1220)
plt.ylim(2300,2400)
plt.savefig("trailed_comet_aperture.jpeg",dpi=900)
plt.show()

#%%
# Measure the total flux in both the aperture and annulus for each star. 
apers = [apertures, annulus_apertures]
phot_table = aperture_photometry(data, apers)

print(phot_table)

# Calculate the mean flux per pixel in the annulus
sky_mean = phot_table['aperture_sum_1'] / annulus_apertures.area

# Multiply this by the number of pixels in the aperture and subtract from the aperture flux measurement.
# Put the result in a new column of the table.
aperture_sky_sum = sky_mean * apertures.area
phot_table['flux'] = phot_table['aperture_sum_0'] - aperture_sky_sum

print(phot_table)

#%%
# Magnitude zero point is arbitrary
ZP = 25
phot_table['mag'] = ZP - 2.5*np.log10(phot_table['flux'])

print(phot_table)

plt.figure()
plt.hist(phot_table['mag'],color="darkviolet")
plt.xlabel("Apparent Magnitude")
plt.ylabel("Number of Stars")
plt.title("Histogram of Apparent Magnitudes for Stars with ZP {}".format(ZP))
plt.savefig("mags_with_zp.jpeg",dpi=900)
plt.show()

#%%

plt.plot(phot_table['flux'],phot_table['mag'],'.',color="darkviolet")
plt.xlabel("Flux")
plt.ylabel("Apparent Magnitude")
plt.title("Flux vs Magnitude of Apertures")
plt.savefig("flux_vs_mag.jpeg",dpi=900)
plt.show()

#%%

plt.hist(phot_table['ycenter'].value,bins=50,color="darkviolet")
plt.xlabel("Flux")
# plt.ylabel("Apparent Magnitude")
plt.title("Flux vs Magnitude")
# plt.savefig("flux_vs_mag.jpeg",dpi=900)
plt.show()
#%%

xcenter_sn_ratios = []
for i in phot_table['xcenter'].value:
    xcenter_sn_ratios.append(i/np.sqrt(i))
   
ycenter_sn_ratios = []
for i in phot_table['ycenter'].value:
    ycenter_sn_ratios.append(i/np.sqrt(i))

plt.hist(xcenter_sn_ratios,label="xcenter")
plt.hist(ycenter_sn_ratios,alpha=0.5,label="ycenter")
plt.title("S/N Ratios of Point Sources in Image")
plt.xlabel("S/N Ratio")
plt.ylabel("frequency")
plt.legend()
plt.savefig("sn_ratio_hist.jpeg",dpi=900)
plt.show()









