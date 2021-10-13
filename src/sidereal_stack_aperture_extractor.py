# -*- coding: utf-8 -*-

###############################################################################
#-------------------SECTION ZERO: IMPORTING PACKAGES--------------------------#
###############################################################################

import numpy as np
from astropy.io import fits
import astropy.units as u
import matplotlib.pyplot as plt
from photutils import DAOStarFinder, aperture_photometry
from astropy.stats import sigma_clipped_stats
from astropy.stats import SigmaClip
from photutils.background import Background2D, MedianBackground
from astropy.wcs import WCS
from photutils.datasets import load_spitzer_image, load_spitzer_catalog
from photutils.aperture import SkyCircularAperture
from astropy.coordinates import SkyCoord

from astropy.visualization import SqrtStretch, simple_norm
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils import CircularAperture, CircularAnnulus

###############################################################################
#----------------------SECTION ONE: INITIALISATION----------------------------#
###############################################################################

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++CHANGES++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Read the image
data, hdr = fits.getdata('2021-02-18-stacked-C2021_A7-300-R-a-3.fit', header=True)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++CHANGES++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#%%
# displaying image
z_min = np.percentile(data,2)
z_max = np.percentile(data,98)
plt.imshow(data,vmin=z_min,vmax=z_max)
plt.colorbar()
plt.xlim(1120,1220)
plt.ylim(2300,2400)
#%%
###############################################################################
#--------------------SECTION TWO: BACKGROUND DETECTION------------------------#
###############################################################################

# # Estimate the sky background level and its standard deviation
mean, median, std = sigma_clipped_stats(data, sigma=3.0)    

# creating a background object
sigma_clip = SigmaClip(sigma=3.)
bkg_estimator = MedianBackground()
bkg = Background2D(data, (50, 50), filter_size=(3, 3),
                   sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
bkg_median = bkg.background_rms_median

#%%
###############################################################################
#---------------------SECTION THREE: SOURCE DETECTION-------------------------#
###############################################################################

# Start up the DAOStarFinder object and detect stars
daofind = DAOStarFinder(fwhm=5.0, threshold=5.*std)    
sources = daofind(data - median)

# Print a table of what we have found
for col in sources.colnames:
    sources[col].info.format = '%.8g'  # for consistent table output
print(sources)

#%%
###############################################################################
#------------------------SECTION FOUR: APERTURES------------------------------#
###############################################################################

# creating list of sources
positions = (sources['xcentroid'], sources['ycentroid'])
the_lst = []
x_vals = sources['xcentroid']
y_vals = sources['ycentroid']
for index in range(len(x_vals)):
    to_return = (x_vals[index],y_vals[index])
    the_lst.append(to_return)

positions = (sources['xcentroid'], sources['ycentroid'])
apertures = CircularAperture(the_lst, r=20.0)
positions = the_lst

# Display the image with a "square root stretch" - this makes fainter stars show up better
z_min = np.percentile(data,2)
z_max = np.percentile(data,98)
norm = ImageNormalize(stretch=SqrtStretch())
plt.figure() #figsize=(10,10)
plt.imshow(data, cmap='Greys', origin='lower',vmin=z_min,vmax=z_max,norm=norm)
apertures.plot(color='blue', lw=1.5, alpha=0.5)

#%%
# Set up a set of circular apertures (one for each position) with a radius of 5 pixels and annuli with
# inner and outer radii of 10 and 15 pixels.
apertures = CircularAperture(positions, r=5)
annulus_apertures = CircularAnnulus(positions, r_in=10, r_out=15)

#%%
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
###############################################################################
#------------------------SECTION FIVE: PHOTOMETRY-----------------------------#
###############################################################################

# Measure the total flux in both the aperture and annulus for each star. 
apers = [apertures, annulus_apertures]
phot_table = aperture_photometry(data - bkg.background, apers, method='subpixel',
                                 subpixels=5)

# Calculate the mean flux per pixel in the annulus
sky_mean = phot_table['aperture_sum_1'] / annulus_apertures.area

# Multiply this by the number of pixels in the aperture and subtract from the aperture flux measurement.
# Put the result in a new column of the table.
aperture_sky_sum = sky_mean * apertures.area
phot_table['flux'] = phot_table['aperture_sum_0'] - aperture_sky_sum

# setting zero point and adjusting magnitudes
# Magnitude zero point is arbitrary
ZP = 30
phot_table['mag'] = ZP - 2.5*np.log10(phot_table['flux'])
print(phot_table)

#%%
# plotting stats
plt.figure()
plt.hist(phot_table['mag'],color="darkviolet")
plt.xlabel("Apparent Magnitude")
plt.ylabel("Number of Stars")
plt.title("Histogram of Apparent Magnitudes for Stars with ZP {}".format(ZP))
plt.savefig("mags_with_zp.jpeg",dpi=900)
plt.show()

plt.plot(phot_table['flux'],phot_table['mag'],'.',color="darkviolet")
plt.xlabel("Flux")
plt.ylabel("Apparent Magnitude")
plt.title("Flux vs Magnitude of Apertures")
plt.savefig("flux_vs_mag.jpeg",dpi=900)
plt.show()

plt.hist(phot_table['ycenter'].value,bins=50,color="darkviolet")
plt.xlabel("Flux")
plt.ylabel("frequency")
plt.title("Flux Values")
plt.savefig("flux_vals.jpeg",dpi=900)
plt.show()

# getting S/N ratios
xcenter_sn_ratios = []
for i in phot_table['xcenter'].value:
    xcenter_sn_ratios.append(i/np.sqrt(i))
   
ycenter_sn_ratios = []
for i in phot_table['ycenter'].value:
    ycenter_sn_ratios.append(i/np.sqrt(i))

plt.hist(xcenter_sn_ratios,label="xcenter",color="deeppink")
plt.hist(ycenter_sn_ratios,alpha=0.5,label="ycenter",color="darkviolet")
plt.title("S/N Ratios of Point Sources in Image")
plt.xlabel("S/N Ratio")
plt.ylabel("frequency")
plt.legend()
plt.savefig("sn_ratio_hist.jpeg",dpi=900)
plt.show()

#%%
###############################################################################
#-------------------------SECTION SIX: FLUXES---------------------------------#
###############################################################################

# need to do MOA-R flux conversion now
#f_MOA_R = (((0.275 * f_r) + (0.623 * f_i) + (0.094 * f_z)) * (f_g / f_i)) ** (-0.036)

phot_table_mag = phot_table['mag']
phot_table_flux = phot_table['flux']

#%%
# annulus stuff
annulus_masks = annulus_apertures.to_mask(method='center')
plt.imshow(annulus_masks[0], interpolation='nearest')
plt.colorbar()
annulus_data = annulus_masks[0].multiply(data)
mask = annulus_masks[0].data
annulus_data_1d = annulus_data[mask > 0]
print(annulus_data_1d.shape)
_, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
print(median_sigclip)  

# background subtraction stuff
bkg_median = []
for mask in annulus_masks:
    annulus_data = mask.multiply(data)
    annulus_data_1d = annulus_data[mask.data > 0]
    _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
    bkg_median.append(median_sigclip)
bkg_median = np.array(bkg_median)

# Spitzer catalog stuff
hdu = load_spitzer_image()  
data = u.Quantity(hdu.data, unit=hdu.header['BUNIT'])  
wcs = WCS(hdu.header)  
catalog = load_spitzer_catalog() 
positions = SkyCoord(catalog['l'], catalog['b'], frame='galactic')  
aperture = SkyCircularAperture(positions, r=4.8 * u.arcsec) 

# error stuff
error = 0.1 * data

#%%
phot_table = aperture_photometry(data, apertures, error=error, wcs=wcs)
phot_table['annulus_median'] = bkg_median
phot_table['aper_bkg'] = bkg_median * apertures.area
phot_table['aper_sum_bkgsub'] = phot_table['aperture_sum'] - phot_table['aper_bkg']
factor = (1.2 * u.arcsec) ** 2 / u.pixel
fluxes_catalog = catalog['f4_5']  
converted_aperture_sum = (phot_table['aperture_sum'] * factor).to(u.mJy / u.pixel) 
# for col in phot_table.colnames:
#     phot_table[col].info.format = '%.8g'  # for consistent table output
print(phot_table)

#%%















