import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import search_around_sky
from astropy import units as u
from astropy.table import Table
import pandas as pd

import time
import pickle
from scipy import spatial
#np.set_printoptions(threshold=np.nan)

#####################################################
"""
# input
folder = int(input("What folder do you want this stored in? (1) For Color Plot, (2) For Error Plot"))
if folder == 1:
    folder = str("ColorPlots")
else:
    folder = str("ErrorPlots")
general_title = str(input("What do you want to title the plots?"))
title_cmodel = str(input("What do you want to name this CMODEL plot png?"))
title_psfmodel = str(input("What do you want to name this PSF plot png?"))

#error = float(input("What error value do you want to consider for this plot?"))
#every_num_of_points = int(input("Every how many points do you want to plot?"))

"""
#####################################################
# data
data = Table.read('Chandra_COSMOS_Legacy_20151120_4d.fits', format='fits')

Chandra_Data = data.to_pandas()

data_COSMOS_DeepDepth = pd.read_csv("Color_COSMOS_DD.csv", dtype = np.float64)


HSC_Data = data_COSMOS_DeepDepth.dropna(axis = 0, how = "any")

#Getting RA and Dec
HSC_RA = HSC_Data["# ra"]
HSC_Dec = HSC_Data["dec"]

Chandra_RA = Chandra_Data["RA_x"]
Chandra_Dec = Chandra_Data["DEC_x"]

"""
# Making the cross match using astropy's SkyCoord
Chandra_catalog = SkyCoord(ra=Chandra_RA*u.degree, dec=Chandra_Dec*u.degree)
HSC_catalog = SkyCoord(ra=HSC_RA*u.degree, dec=HSC_Dec*u.degree)


with open('HSC_skycoords.pkl', 'wb') as f:
    pickle.dump(HSC_catalog, f)
with open('Chandra_skycoords.pkl', 'wb') as f:
    pickle.dump(Chandra_catalog, f)
"""
#####################################################
#getting catalogs in sky coords
with open('HSC_skycoords.pkl', 'rb') as f:
    HSC_catalog = pickle.load(f)
with open('Chandra_skycoords.pkl', 'rb') as f:
    Chandra_catalog = pickle.load(f)

#getting matches
max_sep = 1.0 * u.arcsec
idx1, idx2, sep2d, dist3d = search_around_sky(Chandra_catalog, HSC_catalog, max_sep)
number_matches = Chandra_RA[idx1].shape


#combining both catalogs based on matches
#matches_df = Chandra_Data["RA_x"].iloc[idx1].join(HSC_Data["# ra"].iloc[idx2])


#matches g,r, & i in cmodel and psf, as well as flux_f and hardness ratio
matches_cmodel_gr = HSC_Data["gcmodel_mag"][idx2] - HSC_Data["rcmodel_mag"][idx2]
matches_cmodel_ri = HSC_Data["rcmodel_mag"][idx2] - HSC_Data["icmodel_mag"][idx2]
matches_psf_gr =  HSC_Data["gmag_psf"][idx2] -  HSC_Data["rmag_psf"][idx2]
matches_psf_ri =  HSC_Data["rmag_psf"][idx2] -  HSC_Data["imag_psf"][idx2]
matches_flux_f = Chandra_Data["flux_F"][idx1]
matches_hardness_ratio = Chandra_Data["HR"][idx1]


"""
#####################################################
#plot RA & Dec of Chandra matches v all Chandra data
plt.style.use("seaborn")
plt.title("Chandra Matches v All Chandra Data", weight = "bold", size = 14)
plt.scatter(Chandra_RA, Chandra_Dec, label = "All Chandra Data")
plt.scatter(Chandra_RA[idx1], Chandra_Dec[idx1], alpha = 0.6, label = "Chandra Matches: \n %i matches" % (number_matches[0]))
plt.xlim(149.25,151.75)
plt.xlabel("Right Ascension (J2000)", weight = "bold", size = 12)
plt.ylabel("Declination (J2000)", weight = "bold", size = 12)
plt.legend(loc = "best", prop={'size': 12})
#plt.savefig("/Users/linaflorez/Desktop/ObsCos/COSMOS_research/Plots/Matches/Chandra_RAvDec.pdf")
plt.show()

#####################################################

#plot histogram of flux_F
plt.title("Chandra Matches v All Chandra Data", weight = "bold", size = 14)
plt.hist(Chandra_Data["flux_F"], bins = 80, label = "All Chandra Data")
plt.hist(Chandra_Data["flux_F"][idx1], bins = 80, alpha = 0.5, label = "Chandra Matches: \n %i matches" % (number_matches[0]))
plt.xlabel("Flux F (0.5-10 keV flux [erg/s/cm${^2}$])", weight = "bold", size = 12)
plt.legend(loc = "upper right", prop={'size': 12})
plt.yscale("log")
#plt.savefig("/Users/linaflorez/Desktop/ObsCos/COSMOS_research/Plots/Matches/Chandra_fluxFhist.pdf")
plt.show()

#####################################################
"""

matches_cmodel_gr = matches_cmodel_gr.values
matches_flux_f = matches_flux_f.values
mask = (np.isfinite(matches_cmodel_gr)) & (np.isfinite(matches_flux_f))

"""
plt.title("g-r v flux_F", weight = "bold", size = 14)
plt.semilogy(matches_cmodel_gr, matches_flux_f, 'k.')
fig.colorbar(im1, cax=cax, orientation='vertical')plt.xlabel("g-r", weight = "bold", size = 12)
plt.ylabel("Flux F (0.5-10 keV flux [erg/s/cm${^2}$])", weight = "bold", size = 12)
#plt.savefig("/Users/linaflorez/Desktop/ObsCos/COSMOS_research/Plots/Matches/Chandra_grvfluxF.pdf")
plt.show()
 
#####################################################

#plot r-i v flux_F
plt.title("r-i v flux_F", weight = "bold", size = 14)
plt.semilogy(matches_cmodel_ri, matches_flux_f, 'k.')
plt.xlabel("r-i", weight = "bold", size = 12)
plt.ylabel("Flux F (0.5-10 keV flux [erg/s/cm${^2}$])", weight = "bold", size = 12)
#plt.savefig("/Users/linaflorez/Desktop/ObsCos/COSMOS_research/Plots/Matches/Chandra_rivfluxF.pdf")
plt.show()


#####################################################

#plot g-r v hardness ratio
plt.title("g-r v flux_F", weight = "bold", size = 14)
plt.semilogy(matches_cmodel_gr, matches_hardness_ratio, 'k.')
plt.xlabel("g-r", weight = "bold", size = 12)
plt.ylabel("Hardness Ratio (H-S)/(H+S)", weight = "bold", size = 12)
#plt.savefig("/Users/linaflorez/Desktop/ObsCos/COSMOS_research/Plots/Matches/Chandra_grvHR.pdf")
plt.show()


#####################################################

#plot r-i v hardness ratio
plt.title("r-i v flux_F", weight = "bold", size = 14)
plt.semilogy(matches_cmodel_ri, matches_hardness_ratio, 'k.')
plt.xlabel("r-i", weight = "bold", size = 12)
plt.ylabel("Hardness Ratio (H-S)/(H+S)", weight = "bold", size = 12)
#plt.savefig("/Users/linaflorez/Desktop/ObsCos/COSMOS_research/Plots/Matches/Chandra_rivHR.pdf")
plt.show()
"""

