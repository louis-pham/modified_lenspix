# rotate map using code from george -- used to debug cib lensing; meant to break correlation between primary and kappa
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits 

fileName = '/scratch2/r/rbond/phamloui/lenspix_files/kappa_maps/8Gpc_n2048_nb23_nt18_kap_comb.fits'
y = hp.read_map(fileName)

#get nside                                                                                                                 
hdulist = fits.open(fileName)
nside = hdulist[1].header['NSIDE']
hdulist.close()

pix = np.arange(len(y))
theta,phi = hp.pix2ang(nside,pix)
thetam = np.pi - theta
pixnew = hp.ang2pix(nside,thetam,phi)
y = y[pixnew]

hp.write_map("/scratch2/r/rbond/phamloui/lenspix_files/output/rotated_maps/8Gpc_n2048_nb23_nt18_kap_comb.fits", y)
