import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from matplotlib import cm
from astropy.io import fits

kappaFile = "/scratch2/r/rbond/phamloui/lenspix_files/cl_kappa_restricted.txt"
primaryFile = "/scratch2/r/rbond/phamloui/lenspix_files/FromNERSC/cls/ffp10_scalCls.dat"
kappaMapFile = "/scratch2/r/rbond/phamloui/lenspix_files/kappa_maps/map_kappa_restricted.fits"
primaryMapFile = "/scratch2/r/rbond/phamloui/lenspix_files/FromNERSC/ffp10_map_from_scalCls.fits"
primaryAlmFile = "/scratch2/r/rbond/phamloui/lenspix_files/FromNERSC/ffp10_alm_from_scalCls.fits"
nside = 2048

kappaElls, kappaCl = np.loadtxt(kappaFile, usecols=(0,1), unpack=True)
primaryElls,primaryCl = np.loadtxt(primaryFile, usecols=(0,1), unpack=True)

primaryCl = primaryCl * 2 * np.pi / (primaryElls * (primaryElls+1)) / 1e12

primaryElls = np.insert(primaryElls, 0, 1)
primaryElls = np.insert(primaryElls, 0, 0)
primaryCl = np.insert(primaryCl, 0, 0)
primaryCl = np.insert(primaryCl, 0, 0)

print primaryElls
print primaryCl

kappaMap = hp.synfast(kappaCl, nside)
primaryMap, primaryAlm = hp.synfast(primaryCl, nside, alm=True)

#print kappaCl.shape
#print primaryCl.shape
hot_cmap = cm.get_cmap("hot")
hot_cmap.set_under("w")

hp.mollview(kappaMap,xsize=2048, cmap=hot_cmap)
hp.mollview(primaryMap,xsize=2048, cmap=hot_cmap)

plt.show()

hp.write_map(kappaMapFile, kappaMap)
hp.write_map(primaryMapFile, primaryMap)
hp.write_alm(primaryAlmFile, primaryAlm)
