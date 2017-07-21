import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from matplotlib import cm
import os.path

def gaussFunc(ell, ell_o):
    return np.exp(-(ell ** 2 / ell_o ** 2))

def appendToFilename(filename, newStr):
    return "{0}_{2}.{1}".format(*filename.rsplit('.', 1) + [newStr])

zInit=0.0
zFinal=4.6
n=22
dz=0.2
i=1
nside = 2048
_fwhm = 0.002

while i<n:
    print "i:", i
    zMin = i*dz
    zMax = (i+1)*dz
    zMinCib = (i+1)*dz
    zMaxCib = (i+2)*dz

    primary_file = "/scratch2/r/rbond/phamloui/lenspix_files/cib_v2_unlensed/cib_fullsky_ns2048_zmin%.2f_zmax%.2f_nu217_ns2048_tot.fits" % (zMinCib, zMaxCib)
    primary_smoothed_file = appendToFilename(primary_file, "fwhm_%s" % _fwhm)


    if not os.path.exists(primary_smoothed_file):
        primary_map = hp.read_map(primary_file)
        hp.mollzoom(primary_map, xsize=2048, title="Unlensed %.1f<z<%.1f" % (zMinCib, zMaxCib))
        print "smoothing map..."
        primary_map_smoothed = hp.smoothing(primary_map, fwhm=_fwhm)
        hp.mollzoom(primary_map_smoothed, xsize=2048, title="Unlensed smoothed %.1f<z<%.1f" % (zMinCib, zMaxCib))
        plt.show()
        print "writing smoothed primary to file..."
        hp.write_map(primary_smoothed_file, primary_map_smoothed)
    else:
        print "smoothed map already exists - skipping"

    i += 1
