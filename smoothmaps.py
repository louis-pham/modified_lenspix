import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from matplotlib import cm

def gaussFunc(ell, ell_o):
    return np.exp(-(ell ** 2 / ell_o ** 2))

def appendToFilename(filename, newStr):
    return "{0}_{2}.{1}".format(*filename.rsplit('.', 1) + [newStr])

#primary_file = "/scratch2/r/rbond/phamloui/lenspix_files/FromNERSC/ffp10_unlensed_scl_cmb_000_alm.fits"
kappa_map_file = "/scratch2/r/rbond/phamloui/lenspix_files/kappa_maps/8Gpc_n2048_nb23_nt18_kap_comb.fits"
primary_file = "/scratch2/r/rbond/phamloui/lenspix_files/primary_maps/cib_fullsky_ns2048_zmin0.0_zmax1.245_nu217_13579_normalized_alm.fits"
#kappa_map_file = "/scratch2/r/rbond/phamloui/lenspix_files/kappa_maps/8Gpc_n2048_nb23_nt18_kap_comb.fits"


l_o = 4000.

primary_smooth_file = appendToFilename(primary_file, "gauss_smoothed_ell_o_%.0f" % (l_o))
kappa_smooth_file = appendToFilename(kappa_map_file, "gauss_smoothed_ell_o_%.0f" % (l_o))

print primary_smooth_file
print kappa_smooth_file

primary_alm = hp.read_alm(primary_file)
primary_map = hp.alm2map(primary_alm, 2048) 
lmax = hp.Alm.getlmax(len(primary_alm))
l,m = hp.Alm.getlm(lmax)
#primary_cl = hp.alm2cl(primary_alm)

kappa_map = hp.read_map(kappa_map_file)
kappa_alm = hp.map2alm(kappa_map, lmax=lmax)
#kappa_cl = hp.alm2cl(kappa_alm)

smooth_array = gaussFunc(l, l_o)
#plt.figure()
#plt.plot(l, smooth_array)
#plt.show()

primary_smooth = primary_alm * smooth_array
kappa_smooth = kappa_alm * smooth_array

#primary_smooth_cl = hp.alm2cl(primary_smooth)
#kappa_smooth_cl = hp.alm2cl(kappa_smooth)

primary_smooth_map = hp.alm2map(primary_smooth, 2048)
kappa_smooth_map = hp.alm2map(kappa_smooth, 2048)

hot_cmap = cm.get_cmap("hot")
hot_cmap.set_under("w")

print primary_map.max(), primary_smooth_map.max()
print primary_map.min(), primary_smooth_map.min()
print kappa_map.max(), kappa_smooth_map.max()
print kappa_map.min(), kappa_smooth_map.min()

temp_max = np.maximum(primary_map.max(), primary_smooth_map.max())
temp_min = np.minimum(primary_map.min(), primary_smooth_map.min())

kap_max = np.maximum(kappa_map.max(), kappa_smooth_map.max())
kap_min = np.minimum(kappa_map.min(), kappa_smooth_map.min())

#hp.mollview(kappa_map, title="Kappa", xsize=2048, cmap=hot_cmap, min=kap_min, max=kap_max)
#hp.mollview(kappa_smooth_map, title="Kappa smooth", xsize=2048, cmap=hot_cmap, min=kap_min, max=kap_max)
hp.mollview(primary_map, title="Primary", xsize=2048, cmap=hot_cmap, min=temp_min, max=temp_max)
hp.mollview(primary_smooth_map, title="Primary smooth", xsize=2048, cmap=hot_cmap, min=temp_min, max=temp_max)

#ell = np.arange(len(primary_cl))

#plt.figure()
#plt.title("Primary cls")
#plt.loglog(ell, primary_cl, 'r')
#plt.loglog(ell, primary_smooth_cl, 'b')
#plt.figure()
#plt.title("Kappa cls")
#plt.loglog(ell, kappa_cl, 'r')
#plt.loglog(ell, kappa_smooth_cl, 'b')

plt.show()

print "writing smoothed primary to file..."
hp.write_alm(primary_smooth_file, primary_smooth)
print "writing smoothed kappa to file..."
hp.write_map(kappa_smooth_file, kappa_smooth_map)
