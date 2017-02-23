import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import healpy as hp

nside = 4096

kappa_map_file = "kappa_maps/downsized_kappa_n4096.fits"
phi_map_file = "kappa_maps/z_4.6_n4096_lmax5120_phi_map.fits"
unlensed_file = "FromNERSC/ffp10_unlensed_scl_cmb_000_alm.fits" #alm
lensed_file = "z4.6_new_primary_lmax5120_nside4096_interp1.5_method1_pol_1_lensed.fits" #map

kappa_map = hp.read_map(kappa_map_file)
phi_map = hp.read_map(phi_map_file)

lensed_map, lensed_header = hp.read_map(lensed_file, h=True)
unlensed_alm = hp.read_alm(unlensed_file)
unlensed_map = hp.alm2map(unlensed_alm, nside)
#unlensed_map = hp.read_map(unlensed_file)

#diff_map = lensed_map - unlensed_map

temp_bound = np.maximum(np.absolute(lensed_map).max(), np.absolute(unlensed_map).max())
#diff_bound = np.ceil(max(np.absolute(diff_map)))

hot_cmap = cm.get_cmap("hot")
hot_cmap.set_under("w")
seismic_cmap = cm.get_cmap("seismic")
seismic_cmap.set_under("w")

#hp.cartview(kappa_map, title="Kappa", xsize=2048, cmap=hot_cmap, lonra=[-10,10], latra=[-10,10])
#hp.cartview(phi_map, title="Phi", xsize=2048, cmap=hot_cmap, lonra=[-10,10], latra=[-10,10])
hp.mollview(kappa_map, title="Kappa", xsize=2048, cmap=hot_cmap)
hp.mollview(phi_map, title="Phi", xsize=2048, cmap=hot_cmap)
hp.mollview(unlensed_map, title="Unlensed", xsize=2048, cmap=hot_cmap, min=-temp_bound, max=temp_bound, unit="muK")
hp.mollview(lensed_map, title="Lensed", xsize=2048, cmap=hot_cmap, min=-temp_bound, max=temp_bound, unit="muK")
#hp.cartview(unlensed_map, title="Unlensed", xsize=2048, cmap=hot_cmap, min=-temp_bound, max=temp_bound, unit="muK", lonra=[-10,10], latra=[-10,10])
#hp.cartview(lensed_map, title="Lensed", xsize=2048, cmap=hot_cmap, min=-temp_bound, max=temp_bound, unit="muK", lonra=[-10,10], latra=[-10,10])
#hp.cartview(diff_map, title="Diff map", xsize=2048, cmap=seismic_cmap, min=-diff_bound, max=diff_bound, lonra=[-10,10], latra=[-10,10], unit="muK")
plt.show()
