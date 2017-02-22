import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import healpy as hp

nside = 4096
#kappa_map_file = "kappa_maps/" + str(nside) + "_kappa_map.fits"
kappa_map_file = "kappa_maps/8Gpc_n4096_nb23_nt18_kap_13579_hp.fits"
phi_map_file = "kappa_maps/n" + str(nside) + "_lmax3000_phi_map.fits"
#ksz_file = "kappa_maps/8Gpc_n4096_nb23_nt18_ksz_13579_hp.fits"
unlensed_file = "test_primary.fits"
lensed_file = "test_custom_lmax3000_nside4096_interp1.5_method1_pol_1_lensed.fits" #map

#default lenspix maps
#default_unlensed_file = "test_original_lmax2250_nside4096_interp1.5_method1_pol_1_unlensed.fits"
#default_lensed_file = "test_original_lmax2250_nside4096_interp1.5_method1_pol_1_lensed.fits"

kappa_map = hp.read_map(kappa_map_file)
phi_map = hp.read_map(phi_map_file)
#ksz_map, ksz_h = hp.read_map(ksz_file, h=True)
#print ksz_h
#ksz_map = ksz_map * 2.72e6
#ksz_alm = hp.map2alm(ksz_map, lmax=2250)

#unlensed_alm = hp.read_alm(unlensed_file)
lensed_map, lensed_header = hp.read_map(lensed_file, h=True)
#unlensed_map = hp.alm2map(unlensed_alm.astype('D'), nside)
unlensed_map = hp.read_map(unlensed_file)

#default_unlensed_map = hp.alm2map(hp.read_alm(default_unlensed_file).astype('D'), nside)
#default_lensed_map = hp.read_map(default_lensed_file)

#diff_map = lensed_map - unlensed_map

#print ksz_map
#print (ksz_alm * 2.72e6).max()

#unlensed_and_ksz = hp.alm2map(unlensed_alm + ksz_alm, nside)

temp_bound = np.maximum(np.absolute(lensed_map).max(), np.absolute(unlensed_map).max())
#temp_bound = np.maximum(temp_bound, np.absolute(unlensed_and_ksz).max())

#diff_bound = np.ceil(max(np.absolute(diff_map)))

hot_cmap = cm.get_cmap("hot")
hot_cmap.set_under("w")
seismic_cmap = cm.get_cmap("seismic")
seismic_cmap.set_under("w")
hp.cartview(kappa_map, title="Kappa", xsize=2048, cmap=hot_cmap, lonra=[-10,10], latra=[-10,10])
hp.cartview(phi_map, title="Phi", xsize=2048, cmap=hot_cmap, lonra=[-10,10], latra=[-10,10])

#hp.mollview(phi_map, title="Phi", xsize=2048, cmap=new_cmap)
hp.cartview(unlensed_map, title="Unlensed", xsize=2048, cmap=hot_cmap, min=-temp_bound, max=temp_bound, unit="muK", lonra=[-10,10], latra=[-10,10])
hp.cartview(lensed_map, title="Lensed", xsize=2048, cmap=hot_cmap, min=-temp_bound, max=temp_bound, unit="muK", lonra=[-10,10], latra=[-10,10])
#hp.cartview(default_unlensed_map, title="DEFAULT Unlensed", xsize=2048, cmap=hot_cmap, min=-temp_bound, max=temp_bound, unit="muK", lonra=[-10,10], latra=[-10,10])
#hp.cartview(default_lensed_map, title="DEFAULT Lensed", xsize=2048, cmap=hot_cmap, min=-temp_bound, max=temp_bound, unit="muK", lonra=[-10,10], latra=[-10,10])
#hp.cartview(diff_map, title="Diff map", xsize=2048, cmap=seismic_cmap, min=-diff_bound, max=diff_bound, lonra=[-10,10], latra=[-10,10], unit="muK")
#hp.cartview(unlensed_and_ksz, title="Unlensed and kSZ", xsize=2048, cmap=hot_cmap, min=-temp_bound, max=temp_bound, unit="muK", lonra=[-10,10], latra=[-10,10])
plt.show()
