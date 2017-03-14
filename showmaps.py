import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import healpy as hp

nside = 2048

kappa_map_file = "kappa_maps/8Gpc_n2048_nb23_nt18_kap_comb.fits"
phi_map_file = "kappa_maps/test_cib_lens_n2048_lmax5120_phi_map.fits"
unlensed_file = "kappa_maps/cib_fullsky_ns2048_zmin0.0_zmax1.245_nu217_13579_normalized_alm.fits" #alm
lensed_file = "test_cib_lens_lmax5120_nside2048_interp1.5_method1_pol_1_lensed.fits" #map

kappa_map = hp.read_map(kappa_map_file)
phi_map = hp.read_map(phi_map_file)

lensed_map, lensed_header = hp.read_map(lensed_file, h=True)
unlensed_alm = hp.read_alm(unlensed_file)
unlensed_map = hp.alm2map(unlensed_alm, nside)
#unlensed_map = hp.read_map(unlensed_file)

diff_map = lensed_map - unlensed_map

print diff_map

temp_bound = np.maximum(np.absolute(lensed_map).max(), np.absolute(unlensed_map).max())
diff_bound = max(np.absolute(diff_map))

hot_cmap = cm.get_cmap("hot")
hot_cmap.set_under("w")
seismic_cmap = cm.get_cmap("seismic")
seismic_cmap.set_under("w")

#hp.cartview(kappa_map, title="Kappa", xsize=2048, cmap=hot_cmap, lonra=[-10,10], latra=[-10,10])
#hp.cartview(phi_map, title="Phi", xsize=2048, cmap=hot_cmap, lonra=[-10,10], latra=[-10,10])
#hp.mollview(kappa_map, title="Kappa", xsize=2048, cmap=hot_cmap)
#hp.mollview(phi_map, title="Phi", xsize=2048, cmap=hot_cmap)
#hp.mollview(unlensed_map, title="Unlensed", xsize=2048, cmap=hot_cmap, min=-temp_bound, max=temp_bound, unit="muK")
#hp.mollview(lensed_map, title="Lensed", xsize=2048, cmap=hot_cmap, min=-temp_bound, max=temp_bound, unit="muK")
#hp.mollview(diff_map, title="Diff map", xsize=2048, cmap=seismic_cmap, min=-diff_bound, max=diff_bound, unit="muK") 
#hp.cartview(unlensed_map, title="Unlensed", xsize=2048, cmap=hot_cmap, min=-temp_bound, max=temp_bound, unit="muK", lonra=[-10,10], latra=[-10,10])
#hp.cartview(lensed_map, title="Lensed", xsize=2048, cmap=hot_cmap, min=-temp_bound, max=temp_bound, unit="muK", lonra=[-10,10], latra=[-10,10])
#hp.cartview(diff_map, title="Diff map", xsize=2048, cmap=seismic_cmap, min=-diff_bound, max=diff_bound, lonra=[-10,10], latra=[-10,10], unit="muK")

#manual temp bounds:
hp.mollview(unlensed_map, title="Unlensed manual", min=-14.8063, max=127.516, cmap=hot_cmap, xsize=2048, unit="muK")  
hp.mollview(lensed_map, title="Lensed manual", min=-14.8063, max=127.516, cmap=hot_cmap, xsize=2048, unit="muK")
hp.cartview(unlensed_map, title="Unlensed manual 10deg", min=-4.56, max=38.5, xsize=2048, cmap=hot_cmap, lonra=[-10,10], latra=[-10,10])
hp.cartview(lensed_map, title="Lensed manual 10deg", min=-4.56, max=38.5, xsize=2048, cmap=hot_cmap, lonra=[-10,10], latra=[-10,10])

plt.show()
