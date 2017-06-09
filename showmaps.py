import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import healpy as hp

nside = 2048

kappa_map_file = "/scratch2/r/rbond/phamloui/lenspix_files/lensing/kappa_z4pt6_z0.0_z0.2_nside2048_hp.fits"
phi_alm_file = "/scratch2/r/rbond/phamloui/lenspix_files/output/jun7_cib_test_first_slice_phi.fits"
unlensed_file = "/scratch2/r/rbond/phamloui/lenspix_files/cib/cib_fullsky_ns2048_zmin0.00_zmax0.20_nu217_tot.fits" #alm
lensed_file = "/scratch2/r/rbond/phamloui/lenspix_files/output/jun7_cib_test_first_slice_lensed.fits" #map
#kappa_map_file2 = "/scratch2/r/rbond/phamloui/lenspix_files/output/rotated_maps/8Gpc_n2048_nb23_nt18_kap_comb.fits"
kappa_map = hp.read_map(kappa_map_file)
#kappa_map2 = hp.read_map(kappa_map_file2)
phi_alm = hp.read_alm(phi_alm_file)
phi_map = hp.alm2map(phi_alm, nside)
#phi_map = hp.read_map("/scratch2/r/rbond/phamloui/lenspix_files/output/vanilla_lmax3000_nside2048_interp1.5_method1_pol_1_phi.fits")

lensed_map, lensed_header = hp.read_map(lensed_file, h=True)
lensed_map = np.nan_to_num(lensed_map)
#unlensed_alm = hp.read_alm(unlensed_file)
#unlensed_map = hp.alm2map(unlensed_alm.astype('D'), nside)
unlensed_map = hp.read_map(unlensed_file)
unlensed_alm = hp.map2alm(unlensed_map)

diff_map = lensed_map - unlensed_map

#print diff_map

temp_max = np.maximum(lensed_map.max(), unlensed_map.max())
temp_min = np.minimum(lensed_map.min(), unlensed_map.min())
diff_bound = np.absolute(diff_map.max())

hot_cmap = cm.get_cmap("hot")
hot_cmap.set_under("w")

seismic_cmap = cm.get_cmap("seismic")
seismic_cmap.set_under("w")

#hp.cartview(kappa_map, title="Kappa", xsize=2048, cmap=hot_cmap, lonra=[-10,10], latra=[-10,10])
#hp.cartview(phi_map, title="Phi", xsize=2048, cmap=hot_cmap, lonra=[-10,10], latra=[-10,10])
hp.mollview(kappa_map, title="Kappa", xsize=2048, cmap=hot_cmap)
#hp.mollview(kappa_map2, title="Kappa 2", xsize=2048, cmap=hot_cmap)
hp.mollview(phi_map, title="Phi", xsize=2048, cmap=hot_cmap)
hp.mollview(unlensed_map, title="Unlensed", xsize=2048, cmap=hot_cmap, min=temp_min, max=temp_max, unit="muK")
hp.mollview(lensed_map, title="Lensed", xsize=2048, cmap=hot_cmap, min=temp_min, max=temp_max, unit="muK")
hp.mollview(diff_map, title="Diff map", xsize=2048, cmap=seismic_cmap, min=-diff_bound, max=diff_bound, unit="muK") 
#hp.cartview(unlensed_map, title="Unlensed", xsize=2048, cmap=hot_cmap, min=-temp_bound, max=temp_bound, unit="muK", lonra=[-10,10], latra=[-10,10])
#hp.cartview(lensed_map, title="Lensed", xsize=2048, cmap=hot_cmap, min=-temp_bound, max=temp_bound, unit="muK", lonra=[-10,10], latra=[-10,10])
#hp.cartview(diff_map, title="Diff map", xsize=2048, cmap=seismic_cmap, min=-diff_bound, max=diff_bound, lonra=[-10,10], latra=[-10,10], unit="muK")

#manual temp bounds:
#hp.mollview(unlensed_map, title="Unlensed manual", min=-14.8063, max=127.516, cmap=hot_cmap, xsize=2048, unit="muK")  
#hp.mollview(lensed_map, title="Lensed manual", min=-14.8063, max=127.516, cmap=hot_cmap, xsize=2048, unit="muK")
#hp.cartview(unlensed_map, title="Unlensed manual 10deg", min=-4.56, max=38.5, xsize=2048, cmap=hot_cmap, lonra=[-10,10], latra=[-10,10])
#hp.cartview(lensed_map, title="Lensed manual 10deg", min=-4.56, max=38.5, xsize=2048, cmap=hot_cmap, lonra=[-10,10], latra=[-10,10])

plt.show()
