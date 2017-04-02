import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import camb
import correlations #camb is not importing correlations properly, so do it manually
from savitzky_golay import *

primary_file = "/scratch2/r/rbond/phamloui/lenspix/kappa_maps/cib_fullsky_ns2048_zmin0.0_zmax1.245_nu217_13579_normalized_alm.fits"
#kappa_halo_file = "/scratch2/r/rbond/phamloui/lenspix/kappa_maps/8Gpc_n4096_nb23_nt18_kap_halo.fits"
#kappa_field_file =  "/scratch2/r/rbond/phamloui/lenspix/kappa_maps/8Gpc_n4096_nb23_nt18_kap_field.fits"
phi_alm_file = "/scratch2/r/rbond/phamloui/lenspix/kappa_maps/pycamb_cib_compare_n2048_phi_alm.fits"
lensed_file = "/scratch2/r/rbond/phamloui/lenspix/pycamb_cib_compare_lmax5120_nside2048_interp1.5_method1_pol_1_lensed.fits"

print "Loading primary..."
primary_alm = hp.read_alm(primary_file)
primary_cl = hp.alm2cl(primary_alm)
lmax = hp.Alm.getlmax(len(primary_alm))

#kappa_halo_map = hp.read_map(kappa_halo_file)
#kappa_field_map = hp.read_map(kappa_field_file)

print "Loading phi..."
phi_alm = hp.read_alm(phi_alm_file)
phi_cl = hp.alm2cl(phi_alm)

print "Loading lensed map..."
lensed_map = hp.read_map(lensed_file)
lensed_cl = hp.anafast(lensed_map, lmax=lmax)
nside = hp.get_nside(lensed_map)

ell = np.arange(len(primary_cl))

#include factors
primary_cl = ell * (ell+1) * primary_cl / (2*np.pi)
phi_cl = (ell * (ell + 1))**2 * phi_cl / (2*np.pi)
lensed_cl = ell * (ell + 1) * lensed_cl / (2*np.pi)

#fill EE,BB,TE with zeroes
TEB_cls = np.array([[primary_cl[_l], 0, 0, 0] for _l in ell])

lensed_theory_TEB = correlations.lensed_cls(TEB_cls, phi_cl)
lensed_theory_T = np.array([_l[0] for _l in lensed_theory_TEB])
hp.write_cl("lensed_theoryCls.dat", lensed_theory_T)

#smooth the spectra
primary_cl = savitzky_golay(primary_cl, 75, 3)
lensed_cl = savitzky_golay(lensed_cl, 75, 3)
lensed_theory_T = savitzky_golay(lensed_theory_T, 75, 3)

plt.figure()
plt.title("Power Spectra")
plt.xlabel(r'$l$')
plt.ylabel(r'$l(l+1)*C_l / 2\pi$')

plt.loglog(ell, primary_cl, 'r', label="Unlensed")
plt.loglog(ell, lensed_cl, 'b', label="Lensed")
plt.loglog(ell, lensed_theory_T, '--g', label="Lensed (CAMB)")

plt.grid()
legend2 = plt.legend(loc="lower right", shadow=True)
frame2 = legend2.get_frame()
frame2.set_facecolor('0.90')

plt.figure()
plt.title("Fractional Difference - Unlensed vs Lensed")
plt.xlabel(r"$l$")
plt.ylabel(r'$\Delta C_l / C_l$')
plt.plot((lensed_cl - primary_cl) / primary_cl)
plt.grid()

plt.figure()
plt.title("Fractional Difference - Lensed (pyCAMB) vs Lensed (Lenspix)")
plt.xlabel(r"$l$")
plt.ylabel(r'$\Delta C_l / C_l$')
plt.plot((lensed_cl - lensed_theory_T) / lensed_theory_T)
plt.grid()

plt.show()
