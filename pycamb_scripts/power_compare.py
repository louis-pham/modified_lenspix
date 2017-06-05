import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import camb
import correlations #camb is not importing correlations properly, so do it manually
from savitzky_golay import *

graphTitleRoot = ""

#original cib comparison
#primary_file = "/scratch2/r/rbond/phamloui/lenspix_files/primary_maps/cib_fullsky_ns2048_zmin0.0_zmax1.245_nu217_13579_normalized_alm.fits"
#phi_alm_file = "/scratch2/r/rbond/phamloui/lenspix_files/output/may17_cib_test_n2048_phi_alm.fits"
#lensed_file = "/scratch2/r/rbond/phamloui/lenspix_files/output/may17_cib_test_fullsky_ns2048_lensed.fits"
#theory_cl_file = "/scratch2/r/rbond/phamloui/lenspix_files/may17_lensed_theoryCls_cib_discrepancy.dat"

#rotate test
#primary_file = "/scratch2/r/rbond/phamloui/lenspix_files/primary_maps/cib_fullsky_ns2048_zmin0.0_zmax1.245_nu217_13579_normalized_alm.fits"
#phi_alm_file = "/scratch2/r/rbond/phamloui/lenspix_files/output/may15_test_rotate_n2048_phi_alm.fits"
#lensed_file = "/scratch2/r/rbond/phamloui/lenspix_files/output/may15_test_rotate_cib_fullsky_ns2048_lensed.fits"
#theory_cl_file = "/scratch2/r/rbond/phamloui/lenspix_files/may15_lensed_theoryCls_rotate_test.dat"

#cmb test
#primary_file = "/scratch2/r/rbond/phamloui/lenspix_files/FromNERSC/ffp10_unlensed_scl_cmb_000_alm.fits"
#phi_alm_file = "/scratch2/r/rbond/phamloui/lenspix_files/output/8Gpc_n2048_nb23_nt18_phi_comb.fits"
#lensed_file = "/scratch2/r/rbond/phamloui/lenspix_files/output/may15_NERSC_scl_cmb_lensed.fits"
#theory_cl_file = "/scratch2/r/rbond/phamloui/lenspix_files/may15_lensed_theoryCls.dat"

#vanilla test
#primary_file = "/scratch2/r/rbond/phamloui/lenspix_files/output/vanilla_lmax3000_nside2048_interp1.5_method1_pol_1_unlensed_simulated_alm.fits"
#phi_alm_file = "/scratch2/r/rbond/phamloui/lenspix_files/output/vanilla_lmax3000_nside2048_interp1.5_method1_pol_1_phi.fits"
#lensed_file = "/scratch2/r/rbond/phamloui/lenspix_files/output/vanilla_lmax3000_nside2048_interp1.5_method1_pol_1_lensed_map.fits"
#theory_cl_file = "/scratch2/r/rbond/phamloui/lenspix_files/may18_lensed_theoryCls_vanilla.dat"

#smooth test cmb
#primary_file = "/scratch2/r/rbond/phamloui/lenspix_files/FromNERSC/ffp10_unlensed_scl_cmb_000_alm_gauss_smoothed.fits"
#phi_alm_file = "/scratch2/r/rbond/phamloui/lenspix_files/output/may19_test_smooth_n2048_phi_alm.fits"
#lensed_file = "/scratch2/r/rbond/phamloui/lenspix_files/output/may23_test_smooth_cmb_lensed.fits"
#theory_cl_file = "/scratch2/r/rbond/phamloui/lenspix_files/may23_lensed_theoryCls_smooth_cmb_test.dat"

#smooth test cib
#primary_file = "/scratch2/r/rbond/phamloui/lenspix_files/primary_maps/cib_fullsky_ns2048_zmin0.0_zmax1.245_nu217_13579_normalized_alm_gauss_smoothed.fits"
#phi_alm_file = "/scratch2/r/rbond/phamloui/lenspix_files/output/may19_test_smooth_n2048_phi_alm.fits"
#lensed_file = "/scratch2/r/rbond/phamloui/lenspix_files/output/may19_test_smooth_cib_fullsky_ns2048_lensed.fits"
#theory_cl_file = "/scratch2/r/rbond/phamloui/lenspix_files/may19_lensed_theoryCls_smooth_test.dat"

#lpt1
#primary_file = "/scratch2/r/rbond/phamloui/lenspix_files/FromNERSC/ffp10_unlensed_scl_cmb_000_alm.fits"
#phi_alm_file = "/scratch2/r/rbond/phamloui/lenspix_files/output/8Gpc_n4096_nb18_nt16_field_lpt1_fwh_phi.fits"
#lensed_file = "/scratch2/r/rbond/phamloui/lenspix_files/output/may25_lpt1_ffp10_lensed_scl_cmb_000_alm.fits"
#theory_cl_file = "/scratch2/r/rbond/phamloui/lenspix_files/may25_lpt1_lensed_theoryCls.dat"
#graphTitleRoot = "Julian CMB + LPT1 Run"

#lpt2
#primary_file = "/scratch2/r/rbond/phamloui/lenspix_files/FromNERSC/ffp10_unlensed_scl_cmb_000_alm.fits"
#phi_alm_file = "/scratch2/r/rbond/phamloui/lenspix_files/output/8Gpc_n4096_nb18_nt16_field_lpt2_fwh_phi.fits"
#lensed_file = "/scratch2/r/rbond/phamloui/lenspix_files/output/may25_lpt2_ffp10_lensed_scl_cmb_000_alm.fits"
#theory_cl_file = "/scratch2/r/rbond/phamloui/lenspix_files/may25_lpt2_lensed_theoryCls.dat"
#graphTitleRoot = "Julian CMB + LPT2 Run"

#smooth test cib #2
#primary_file = "/scratch2/r/rbond/phamloui/lenspix_files/primary_maps/cib_fullsky_ns2048_zmin0.0_zmax1.245_nu217_13579_normalized_alm_gauss_smoothed_ell_o_3000.fits"
#phi_alm_file = "/scratch2/r/rbond/phamloui/lenspix_files/output/may26_test_smooth_n2048_phi_alm_ell_o_3000.fits"
#lensed_file = "/scratch2/r/rbond/phamloui/lenspix_files/output/may26_test_smooth_ell_o_3000_cib_fullsky_ns2048_lensed.fits"
#theory_cl_file = "/scratch2/r/rbond/phamloui/lenspix_files/may26_lensed_theoryCls_smooth_ell_o_3000.dat"
#graphTitleRoot = "Smoothed CIB 2"

#smooth test cib #3
#primary_file = "/scratch2/r/rbond/phamloui/lenspix_files/primary_maps/cib_fullsky_ns2048_zmin0.0_zmax1.245_nu217_13579_normalized_alm_gauss_smoothed_ell_o_5000.fits"
#phi_alm_file = "/scratch2/r/rbond/phamloui/lenspix_files/output/may26_test_smooth_n2048_phi_alm_ell_o_5000.fits"
#lensed_file = "/scratch2/r/rbond/phamloui/lenspix_files/output/may26_test_smooth_ell_o_5000_cib_fullsky_ns2048_lensed.fits"
#theory_cl_file = "/scratch2/r/rbond/phamloui/lenspix_files/may26_lensed_theoryCls_smooth_ell_o_5000.dat"
#graphTitleRoot = "Smoothed CIB 3"

#only smooth kappa 3000:
#primary_file = "/scratch2/r/rbond/phamloui/lenspix_files/FromNERSC/ffp10_unlensed_scl_cmb_000_alm.fits"
#phi_alm_file = "/scratch2/r/rbond/phamloui/lenspix_files/output/may29_only_kappa_smooth_ell_o_3000_n2048_phi_alm.fits"
#lensed_file = "/scratch2/r/rbond/phamloui/lenspix_files/output/may29_only_kappa_smooth_ell_o_3000_cmb_lensed.fits"
#theory_cl_file = "/scratch2/r/rbond/phamloui/lenspix_files/may29_lensed_theoryCls_only_kappa_smooth_ell_o_3000.dat"

#only smooth kappa 5000:
#primary_file = "/scratch2/r/rbond/phamloui/lenspix_files/FromNERSC/ffp10_unlensed_scl_cmb_000_alm.fits"
#phi_alm_file = "/scratch2/r/rbond/phamloui/lenspix_files/output/may29_only_kappa_smooth_ell_o_5000_n2048_phi_alm.fits"
#lensed_file = "/scratch2/r/rbond/phamloui/lenspix_files/output/may29_only_kappa_smooth_ell_o_5000_cmb_lensed.fits"
#theory_cl_file = "/scratch2/r/rbond/phamloui/lenspix_files/may29_lensed_theoryCls_only_kappa_smooth_ell_o_5000.dat"

#sim vs theory:
primary_file = "/scratch2/r/rbond/phamloui/lenspix_files/FromNERSC/ffp10_alm_from_scalCls.fits"
phi_alm_file = "/scratch2/r/rbond/phamloui/lenspix_files/output/may29_n2048_phi_restricted_alm.fits"
lensed_file = "/scratch2/r/rbond/phamloui/lenspix_files/output/may29_kappa_restricted_ffp10_scalCls_lensed.fits"
kappa_cl_file = "/scratch2/r/rbond/phamloui/lenspix_files/cl_kappa_restricted.txt"
primary_theory_cl_file = "/scratch2/r/rbond/phamloui/lenspix_files/FromNERSC/cls/ffp10_scalCls.dat"

print "Loading primary..."
primary_alm = hp.read_alm(primary_file)
primary_theory_ell, primary_theory_cl = np.loadtxt(primary_theory_cl_file, usecols(0,1), unpack=True)
lmax = hp.Alm.getlmax(len(primary_alm))
#kappa_halo_map = hp.read_map(kappa_halo_file)
#kappa_field_map = hp.read_map(kappa_field_file)

print "Loading phi..."
phi_alm = hp.read_alm(phi_alm_file)
phi_cl = hp.alm2cl(phi_alm)

#vanilla lenspix has phi map
#phi_map = hp.read_map(phi_alm_file)
#phi_cl = hp.anafast(phi_map, lmax=lmax)

print "Loading lensed map..."
lensed_map = hp.read_map(lensed_file)
lensed_cl = hp.anafast(lensed_map, lmax=lmax)
nside = hp.get_nside(lensed_map)

ell = np.arange(len(primary_cl)) 
#include factors
primary_cl = ell * (ell+1) * primary_cl / (2*np.pi)
phi_cl = (ell * (ell + 1))**2 * phi_cl / (2*np.pi)
lensed_cl = ell * (ell + 1) * lensed_cl / (2*np.pi)

#set l=0,1 to 0 because reasons
primary_cl[0] = 0
primary_cl[1] = 0 
phi_cl[0] = 0
phi_cl[1] = 0

#print primary_cl
#print phi_cl
#print lensed_cl

#smooth the spectra
primary_cl = savitzky_golay(primary_cl, 75, 3)
phi_cl = savitzky_golay(phi_cl, 75, 3)
lensed_cl = savitzky_golay(lensed_cl, 75, 3)

#print "primary, phi, lensed:"
#print primary_cl
#print phi_cl
#print lensed_cl

#fill EE,BB,TE with zeroes
TEB_cls = np.array([[primary_cl[_l], 0, 0, 0] for _l in ell])

lensed_theory_TEB = correlations.lensed_cls(TEB_cls, phi_cl)
lensed_theory_T = np.array([_l[0] for _l in lensed_theory_TEB])
#print "writing camb cl..."
#hp.write_cl(theory_cl_file, lensed_theory_T)

#print "Loading camb cl"
#lensed_theory_T = hp.read_cl(theory_cl_file)

#print "lensed_theory_T:"
#print lensed_theory_T

#fractional differences
lenspix_vs_primary = (lensed_cl - primary_cl) / primary_cl
camb_vs_lenspix = (lensed_cl - lensed_theory_T) / lensed_theory_T
camb_vs_primary = (lensed_theory_T - primary_cl) / primary_cl

#print lensed_vs_primary
#print camb_vs_lenspix
#print lenspix_vs_primary

plt.figure()
plt.title(graphTitleRoot + " - Power Spectra", fontsize=18)
plt.xlabel(r'$l$', fontsize=30)
plt.ylabel(r'$l(l+1)*C_l / 2\pi$', fontsize=30)

plt.loglog(ell[1:], primary_cl[1:], 'b', label="Unlensed")
plt.loglog(ell[1:], lensed_cl[1:], 'r', label="Simulation")
plt.loglog(ell[1:], lensed_theory_T[1:], 'g', label="Theory")

plt.grid()
legend2 = plt.legend(loc="lower left", shadow=True)
frame2 = legend2.get_frame()
frame2.set_facecolor('0.90')

plt.figure()
plt.title(graphTitleRoot + " - Fractional Differences", fontsize=18)
plt.xlabel(r"$l$", fontsize=30)
plt.ylabel(r'$\Delta C_l / C_l$', fontsize=30)
plt.plot(ell[1:], lenspix_vs_primary[1:], 'r', label="Simulation vs Unlensed")
plt.plot(ell[1:], camb_vs_primary[1:], 'g', label="Theory vs Unlensed")
#plt.ylim([-0.1, 0.1])
legend3 = plt.legend(loc="lower left", shadow=True)
frame3 = legend3.get_frame()
frame3.set_facecolor('0.90')
plt.grid()
#apr3_fractional_overlay.pdf

plt.figure()
plt.title(graphTitleRoot + " - Theory vs Simulation Fractional Differences", fontsize=18)
plt.xlabel(r"$l$", fontsize=30)
plt.ylabel(r"$\Delta C_l / C_l$", fontsize=30)
plt.plot(ell[1:], camb_vs_lenspix[1:], 'm', label="Theory vs Simulation")
#plt.ylim([-0.1,0.1])
plt.grid()

plt.show()
