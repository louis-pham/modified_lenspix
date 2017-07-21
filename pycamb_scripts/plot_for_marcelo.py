import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
from scipy import interpolate
import camb
import correlations #camb is not importing correlations properly, so do it manually  
from savitzky_golay import *

lens_lmax = 4000

primary_cl_file = "/scratch2/r/rbond/phamloui/lenspix_files/FromNERSC/cls/ffp10_scalCls.dat"

#interp_factor=3, marcelo+camb  :SIMULATION:
primary_file = "/scratch2/r/rbond/phamloui/lenspix_files/FromNERSC/ffp10_unlensed_scl_cmb_000_alm.fits"
lensed_file = "/scratch2/r/rbond/phamloui/lenspix_files/output/jun10_nonlinear3_interp_factor3_julian_cmb.fits"

#peakpatch+camb :ANALYTICAL:
lensed_cl_file = "/scratch2/r/rbond/phamloui/lenspix_files/jun26_peakpatch_and_camb_convolved_lensed_%s_for_marcelo.dat"
kappa_map_file = "/scratch2/r/rbond/phamloui/lenspix_files/kappa_maps/jun1_nonlinear3_kappa_for_julian_cmb.fits"

#CAMB lensed :THEORETICAL:
theoretical_lensed_cl_file = "/scratch2/r/rbond/phamloui/lenspix_files/FromNERSC/cls/ffp10_lensedCls.dat"

print "Loading primary alm and cl..."
primary_alm = hp.read_alm(primary_file)
unlensed_TT_ell, unlensed_CL_TT = np.loadtxt(primary_cl_file, usecols=(0,1), unpack=True) # l(l+1)/2 factor included 
unlensed_TT_ell = np.insert(unlensed_TT_ell, 0, 1); unlensed_TT_ell = np.insert(unlensed_TT_ell, 0, 0)
unlensed_CL_TT = np.insert(unlensed_CL_TT, 0, 0); unlensed_CL_TT = np.insert(unlensed_CL_TT, 0, 0)
unlensed_CL_TT = unlensed_CL_TT / 1e12

print "Loading theoretical lensed cl..."
lensed_theory_TT_ell, lensed_theory_CL_TT = np.loadtxt(theoretical_lensed_cl_file, usecols=(0,1), unpack=True)
lensed_theory_TT_ell = np.insert(lensed_theory_TT_ell, 0, 1); lensed_theory_TT_ell = np.insert(lensed_theory_TT_ell, 0, 0)
lensed_theory_CL_TT = np.insert(lensed_theory_CL_TT, 0, 0); lensed_theory_CL_TT = np.insert(lensed_theory_CL_TT, 0, 0)
lensed_theory_CL_TT= lensed_theory_CL_TT / 1e12

lmax = hp.Alm.getlmax(len(primary_alm)) if hp.Alm.getlmax(len(primary_alm)) < len(unlensed_TT_ell)-1 else len(unlensed_TT_ell)-1 #change this later
print lmax
unlensed_TT_ell = unlensed_TT_ell[0:lmax+1]
unlensed_CL_TT = unlensed_CL_TT[0:lmax+1]
lensed_theory_TT_ell = lensed_theory_TT_ell[0:lmax+1]
lensed_theory_CL_TT = lensed_theory_CL_TT[0:lmax+1]

# load peakpatch+camb
print "Loading kappa map..."
kappa_map = hp.read_map(kappa_map_file)
kappa_cl = hp.anafast(kappa_map, lmax=lens_lmax)
kappa_ell = np.arange(0, lmax+1)
kappa_cl = savitzky_golay(kappa_cl,75,3) #smooth cl for theory approximation                                                      
kappa_cl = np.pad(kappa_cl, (0, lmax-lens_lmax), 'constant', constant_values=(0,0))
analytic_KK_ell = kappa_ell #too lazy to refactor, just replace variables here                          
analytic_CL_KK = kappa_cl

print "Loading lensed map (sim)..."
lensed_T_map, lensed_Q_map, lensed_U_map = hp.read_map(lensed_file, field=(0,1,2))
sim_CLs_lensed = hp.anafast([lensed_T_map, lensed_Q_map, lensed_U_map], lmax=lmax)
sim_CL_TT_lensed, sim_CL_EE_lensed, sim_CL_BB_lensed = sim_CLs_lensed[0], sim_CLs_lensed[1], sim_CLs_lensed[2]
nside = hp.get_nside(lensed_T_map)
sim_TT_ell = np.arange(0, lmax+1)

print "converting kappa to phi..." #for camb convolution   
analytic_CL_PP = analytic_CL_KK / (analytic_KK_ell*(analytic_KK_ell+1)/2.0)**2
analytic_PP_ell = np.arange(0, lmax+1)
print analytic_CL_PP.shape

print "include ell factors for simulation power spectra..."
#unlensed_CL_TT = unlensed_TT_ell * (unlensed_TT_ell+1) * unlensed_CL_TT / (2*np.pi)        
analytic_CL_PP = (analytic_PP_ell * (analytic_PP_ell + 1))**2 * analytic_CL_PP / (2*np.pi)
sim_CL_TT_lensed = sim_TT_ell * (sim_TT_ell + 1) * sim_CL_TT_lensed / (2*np.pi)# * (sim_TT_ell*(sim_TT_ell+1)/4.0)
print "smoothing power spectrum from map..."
sim_CL_TT_lensed = savitzky_golay(sim_CL_TT_lensed, 75, 3)
TEB_cls = np.array([[unlensed_CL_TT[_l], 0, 0, 0] for _l in unlensed_TT_ell])
print "starting lensing convolution..."                                                                            
lensed_analytic_TEB = correlations.lensed_cls(TEB_cls, analytic_CL_PP)    
lensed_analytic_CL_TT = np.array([_cl[0] for _cl in lensed_analytic_TEB])        
print "writing camb cl..."                                                                                                      
hp.write_cl(lensed_cl_file % ('TT'), lensed_analytic_CL_TT)                                                                 

# print "reading lensed analytic cl..."
# lensed_analytic_CL_TT = hp.read_cl(lensed_cl_file % ('TT'))

lensed_analytic_TT_ell = np.arange(len(lensed_analytic_CL_TT))

#fractional differences   
print "calculating differences..."
lenspix_vs_primary_TT = (sim_CL_TT_lensed - unlensed_CL_TT) / unlensed_CL_TT
analytic_vs_primary_TT = (lensed_analytic_CL_TT - unlensed_CL_TT) / unlensed_CL_TT
theory_vs_primary_TT = (lensed_theory_CL_TT - unlensed_CL_TT) / unlensed_CL_TT

lenspix_vs_analytic_TT = (sim_CL_TT_lensed - lensed_analytic_CL_TT)/ lensed_analytic_CL_TT
lenspix_vs_theory_TT = (sim_CL_TT_lensed - lensed_theory_CL_TT)/ lensed_theory_CL_TT

plt.figure()
plt.title("TT Power Spectra", fontsize=18)
plt.xlabel(r'$l$', fontsize=30)
plt.ylabel(r'$l(l+1)*C_l / 2\pi$', fontsize=30)
plt.loglog(unlensed_TT_ell[1:], unlensed_CL_TT[1:], '#484349', label="Unlensed")
plt.loglog(sim_TT_ell[1:], sim_CL_TT_lensed[1:], '#B9314F', label="Simulation")
plt.loglog(lensed_analytic_TT_ell[1:], lensed_analytic_CL_TT[1:], '#19647E', label="Analytic")
plt.loglog(lensed_theory_TT_ell[1:], lensed_theory_CL_TT[1:], '#00CC66', label="Theoretical")
legend2 = plt.legend(loc="lower left", shadow=True)
frame2 = legend2.get_frame()
frame2.set_facecolor('0.90')

plt.figure()
plt.title("TT Fractional Differences vs. Unlensed", fontsize=14)
plt.xlabel(r"$l$", fontsize=30)
plt.ylabel(r'$\Delta C_l / C_l$', fontsize=30)
plt.plot(unlensed_TT_ell[1:], lenspix_vs_primary_TT[1:], '#B9314F', label="Simulation vs Unlensed")
plt.plot(unlensed_TT_ell[1:], analytic_vs_primary_TT[1:], '#19647E', label="Analytic vs Unlensed")
plt.plot(unlensed_TT_ell[1:], theory_vs_primary_TT[1:], '#00CC66', label="Theoretical vs Unlensed")
legend3 = plt.legend(loc="lower left", shadow=True)
frame3 = legend3.get_frame()
frame3.set_facecolor('0.90')

plt.figure()
plt.title("Analytic vs. Simulation - TT Fractional Differences", fontsize=14)
plt.xlabel(r"$l$", fontsize=30)
plt.ylabel(r'$\Delta C_l / C_l$', fontsize=30)
plt.ylim([-0.1, 0.1])
plt.plot(unlensed_TT_ell[1:], lenspix_vs_analytic_TT[1:], '#008BF8', label="Analytic vs Simulation")
plt.plot(unlensed_TT_ell[1:], lenspix_vs_theory_TT[1:], '#F5B700', label="Theory vs Simulation")
plt.grid()

legend4 = plt.legend(loc="lower left", shadow=True)
frame4 = legend4.get_frame()
frame4.set_facecolor('0.90')

plt.show()
