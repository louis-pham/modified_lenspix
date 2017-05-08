import matplotlib.pyplot as plt
import numpy as np
import healpy as hp

#unlensed_file = "test_original_lmax6250_nside4096_interp1.5_method1_pol_1_unlensed_simulated.dat"
#lensed_file = "test_original_lmax6250_nside4096_interp1.5_method1_pol_1.dat"

unlensed_file = "pycamb_cib_compare_lmax5120_nside2048_interp1.5_method1_pol_1_unlensed_power.dat"
lensed_file = "pycamb_cib_compare_lmax5120_nside2048_interp1.5_method1_pol_1_lensed_power.dat"

ells, unlensed_cls = np.loadtxt(unlensed_file, usecols=(0,1), unpack=True)
lensed_cls = np.loadtxt(lensed_file, usecols=(1,), unpack=True)

#theory_ells, theoretical_cls = np.loadtxt("sample_lenspotentialCls.dat", usecols=(0,1), unpack=True)
#theory_ells2, theoretical_lensed_cls = np.loadtxt("sample_lensedCls.dat", usecols=(0,1), unpack=True)

#print ells
#print unlensed_cls
#print lensed_cls

#plt.figure()
plt.xlabel(r'$l$')
plt.ylabel(r'$l(l+1)*C_l / 2\pi$')
plt.title("CIB Power Spectra")
#plt.loglog(theory_ells, theoretical_cls, '--g', label="Theory expectation")
#plt.loglog(theory_ells2, theoretical_lensed_cls, '--m', label="Theory exp. lensed")
plt.loglog(ells, unlensed_cls, 'r', label="Unlensed", linewidth=1.0)
plt.loglog(ells, lensed_cls, 'b', label="Lensed")
legend = plt.legend(loc='lower right', shadow=True)
frame = legend.get_frame()
frame.set_facecolor('0.90')

#print len(theoretical_cls[:2750])
#print len(theoretical_lensed_cls)

deltas_abs = lensed_cls - unlensed_cls
deltas = deltas_abs/unlensed_cls
#deltas_theory_abs = theoretical_lensed_cls[:2750] - theoretical_cls[:2750]
#deltas_theory = deltas_theory_abs / theoretical_cls[:2750]
#print deltas_abs[:200]
#print deltas_abs[-200:]
#print deltas[-200:]

plt.figure()
plt.title("Relative Differences")
plt.xlabel(r'$l$')
plt.ylabel(r'$\Delta C_l / C_l$')
plt.plot(ells, deltas, label="Lenspix delta")
#plt.plot(theory_ells2, deltas_theory, '--r', linewidth=3.0, label="Theory delta")
legend2 = plt.legend(loc="lower right", shadow=True)
frame2 = legend2.get_frame()
frame2.set_facecolor('0.90')

#deltas_pow = np.absolute(deltas) ** 2
#deltas_theory_pow = np.absolute(deltas_theory) ** 2
#plt.figure()
#plt.xlabel('l')
#plt.ylabel(r'$|\Delta C_l / C_l|^2$')
#plt.plot(ells, deltas_pow, label="Lenspix delta")
#plt.plot(theory_ells2, deltas_theory_pow, '--r', label="Theory delta")               
#legend3 = plt.legend(loc="lower right", shadow=True)                                   
#frame3 = legend3.get_frame()                                                            
#frame3.set_facecolor('0.90') 

#plt.figure() 
#plt.xlabel('l')                                                                       
#plt.ylabel(r'$\Delta C_l$')                                  
#plt.plot(ells, deltas_abs, label="Lenspix delta")
#plt.plot(ells, np.zeros(len(ells)))
plt.show()
