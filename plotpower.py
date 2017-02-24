import matplotlib.pyplot as plt
import numpy as np
import healpy as hp

#unlensed_file = "test_original_lmax6250_nside4096_interp1.5_method1_pol_1_unlensed_simulated.dat"
#lensed_file = "test_original_lmax6250_nside4096_interp1.5_method1_pol_1.dat"

unlensed_file = "z4.6_new_primary_lmax5120_nside4096_interp1.5_method1_pol_1_unlensed_power.dat"
lensed_file = "z4.6_new_primary_lmax5120_nside4096_interp1.5_method1_pol_1_lensed_power.dat"

ells, unlensed_cls = np.loadtxt(unlensed_file, usecols=(0,1), unpack=True)
lensed_cls = np.loadtxt(lensed_file, usecols=(1,), unpack=True)

theory_ells, theoretical_cls = np.loadtxt("sample_lenspotentialCls.dat", usecols=(0,1), unpack=True)
theory_ells2, theoretical_lensed_cls = np.loadtxt("sample_lensedCls.dat", usecols=(0,1), unpack=True)

#print ells
#print unlensed_cls
#print lensed_cls

#plt.figure()
plt.xlabel('l')
plt.ylabel('l(l+1)*C_l / 2pi')
#plt.loglog(theory_ells, theoretical_cls, '--g', label="Theory expectation")
#plt.loglog(theory_ells2, theoretical_lensed_cls, '--m', label="Theory exp. lensed")
plt.loglog(ells, unlensed_cls, ':r', label="Unlensed")
plt.loglog(ells, lensed_cls, 'b', label="Lensed")
legend = plt.legend(loc='lower right', shadow=True)
frame = legend.get_frame()
frame.set_facecolor('0.90')

#print len(theoretical_cls[:2750])
#print len(theoretical_lensed_cls)

deltas = (unlensed_cls - lensed_cls)/unlensed_cls
deltas_theory = (theoretical_cls[:2750] - theoretical_lensed_cls) / theoretical_cls[:2750]

#print deltas[3000:]
plt.figure()
plt.xlabel('l')
plt.ylabel('delta C_l / C_l')
plt.plot(ells, deltas, label="Lenspix delta")
plt.plot(theory_ells2, deltas_theory, '--r', label="Theory delta")
legend2 = plt.legend(loc="lower right", shadow=True)
frame2 = legend2.get_frame()
frame2.set_facecolor('0.90')
plt.show()
