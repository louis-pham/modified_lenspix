import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import os.path

unlensed = [hp.alm2map(hp.read_alm("/scratch2/r/rbond/phamloui/lenspix_files/FromNERSC/ffp10_unlensed_scl_cmb_000_alm.fits", hdu=i), 2048, lmax=5120) for i in [1,2,3]]
T,Q,U = hp.read_map("/scratch2/r/rbond/phamloui/lenspix_files/output/jun10_nonlinear3_interp_factor3_julian_cmb.fits", field=(0,1,2))
lensed = [T,Q,U]

plt.figure()
plt.title("TEB Power Spectra")
plt.xlabel(r"$l$", fontsize=30)
plt.ylabel(r"$l(l+1)C_l/2\pi$", fontsize=30)
unlensed_line = ['r--','g--','m--']
lensed_line = ['b','y','k']
labels = ['T', 'E', 'B']
for i in [0,1,2]:
    print i
    unlensed_map = unlensed[i]
    lensed_map = lensed[i]
    unl_filename = "/scratch2/r/rbond/phamloui/lenspix_files/output/test_jun10_cl_unlensed_%s.fits" % (labels[i])
    len_filename = "/scratch2/r/rbond/phamloui/lenspix_files/output/test_jun10_cl_lensed_%s.fits" % (labels[i])

    if os.path.exists(unl_filename):
        print "reading unlensed cl..."
        unlensed_cl = hp.read_cl(unl_filename)
    else:
        print "unlensed anafast..."
        unlensed_cl = hp.anafast(unlensed_map, lmax=5120)    
        hp.write_cl(unl_filename, unlensed_cl)

    if os.path.exists(len_filename):
        print "reading lensed cl..."
        lensed_cl = hp.read_cl(len_filename)
    else:
        print "lensed anafast..."
        lensed_cl = hp.anafast(lensed_map, lmax=5120)    
        hp.write_cl(len_filename, lensed_cl)
    
    ell = np.arange(unlensed_cl.shape[0])

    unlensed_cl = ell * (ell+1) * unlensed_cl / (2 * np.pi)
    lensed_cl = ell * (ell+1) * lensed_cl / (2 * np.pi)

    #hp.mollview(unlensed_map, title="Unlensed %s" % (labels[i]), xsize=2048, unit="muK")
    plt.loglog(ell, unlensed_cl, unlensed_line[i], label="Unlensed %s" % (labels[i]))
    #hp.mollview(lensed_map, title="Lensed %s" % (labels[i]), xsize=2048, unit="muK")
    plt.loglog(ell, lensed_cl, lensed_line[i], label="Lensed %s" % (labels[i]))

legend2 = plt.legend(loc="lower left", shadow=True)
frame2 = legend2.get_frame()
frame2.set_facecolor('0.90')

plt.show()
