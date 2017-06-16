import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import os.path

print "Getting unlensed cls..."
unlensedCls = hp.alm2cl([hp.read_alm("/scratch2/r/rbond/phamloui/lenspix_files/FromNERSC/ffp10_unlensed_scl_cmb_000_alm.fits", hdu=i) for i in [1,2,3]])
lmax = unlensedCls[0].shape[0]-1
ell = np.arange(0,lmax+1)
#[unlensedTT, unlensedEE, unlensedBB] = [unlensedCls[0], unlensedCls[1], unlensedCls[2]]

print "Getting lensed cls..."
lensedT,lensedQ,lensedU = hp.read_map("/scratch2/r/rbond/phamloui/lenspix_files/output/jun10_nonlinear3_interp_factor3_julian_cmb.fits", field=(0,1,2))
filenameBase = '/scratch2/r/rbond/phamloui/lenspix_files/output/jun16_julian_lensed_%s_cl.dat'
spectrumIds = ['TT','EE','BB']
if os.path.exists(filenameBase % ('TT')): #load spectra if saved
    print "Reading lensed cls from file..."
    lensedCls = [hp.read_cl(filenameBase % (_id)) for _id in spectrumIds]
else:
    lensedCls = hp.anafast([lensedT, lensedQ, lensedU], lmax=lmax)
    print "Writing lensed cls to file..."
    for _id in spectrumIds:
        curCl = lensedCls[spectrumIds.index(_id)]
        hp.write_cl(filenameBase % (_id), curCl)
#[lensedTT, lensedEE, lensedBB] = [lensedCls[0], lensedCls[1], lensedCls[2]]

print "Plotting..."
plt.figure()
plt.title("TEB Power Spectra")
plt.xlabel(r"$l$", fontsize=30)
plt.ylabel(r"$l(l+1)C_l/2\pi$", fontsize=30)
unlensed_linecol = ['#053C5E','#FFCF99','#000000'] # BB is all zero so colour doesnt matter
lensed_linecol = ['#89BD9E','#DB222A','#3454D1']
labels = ['TT', 'EE', 'BB']
for i in [0,1,2]:
    print i
    curUnlensed = unlensedCls[i]
    curLensed = lensedCls[i]
    curUnlensed = ell * (ell + 1) * curUnlensed / (2 * np.pi)
    curLensed = ell * (ell + 1) * curLensed / (2 * np.pi)
    plt.loglog(ell, curUnlensed, c=unlensed_linecol[i], label="Unlensed %s" % (labels[i]))
    plt.loglog(ell, curLensed, c=lensed_linecol[i], label="Lensed %s" % (labels[i]))

legend2 = plt.legend(loc="lower left", shadow=True)
frame2 = legend2.get_frame()
frame2.set_facecolor('0.90')

plt.show()
