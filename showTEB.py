# plot unlensed and lensed TEB for CMB
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import os.path

print "Getting unlensed cls..."
unlensedCls = hp.alm2cl([hp.read_alm("/scratch2/r/rbond/phamloui/lenspix_files/FromNERSC/ffp10_unlensed_scl_cmb_000_alm.fits", hdu=i) for i in [1,2,3]])
lmax = unlensedCls[0].shape[0]-1
print lmax
ell = np.arange(0,lmax+1)
#[unlensedTT, unlensedEE, unlensedBB] = [unlensedCls[0], unlensedCls[1], unlensedCls[2]]

print "Getting lensed TQU map..."
lensedT,lensedQ,lensedU = hp.read_map("/scratch2/r/rbond/phamloui/lenspix_files/output/julian_lensed/jun10_nonlinear3_interp_factor3_julian_cmb.fits", field=(0,1,2))
filenameBase = '/scratch2/r/rbond/phamloui/lenspix_files/output/julian_lensed/jun16_nonlinear3_interp_factor3_julian_lensed_%s_cl.dat'

print "Getting lensed theory..."
theoryClFile = "/scratch2/r/rbond/phamloui/lenspix_files/FromNERSC/cls/ffp10_lensedCls.dat"
theory_TT_ell, theory_CL_TT, theory_CL_EE, theory_CL_BB = np.loadtxt(theoryClFile, usecols=(0,1,2,3), unpack=True)
theory_TT_ell = np.insert(theory_TT_ell, 0, 1); theory_TT_ell = np.insert(theory_TT_ell, 0, 0)
theory_CL_TT = np.insert(theory_CL_TT, 0, 0); theory_CL_TT = np.insert(theory_CL_TT, 0, 0)
theory_CL_EE = np.insert(theory_CL_EE, 0, 0); theory_CL_EE = np.insert(theory_CL_EE, 0, 0)
theory_CL_BB = np.insert(theory_CL_BB, 0, 0); theory_CL_BB = np.insert(theory_CL_BB, 0, 0)
theory_CL_TT = theory_CL_TT / 1e12
theory_CL_EE = theory_CL_EE / 1e12
theory_CL_BB = theory_CL_BB / 1e12
theoryCls = [theory_CL_TT, theory_CL_EE, theory_CL_BB]
# print "lensed nans:", lensedT[np.isnan(lensedT)].shape[0]

spectrumIds = ['TT','EE','BB']
idIndices = np.arange(0,len(spectrumIds))

if os.path.exists(filenameBase % ('TT')): #load spectra if saved
    print "Reading lensed cls from file..."
    lensedCls = [hp.read_cl(filenameBase % (_id)) for _id in spectrumIds]
else:
    print "Getting lensed cls from map..."
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
lensed_theory_linecol = ['#89BD9E','#DB222A','#3454D1']
labels = ['TT', 'EE', 'BB']
for i in idIndices:
    print i
    curUnlensed = unlensedCls[i]
    curLensed = lensedCls[i]
    curTheoryLensed = theoryCls[i]
    curUnlensed = ell * (ell + 1) * curUnlensed / (2 * np.pi)
    curLensed = ell * (ell + 1) * curLensed / (2 * np.pi)
    # curTheoryLensed = theory_TT_ell * (theory_TT_ell + 1) * curTheoryLensed / (2 * np.pi)
    plt.loglog(ell, curUnlensed, c=unlensed_linecol[i], label="Unlensed %s" % (labels[i]))
    plt.loglog(ell, curLensed, c=lensed_linecol[i], label="Lensed %s" % (labels[i]))
    plt.loglog(theory_TT_ell,  curTheoryLensed, c=lensed_linecol[i], ls='--', label="Theory Lensed %s" % (labels[i]))

legend2 = plt.legend(loc="lower left", shadow=True)
frame2 = legend2.get_frame()
frame2.set_facecolor('0.90')

plt.show()
