import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
from astropy.io import fits
import os

zInit=0.0
zFinal=4.6
n=22
dz=0.2
i=0

myCmap = cm.get_cmap("hot")
myCmap.set_under("w")
myCmap.set_over("w")

seismic_cmap = cm.get_cmap("seismic")
seismic_cmap.set_under("w")
seismic_cmap.set_over("w")

customCDict = {'red': ((0.0, 0.0, 0.0), 
                       (0.02, 0.3, 0.3),
                       (0.3, 1.0, 1.0),
                       (1.0, 1.0, 1.0)),
               'green': ((0.0, 0.0, 0.0),
                         (0.3, 0.0, 0.0), 
                         (0.7, 1.0, 1.0),
                         (1.0, 1.0, 1.0)),
               'blue': ((0.0, 0.0, 0.0), 
                        (0.7, 0.0, 0.0),
                        (1.0, 1.0, 1.0))
               }

customCmap = colors.LinearSegmentedColormap('hot2', customCDict, N=750, gamma=1.0)

unlensedAlmBaseFilename = "/scratch2/r/rbond/phamloui/lenspix_files/cib_v2_unlensed/cib_fullsky_ns2048_zmin%.2f_zmax%.2f_nu217_ns2048_tot_alm.fits"
lensedMapBaseFilename = "/scratch2/r/rbond/phamloui/lenspix_files/cib_v2_lensed/lensed_cib_fullsky_ns2048_zmin%.2f_zmax%.2f_nu217_ns2048_tot.fits"
kappaMapBaseFilename = "/scratch2/r/rbond/phamloui/lenspix_files/cib_v2_kappa/8Gpc_n2048_nb18_nt16_kap_sis_2_ns2048_zmin0.0_zmax%s_hp.fits" #placeholder %s for later format decision
phiAlmBaseFilename = "/scratch2/r/rbond/phamloui/lenspix_files/cib_v2_phi/8Gpc_n2048_nb18_nt16_phi_sis_2_ns2048_zmin0.0_zmax%s_hp.fits"

mapsSaved = True
savedUnlensedMapBaseFilename = "/scratch2/r/rbond/phamloui/lenspix_files/jun23_summed_unlensed_cib_v2_zmax_%.2f.fits"
savedLensedMapBaseFilename = "/scratch2/r/rbond/phamloui/lenspix_files/jun23_summed_lensed_cib_v2_zmax_%.2f.fits"
savedKappaMapBaseFilename = "/scratch2/r/rbond/phamloui/lenspix_files/jun23_summed_kappa_zmax_%.2f.fits"
savedPhiAlmBaseFilename = "/scratch2/r/rbond/phamloui/lenspix_files/jun23_summed_phi_zmax_%.2f.fits"

hdulist = fits.open(lensedMapBaseFilename % (0.2, 0.4))
nside = hdulist[1].header['NSIDE']
tUnit = "" # hdulist[1].header['TUNIT1']
#print tUnit
almSize = hp.read_alm(phiAlmBaseFilename % ('%d') % (1)).shape[0]
currLensedMap = np.zeros(hp.nside2npix(nside)) #initial maps to add to
currUnlensedMap = np.zeros(hp.nside2npix(nside))
currKappaMap = np.zeros(hp.nside2npix(nside))
currPhiAlm = np.zeros(almSize)
currPhiMap = np.zeros(hp.nside2npix(nside))
plotArray = []
unlenPlotArray = []
kappaArray = []
phiArray = []

if mapsSaved:
    for _z in [0.4, 1.8, 3.2, 4.6]:
        plotArray.append((_z, hp.read_map(savedLensedMapBaseFilename % _z)))
        unlenPlotArray.append((_z, hp.read_map(savedUnlensedMapBaseFilename % _z)))
        kappaArray.append((_z, hp.read_map(savedKappaMapBaseFilename % _z)))
        phiArray.append((_z, hp.read_alm(savedPhiAlmBaseFilename % _z)))
while i<n and not mapsSaved:
    print "i:", i
    zMin = i*dz
    zMax = (i+1)*dz
    zMinCib = (i+1)*dz
    zMaxCib = (i+2)*dz
    # load maps
    nextUnlensedAlmFilename = unlensedAlmBaseFilename % (zMinCib, zMaxCib)
    nextLensedMapFilename = lensedMapBaseFilename % (zMinCib, zMaxCib)
    if zMax%1 == 0: #phi files have weird formatting
        nextKappaMapFilename = kappaMapBaseFilename % ('%d') % (int(zMax))
        nextPhiAlmFilename = phiAlmBaseFilename % ('%d') % (int(zMax))
    else:
        nextKappaMapFilename = kappaMapBaseFilename % ('%.1f') % (zMax)
        nextPhiAlmFilename = phiAlmBaseFilename % ('%.1f') % (zMax)
    
    nextUnlensedAlm = hp.read_alm(nextUnlensedAlmFilename)
    print "num of nans (unlensed alm):", nextUnlensedAlm[np.isnan(nextUnlensedAlm)].shape[0]
    nextUnlensedMap = hp.alm2map(nextUnlensedAlm, 2048, verbose=False)
    print "num of nans (unlensed map):", nextUnlensedMap[np.isnan(nextUnlensedMap)].shape[0]
    # nextUnlensedMap = np.nan_to_num(nextUnlensedMap)

    nextLensedMap = hp.read_map(nextLensedMapFilename, verbose=False)
    print "num of nans (lensed map):", nextLensedMap[np.isnan(nextLensedMap)].shape[0]
    nextLensedMap = np.nan_to_num(nextLensedMap)

    nextKappaMap = hp.read_map(nextKappaMapFilename, verbose=False)
    print "num of nans (kappa map):", nextKappaMap[np.isnan(nextKappaMap)].shape[0]

    nextPhiAlm = hp.read_alm(nextPhiAlmFilename)
    print "num of nans (phi alm):", nextPhiAlm[np.isnan(nextPhiAlm)].shape[0]
    # nextPhiMap = hp.alm2map(nextPhiAlm, 2048)
    # print "num of nans (phi map):", nextPhiMap[np.isnan(nextPhiMap)].shape[0]
    
    currLensedMap += nextLensedMap
    currUnlensedMap += nextUnlensedMap
    currKappaMap += nextKappaMap
    currPhiAlm += nextPhiAlm
    # currPhiMap += nextPhiMap
    if i%7 == 0: #save every 7 iterations (should be 4 maps total)
        diffMap = nextLensedMap - nextUnlensedMap
        diffRange = np.absolute(diffMap).max()
        # hp.mollview(nextUnlensedMap, xsize=2048, title="Unlensed slice @ %.2f<zmax<%.2f" % (zMinCib, zMaxCib), cmap=myCmap)
        # hp.mollview(nextLensedMap, xsize=2048, title="Lensed slice @ %.2f<zmax<%.2f" % (zMinCib, zMaxCib), cmap=myCmap)
        # hp.mollview(nextKappaMap, xsize=2048, title="Kappa slice @ 0.00<zmax<%.2f" % zMax, cmap=myCmap)
        # hp.mollview(nextPhiMap, xsize=2048, title="Phi slice @ 0.00<zmax<%.2f" % zMax, cmap=myCmap)
        # hp.mollview(diffMap, xsize=2048, title="Diff slice @ %.2f<zmax<%.2f" % (zMinCib, zMaxCib), cmap=seismic_cmap, min=-diffRange, max=diffRange)
        # plt.show()
        print "Saving current step..."
        plotArray.append((zMaxCib, np.copy(currLensedMap)))
        unlenPlotArray.append((zMaxCib, np.copy(currUnlensedMap)))
        hp.write_map(savedLensedMapBaseFilename % zMaxCib, currLensedMap)
        hp.write_map(savedUnlensedMapBaseFilename % zMaxCib, currUnlensedMap)
        hp.write_map(savedKappaMapBaseFilename % zMaxCib, currKappaMap)
        hp.write_alm(savedPhiAlmBaseFilename % zMaxCib, currPhiAlm)
        # hp.mollview(currPhiMap, xsize=2048, title="Phi summed to zmax=%.2f" % zMaxCib, cmap=myCmap)
    i += 1

#hp.mollview(currLensedMap, xsize=2048, title="Step %d zmax=%.2f" % (i, zMaxCib))

#get colorbar ranges
maxTemp = max(plotArray[0][1].max(),plotArray[1][1].max(),plotArray[2][1].max(),plotArray[3][1].max(),unlenPlotArray[0][1].max(),unlenPlotArray[1][1].max(),unlenPlotArray[2][1].max(),unlenPlotArray[3][1].max())
minTemp = min(plotArray[0][1].min(),plotArray[1][1].min(),plotArray[2][1].min(),plotArray[3][1].min(),unlenPlotArray[0][1].min(),unlenPlotArray[1][1].min(),unlenPlotArray[2][1].min(),unlenPlotArray[3][1].min())
#print maxTemp
#print minTemp

# anotherMax = np.maximum(maxTemp, -minTemp) #test another colorbar range

divisor = 1. #play around with colorbar ranges

for i in np.arange(0,4):
    z = plotArray[i][0]
    lenmap = plotArray[i][1]
    # lenmap = hp.smoothing(lenmap, fwhm=0.001)
    unlenmap = unlenPlotArray[i][1]
    # unlenmap = hp.smoothing(unlenmap, fwhm=0.001)
    kappamap = kappaArray[i][1]
    # print kappamap.max()
    # print kappamap.min()
    phimap = hp.alm2map(phiArray[i][1],2048)
    # print phimap.max()
    # print phimap.min()

    diffMap = lenmap - unlenmap
    diffRange = np.absolute(diffMap).max()
    diffScheme = colors.SymLogNorm(linthresh=2.3e-20, linscale=0.3, vmin=-diffRange/divisor, vmax=diffRange/divisor)
    normScheme = colors.SymLogNorm(linthresh=4.1e-20, linscale=0.03, vmin=minTemp/divisor, vmax=maxTemp/divisor)
    # kappamap = hp.read_map(kappaMapBaseFilename % ('.1f') % (z))
    
    # hp.mollview(unlenmap, xsize=2048, title="Unlensed summed to zmax=%.2f" % (z), cmap=myCmap, min=minTemp/divisor, max=maxTemp/divisor, unit=tUnit, norm=normScheme)
    hp.mollview(kappamap, xsize=2048, title="Kappa summed to zmax=%.2f" % (z), cmap=myCmap)
    # hp.mollview(phimap, xsize=2048, title="Phi summed to zmax=%.2f" % (z), cmap=myCmap)
    # hp.mollview(lenmap, xsize=2048, title="Lensed summed to zmax=%.2f" % (z), cmap=myCmap, min=minTemp/divisor, max=maxTemp/divisor, unit=tUnit, norm=normScheme)
    # hp.mollview(diffMap, xsize=2048, title="Diff of summed to zmax=%.2f"%(z), cmap=seismic_cmap, min=-diffRange/divisor, max=diffRange/divisor, unit=tUnit, norm=diffScheme)
    
# plot lensed map as if it were single redshift
# summedLensedMap = hp.read_map("/scratch2/r/rbond/phamloui/lenspix_files/output/jun14_cib_test_summed_slice_lensed.fits")
# hp.mollview(summedLensedMap, xsize=2048, title="Lensed single redhsift zmax=4.60", cmap=myCmap, min=minTemp/divisor, max=maxTemp/divisor, unit=tUnit)
# summedLensedMap = np.nan_to_num(summedLensedMap)
# summedLensedCls = hp.anafast(summedLensedMap)
# summedLensedElls = np.arange(summedLensedCls.shape[0])
# summedLensedCls = summedLensedElls * (summedLensedElls + 1) * summedLensedCls / (2*np.pi)

# plt.figure()
# plt.title("Cib Power Spectra")
# plt.xlabel(r"$l$")
# plt.ylabel(r"$l(l+1)C_l/2\pi$")
# colours = ['#EA638C','#4B3F72','#0AFFED','#1F2041']
# plt.xlabel(r'$l$')
# plt.ylabel(r"$l(l+1)C_l/2\pi$")

# testcmbkappamap = hp.read_map("/scratch2/r/rbond/phamloui/lenspix_files/kappa_maps/jun1_nonlinear3_kappa_for_julian_cmb.fits")
# testcmbphialm = hp.read_alm("/scratch2/r/rbond/phamloui/lenspix_files/output/jun10_nonlinear3_interp_factor3_phi_for_julian_cmb.fits")
# testcmbphicl = hp.alm2cl(testcmbphialm)[1:] 
# testell = np.arange(0, testcmbphicl.shape[0]+1)[1:]

# cmb_deflection = np.sqrt(np.sum((2 * testell + 1) / (4 * np.pi) * testell * (testell + 1) * testcmbphicl)) * (180*60/np.pi)
# print "CMB RMS deflection angle:", cmb_deflection

# kapClFile = "/scratch2/r/rbond/phamloui/lenspix_files/jun22_summed_kappa_v2_zmax_%.2f_cl.dat"
# phiClFile = "/scratch2/r/rbond/phamloui/lenspix_files/jun22_summed_phi_v2_zmax_%.2f_cl.dat"

# print "getting power spectra of slices"
# for i in np.arange(0,4):
#     print "i:",i
#     z = plotArray[i][0]
#     curmaplen = plotArray[i][1]
#     # curmapunl = unlenPlotArray[i][1]
#     # curmapkap = kappaArray[i][1]
#     # curalmphi = phiArray[i][1]

#     curcllen = hp.anafast(curmaplen)
#     # if not os.path.exists(kapClFile % z): 
#     #     curclkap = hp.anafast(curmapkap)
#     #     hp.write_cl(kapClFile % z, curclkap)
#     # else:
#     #     curclkap = hp.read_cl(kapClFile % z)
#     # if not os.path.exists(phiClFile % z): 
#     #     curclphi = hp.alm2cl(curalmphi)
#     #     hp.write_cl(phiClFile % z, curclphi)
#     # else:
#     #     curclphi = hp.read_cl(phiClFile % z)

#     ell = np.arange(curcllen.shape[0])

#     # RMS_deflection = np.sqrt(np.sum((2 * ell + 1) / (4 * np.pi) * ell * (ell + 1) * curclphi)) * (180*60/np.pi)
#     # print "RMS deflection angle @ zmax=%.2f:" % z, RMS_deflection
#     # print "CIB Phi RMS vs CMB Phi RMS:", RMS_deflection / cmb_deflection
#     curcllen = ell * (ell+1) * curcllen / (2 * np.pi)
#     # curclkap = ell * (ell+1) * curclkap / (2 * np.pi)
#     # curclphi = ell * (ell+1) * curclphi / (2 * np.pi)
#     plt.loglog(ell, curcllen, c=colours[i], label="zmax %.2f" % (z))
#     # plt.loglog(ell, curclkap, c=colours[i], label="zmax %.2f" % (z))
#     # plt.loglog(ell, curclphi, c=colours[i], label="phi zmax %.2f" % (z))
#     # if i == 3:
#     #     finalCl = curcl

# plt.loglog(summedLensedElls, summedLensedCls, c="#F9564F", label="zmax 4.60 single redshift")

# legend = plt.legend(loc="lower right", shadow=True)
# frame = legend.get_frame()
# frame.set_facecolor('0.90')

# finalDiff = (finalCl - summedLensedCls) / summedLensedCls
# plt.figure()
# plt.title("Power Spectra Difference at zmax=4.60")
# plt.xlabel(r"$l$")
# plt.ylabel(r"$\Delta C_l/C_l$")
# plt.plot(summedLensedElls, finalDiff, lw="1", c="#a6cee3", ls="-")

plt.show()
