import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
from astropy.io import fits

zInit=0.0
zFinal=4.6
n=22
dz=0.2
i=0

# if you've run the script before, set this to True to avoid doing summation again
mapsSaved = False

myCmap = cm.get_cmap("hot")
myCmap.set_under("w")

seismic_cmap = cm.get_cmap("seismic")
seismic_cmap.set_under("w")

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

customCmap = colors.LinearSegmentedColormap('hot2', customCDict, N=500, gamma=1.0)

# zMin = i*dz
# zMax = (i+1)*dz
# zMinCib = (i+1)*dz
# zMaxCib = (i+2)*dz

hdulist = fits.open("/scratch2/r/rbond/phamloui/lenspix_files/cib/lensed_cib_fullsky_ns2048_zmin0.20_zmax0.40_nu217_tot.fits")
nside = hdulist[1].header['NSIDE']
tUnit = hdulist[1].header['TUNIT1']
#print tUnit
currLensedMap = np.zeros(hp.nside2npix(nside)) #initial maps to add to
currUnlensedMap = np.zeros(hp.nside2npix(nside))

plotArray = []
unlenPlotArray = []

if mapsSaved:
    for _z in [0.4, 1.8, 3.2, 4.6]:
        plotArray.append((_z, hp.read_map("/scratch2/r/rbond/phamloui/lenspix_files/test_jun10_summed_lensed_cib_to_zmax_%.2f.fits" % _z)))
        unlenPlotArray.append((_z, hp.read_map("/scratch2/r/rbond/phamloui/lenspix_files/test_jun10_summed_unlensed_cib_to_zmax_%.2f.fits" % _z)))

while i<n and not mapsSaved:
    print "i:", i
    zMin = i*dz
    zMax = (i+1)*dz
    zMinCib = (i+1)*dz
    zMaxCib = (i+2)*dz
    # load maps
    nextUnlensedAlmFilename = "/scratch2/r/rbond/phamloui/lenspix_files/cib/cib_fullsky_ns2048_zmin%.2f_zmax%.2f_nu217_tot_alm.fits" % (zMinCib, zMaxCib)
    nextLensedMapFilename = "/scratch2/r/rbond/phamloui/lenspix_files/cib/lensed_cib_fullsky_ns2048_zmin%.2f_zmax%.2f_nu217_tot.fits" % (zMinCib, zMaxCib)
    if zMax%1 == 0: #phi files have weird formatting
        nextPhiAlmFilename = "/scratch2/r/rbond/phamloui/lenspix_files/lensing/phi_z4pt6_z0.0_z%d_nside2048_hp.fits" % (int(zMax))
    else:
        nextPhiAlmFilename = "/scratch2/r/rbond/phamloui/lenspix_files/lensing/phi_z4pt6_z0.0_z%.1f_nside2048_hp.fits" % (zMax)
    
    nextUnlensedAlm = hp.read_alm(nextUnlensedAlmFilename)
    nextUnlensedMap = hp.alm2map(nextUnlensedAlm, 2048)
    nextUnlensedMap = np.nan_to_num(nextUnlensedMap)

    nextLensedMap = hp.read_map(nextLensedMapFilename)
    nextLensedMap = np.nan_to_num(nextLensedMap)

    #nextPhiAlm = hp.read_alm(nextPhiAlmFilename)

    currLensedMap += nextLensedMap
    currUnlensedMap += nextUnlensedMap
    if i%7 == 0: #save every 7 iterations (should be 4 maps total)
        print "Saving current step..."
        plotArray.append((zMaxCib, np.copy(currLensedMap)))
        unlenPlotArray.append((zMaxCib, np.copy(currUnlensedMap)))
        hp.write_map("/scratch2/r/rbond/phamloui/lenspix_files/test_jun10_summed_lensed_cib_to_zmax_%.2f.fits" % zMaxCib, currLensedMap)
        hp.write_map("/scratch2/r/rbond/phamloui/lenspix_files/test_jun10_summed_unlensed_cib_to_zmax_%.2f.fits" % zMaxCib, currUnlensedMap)
        #hp.mollview(currLensedMap, xsize=2048, title="Step %d zmax=%.2f" % (i, zMaxCib), cmap=myCmap)
    i += 1

#hp.mollview(currLensedMap, xsize=2048, title="Step %d zmax=%.2f" % (i, zMaxCib))

# get colorbar ranges
maxTemp = max(plotArray[0][1].max(),plotArray[1][1].max(),plotArray[2][1].max(),plotArray[3][1].max(),unlenPlotArray[0][1].max(),unlenPlotArray[1][1].max(),unlenPlotArray[2][1].max(),unlenPlotArray[3][1].max())
minTemp = min(plotArray[0][1].min(),plotArray[1][1].min(),plotArray[2][1].min(),plotArray[3][1].min(),unlenPlotArray[0][1].min(),unlenPlotArray[1][1].min(),unlenPlotArray[2][1].min(),unlenPlotArray[3][1].min())]
#print maxTemp
#print minTemp

anotherMax = np.maximum(maxTemp, -minTemp) #test another colorbar range

normScheme = colors.SymLogNorm(linthresh=1e-20, linscale=0.03, vmin=-anotherMax, vmax=anotherMax)

for i in np.arange(3,4):
    z = plotArray[i][0]
    lenmap = plotArray[i][1]
    unlenmap = unlenPlotArray[i][1]
    diffMap = lenmap - unlenmap
    diffRange = np.absolute(diffMap).max()
    diffScheme = colors.SymLogNorm(linthresh=1e-20, linscale=0.03, vmin=-diffRange, vmax=diffRange)
    hp.mollview(lenmap, xsize=2048, title="Lensed summed to zmax=%.2f" % (z), cmap=myCmap, min=minTemp, max=maxTemp, norm=normScheme, unit=tUnit)
    hp.mollview(unlenmap, xsize=2048, title="Unlensed summed to zmax=%.2f" % (z), cmap=myCmap, min=minTemp, max=maxTemp, norm=normScheme, unit=tUnit)
    hp.mollview(diffMap, xsize=2048, title="Diff of summed to zmax=%.2f"%(z), cmap=seismic_cmap, min=-diffRange, max=diffRange, norm=normScheme, unit=tUnit)
plt.show()

# plt.figure()
# plt.title("CIB Power Spectra")
# plt.xlabel(r"$l$")
# plt.ylabel(r"$l(l+1)C_l/2\pi$")
# colours = ['r','b','m','k']
# plt.xlabel(r'$C_l$')
# plt.ylabel('Count')

# for i in np.arange(3,4):
    # z = plotArray[i][0]
    # curmap = plotArray[i][1]
    # curmapunl = unlenPlotArray[i][1]
    # curcl = hp.anafast(curmap)
    # ell = np.arange(curcl.shape[0])
    # curcl = ell * (ell+1) * curcl / (2 * np.pi)
    # plt.loglog(ell, curcl, colours[i], label="zmax %.2f" % (z))
    
# legend = plt.legend(loc="lower right", shadow=True)
# frame = legend.get_frame()
# frame.set_facecolor('0.90')
# plt.show()
