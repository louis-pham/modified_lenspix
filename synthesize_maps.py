import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from matplotlib import cm
from astropy.io import fits

zInit=0.0
zFinal=4.6
n=22
dz=0.2
i=0

hot_cmap = cm.get_cmap("hot")
hot_cmap.set_under("w")

seismic_cmap = cm.get_cmap("seismic")
seismic_cmap.set_under("w")

zMin = i*dz
zMax = (i+1)*dz
zMinCib = (i+1)*dz
zMaxCib = (i+2)*dz


hdulist = fits.open("/scratch2/r/rbond/phamloui/lenspix_files/cib/lensed_cib_fullsky_ns2048_zmin0.20_zmax0.40_nu217_tot.fits")
nside = hdulist[1].header['NSIDE']

currLensedMap = np.zeros(hp.nside2npix(nside))
currUnlensedMap = np.zeros(hp.nside2npix(nside))

#currUnlensedAlmFilename = "/scratch2/r/rbond/phamloui/lenspix_files/cib/cib_fullsky_ns2048_zmin%.2f_zmax%.2f_nu217_tot_alm.fits" % (zMin, zMax)
#currLensedMapFilename = "/scratch2/r/rbond/phamloui/lenspix_files/cib/lensed_cib_fullsky_ns2048_zmin%.2f_zmax%.2f_nu217_tot.fits" % (zMinCib, zMaxCib)
#if zMax%1 == 0: #phi files have weird formatting
#    currPhiAlmFilename = "/scratch2/r/rbond/phamloui/lenspix_files/lensing/phi_z4pt6_z0.0_z%d_nside2048_hp.fits" % (int(zMax))
#else:
#    currPhiAlmFilename = "/scratch2/r/rbond/phamloui/lenspix_files/lensing/phi_z4pt6_z0.0_z%.1f_nside2048_hp.fits" % (zMax)
    
#print "unlensed"
#currUnlensedAlm = hp.read_alm(currUnlensedAlmFilename)
#print currUnlensedAlm.shape
#print "lensed"
#currLensedMap = np.nan_to_num(hp.read_map(currLensedMapFilename))
#print currLensedMap.shape
#print "phi"
#currPhiAlm = hp.read_alm(currPhiAlmFilename)
#print "lmax:", hp.Alm.getlmax(len(currUnlensedAlm))
#i = 1

plotArray = []
unlenPlotArray = []
while i<n:
    print "i:", i
    zMin = i*dz
    zMax = (i+1)*dz
    zMinCib = (i+1)*dz
    zMaxCib = (i+2)*dz
    nextUnlensedAlmFilename = "/scratch2/r/rbond/phamloui/lenspix_files/cib/cib_fullsky_ns2048_zmin%.2f_zmax%.2f_nu217_tot_alm.fits" % (zMinCib, zMaxCib)
    nextLensedMapFilename = "/scratch2/r/rbond/phamloui/lenspix_files/cib/lensed_cib_fullsky_ns2048_zmin%.2f_zmax%.2f_nu217_tot.fits" % (zMinCib, zMaxCib)
    if zMax%1 == 0:
        nextPhiAlmFilename = "/scratch2/r/rbond/phamloui/lenspix_files/lensing/phi_z4pt6_z0.0_z%d_nside2048_hp.fits" % (int(zMax))
    else:
        nextPhiAlmFilename = "/scratch2/r/rbond/phamloui/lenspix_files/lensing/phi_z4pt6_z0.0_z%.1f_nside2048_hp.fits" % (zMax)
    
    #print "unlensed"
    nextUnlensedAlm = hp.read_alm(nextUnlensedAlmFilename)
    nextUnlensedMap = hp.alm2map(nextUnlensedAlm, 2048)
    nextUnlensedMap = np.nan_to_num(nextUnlensedMap)
    #print nextUnlensedAlm.shape
    print "lensed"
    nextLensedMap = hp.read_map(nextLensedMapFilename)
    nextLensedMap = np.nan_to_num(nextLensedMap)
    #print nextLensedMap.shape
    #print "phi"
    #nextPhiAlm = hp.read_alm(nextPhiAlmFilename)
    #print "lmax:", hp.Alm.getlmax(len(nextUnlensedAlm))

    currLensedMap += nextLensedMap
    currUnlensedMap += nextUnlensedMap
    if i%7 == 0:
        print "Saving current step..."
        plotArray.append((zMaxCib, np.copy(currLensedMap)))
        unlenPlotArray.append((zMaxCib, np.copy(currUnlensedMap)))
        #hp.mollview(currLensedMap, xsize=2048, title="Step %d zmax=%.2f" % (i, zMaxCib), cmap=hot_cmap)
    i += 1

#hp.mollview(currLensedMap, xsize=2048, title="Step %d zmax=%.2f" % (i, zMaxCib))

maxTemp = max(plotArray[0][1].max(),plotArray[1][1].max(),plotArray[2][1].max(),plotArray[3][1].max(),unlenPlotArray[0][1].max(),unlenPlotArray[1][1].max(),unlenPlotArray[2][1].max(),unlenPlotArray[3][1].max())
minTemp = min(plotArray[0][1].min(),plotArray[1][1].min(),plotArray[2][1].min(),plotArray[3][1].min(),unlenPlotArray[0][1].min(),unlenPlotArray[1][1].min(),unlenPlotArray[2][1].min(),unlenPlotArray[3][1].min())

#print plotArray[0][1]
#print plotArray[1][1]
#print maxTemp
#print minTemp

for i in np.arange(0,4):
    z = plotArray[i][0]
    lenmap = plotArray[i][1]
    unlenmap = unlenPlotArray[i][1]
    diffMap = lenmap - unlenmap
    diffRange = np.absolute(diffMap).max()
    hp.mollview(lenmap, xsize=2048, title="Lensed summed to zmax=%.2f" % (z), cmap=hot_cmap, min=minTemp, max=maxTemp)
    hp.mollview(unlenmap, xsize=2048, title="Unlensed summed to zmax=%.2f" % (z), cmap=hot_cmap, min=minTemp, max=maxTemp)
    hp.mollview(diffMap, xsize=2048, title="Diff of summed to zmax=%.2f"%(z), cmap=seismic_cmap, min=-diffRange, max=diffRange)
plt.show()
