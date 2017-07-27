# file for looping over CIB maps - mainly used for debugging and plotting
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import colormaps as cmaps
import os.path
import imageio
# from kappa2phi import kap2phi

def appendToFilename(filename, newStr):
    return "{0}_{2}.{1}".format(*filename.rsplit('.', 1) + [newStr])

def printNanAmount(_map, s=''):
    if s:
        s = ' (%s)' % s
    print "num of nans%s:" % s, _map[np.isnan(_map)].shape[0]
    print _map.min()
    print _map.max()

zInit=0.0
zFinal=4.6
n=2
dz=0.2
i=1
nside = 2048
zeroMap = np.zeros(12 * (nside**2))

# manually load colormaps from later version of numpy
plt.register_cmap(name='viridis', cmap=cmaps.viridis)
plt.register_cmap(name='magma', cmap=cmaps.magma)
plt.register_cmap(name='inferno', cmap=cmaps.inferno)
plt.register_cmap(name='plasma', cmap=cmaps.plasma)
hotCmap = cm.get_cmap('hot')
hotCmap.set_under('w')
viridisCmap = cm.get_cmap('viridis')
viridisCmap.set_under('w')
magmaCmap = cm.get_cmap('magma')
magmaCmap.set_under('w')
seismicCmap = cm.get_cmap('seismic')
seismicCmap.set_under('w')

cibCmap = viridisCmap
lensingCmap = hotCmap
diffCmap = seismicCmap

_fwhm = 0.0001 #how much to smooth by
rotAngle = (2.356, 0.812, 0.) #change map centre
_range = [-10.,10.]

# took min and max of all maps and saved them in these vars
# diffAbsMax = 5.98130261969e-21
# cibAbsMin = 4.54899993655e-24
# cibAbsMax = 5.22722053948e-20

# curSummedKappa = np.copy(zeroMap) # used for map summing
while i<n:
    print "i:", i
    # if i%3 != 0 and i != 1:
    # if i%3 == 0 or i == 1:
    #     i += 1
    #     continue

    zMin = i*dz
    zMax = (i+1)*dz
    zMinCib = (i+1)*dz
    zMaxCib = (i+2)*dz

    mapFilename = "/scratch2/r/rbond/phamloui/lenspix_files/cib_v2_unlensed/cib_fullsky_ns2048_zmin%.2f_zmax%.2f_nu217_ns2048_tot_fwhm_0.0035.fits" % (zMinCib, zMaxCib)
    lensedMapFilename = "/scratch2/r/rbond/phamloui/lenspix_files/cib_v2_unlensed/cib_fullsky_ns2048_zmin%.2f_zmax%.2f_nu217_ns2048_tot.fits" % (zMinCib, zMaxCib)
    # lensedMapFilename = "/scratch2/r/rbond/phamloui/lenspix_files/cib_lensed_kap0.2/jul13_cib_fullsky_ns2048_zmin%.2f_zmax%.2f_nu217_ns2048_tot_fwhm_0.0035.fits" % (zMinCib, zMaxCib)
    almFilename = appendToFilename(mapFilename, 'alm')
    if zMax%1 == 0: #phi files have weird formatting                                                           
        kappaMapFilename = "/scratch2/r/rbond/phamloui/lenspix_files/cib_kappa_zmin0.2/8Gpc_n4096_nb18_nt16_kap_sis_2_ns2048_zmin0.2_zmax%d_hp.fits" % (int(zMax))
        phiAlmFilename = "/scratch2/r/rbond/phamloui/lenspix_files/cib_phi_zmin0.2/8Gpc_n4096_nb18_nt16_phi_sis_2_ns2048_zmin0.2_zmax%d_hp.fits" % (int(zMax))
        oldPhiFile = "/scratch2/r/rbond/phamloui/lenspix_files/cib_v2_phi/8Gpc_n2048_nb18_nt16_phi_sis_2_ns2048_zmin0.0_zmax%d_hp.fits" % (int(zMax))
    else:
        kappaMapFilename = "/scratch2/r/rbond/phamloui/lenspix_files/cib_kappa_zmin0.2/8Gpc_n4096_nb18_nt16_kap_sis_2_ns2048_zmin0.2_zmax%.1f_hp.fits" % (zMax)
        phiAlmFilename = "/scratch2/r/rbond/phamloui/lenspix_files/cib_phi_zmin0.2/8Gpc_n4096_nb18_nt16_phi_sis_2_ns2048_zmin0.2_zmax%.1f_hp.fits" % (zMax)
        oldPhiFile = "/scratch2/r/rbond/phamloui/lenspix_files/cib_v2_phi/8Gpc_n2048_nb18_nt16_phi_sis_2_ns2048_zmin0.0_zmax%.1f_hp.fits" % (zMax)
    phiMapFilename = appendToFilename(phiAlmFilename, 'map')
    oldPhiMapFile = appendToFilename(oldPhiFile, 'map')
    # unlensedSmoothedFilename = appendToFilename(mapFilename, 'fwhm_%s'%_fwhm)
    # lensedSmoothedFilename = appendToFilename(lensedMapFilename, 'fwhm_%s'%_fwhm)
    # # for primary map -> alm
    print "read unlensed"
    unlensedMap = hp.read_map(mapFilename)
    # printNanAmount(unlensedMap, 'unlensed')
    # # unlensedMapRaw = np.copy(unlensedMap)
    # if os.path.exists(unlensedSmoothedFilename):
    #     unlensedMap = hp.read_map(unlensedSmoothedFilename)
    # else:
    #     unlensedMap = hp.smoothing(unlensedMap, fwhm=_fwhm)
    #     hp.write_map(unlensedSmoothedFilename, unlensedMap)
    # unlensedAlm = hp.map2alm(unlensedMap)
    # hp.write_alm(almFilename, unlensedAlm)

    # unlensedAlm = hp.read_alm(appendToFilename(appendToFilename(mapFilename, 'd_1e3'), 'alm'))
    # printNanAmount(unlensedAlm, 'unlensed alm d1e3')
    
    # #for kappa map -> phi alm
    # print "read kappa"
    # kappaMap = hp.read_map(kappaMapFilename)
    # printNanAmount(unlensedMap, 'kappa')
    # kappaAlm = hp.map2alm(kappaMap)
    # # kap2phi(zeroMap, kappaMap, unlensedPrimary, phiAlmFilename)
    # lmax = hp.Alm.getlmax(kappaAlm.shape[0])
    # # convert to phi (grav potential)        
    # print "----Converting kappa to phi..."
    # l,m = hp.Alm.getlm(lmax)
    # phiAlm = kappaAlm * (2.0 / (l*(l+1.0)))
    # phiAlm[l==0] = 0
    # print "----Writing phi alm to file..." 
    # hp.write_alm(phiAlmFilename, phiAlm)
    # print "----Done."

    # phiMap = hp.read_map(phiMapFilename)
    # printNanAmount(phiMap, 'phi')
    
    print 'read lensed'
    lensedMap = hp.read_map(lensedMapFilename)
    printNanAmount(lensedMap, "%.1f<z<%.1f" % (zMinCib, zMaxCib))
    lensedUnseen = np.copy(lensedMap)
    lensedUnseen[np.isnan(lensedUnseen)] = hp.UNSEEN
    # printNanAmount(lensedMap, 'lensed d_1e3')
    lensedNanToZero = np.nan_to_num(lensedMap)

    # if os.path.exists(lensedSmoothedFilename):
    #     lensedMap = hp.read_map(lensedSmoothedFilename)
    # else:
    #     lensedMap = hp.smoothing(lensedMap, fwhm=_fwhm)
    #     hp.write_map(lensedSmoothedFilename, lensedMap)
    
    diffMap = lensedMap - unlensedMap
    diffUnseen = np.copy(diffMap)
    diffUnseen[np.isnan(diffMap)] = hp.UNSEEN
    diffNanToZero = np.nan_to_num(diffMap)

    # maxTemp = unlensedMap.max() / divisor
    # minTemp = unlensedMap.min() / divisor   
    unlensedMin = unlensedMap.min()
    unlensedMax = unlensedMap.max()
    lensedMin = lensedNanToZero.min()
    lensedMax = lensedNanToZero.max()

    maxTemp = max(lensedMax, unlensedMax)
    minTemp = min(lensedMin, unlensedMin)

    # if minTemp > cibAbsMin:
    #     cibAbsMin = minTemp
    # if maxTemp > cibAbsMax:
    #     cibAbsMax = maxTemp
    # maxTemp = cibAbsMax
    # minTemp = cibAbsMin

    diffMax = np.absolute(diffNanToZero).max()
    # diffMax = diffAbsMax
    # if diffMax > diffAbsMax:
    #     diffAbsMax = diffMax

    # if os.path.exists(phiMapFilename):
    #     print "read phi map"
    #     phiMap = hp.read_map(phiMapFilename)
    # else:
    #     print "read phi alm"
    #     phiAlm = hp.read_alm(phiAlmFilename)
    #     phiMap = hp.alm2map(phiAlm, 2048)
    #     hp.write_map(phiMapFilename, phiMap)
    
    # diffScheme = colors.SymLogNorm(linthresh=9.e-5, linscale=0.03, vmin=-diffMax, vmax=diffMax)
    # normScheme = colors.SymLogNorm(linthresh=1.5e-8, linscale=1., vmin=minTemp, vmax=maxTemp)
    
    print 'plot unlensed'
    hp.mollzoom(unlensedMap, xsize=2048, cmap=cibCmap, title="Unlensed @ %.2f<zmax<%.2f smoothed fwhm=0.0035" % (zMinCib,zMaxCib), min=unlensedMin, max=unlensedMax, rot=rotAngle)#, norm=normScheme)
    # hp.set_g_clim(minTemp, maxTemp) #sets min/max of zoomed portion
    hp.set_g_clim(unlensedMin, unlensedMax)
    # plt.savefig('/scratch2/r/rbond/phamloui/pics/jul17_unlensed_fwhm_0.0035/unlensed_zmin%.1f_zmax%.1f.png' % (zMinCib, zMaxCib), dpi=220) 
    print 'plot lensed'
    hp.mollzoom(lensedUnseen, xsize=2048, cmap=cibCmap, title="Unlensed @ %.2f<zmax<%.2f" % (zMinCib,zMaxCib), min=lensedMin, max=lensedMax, rot=rotAngle)#, norm=normScheme)
    # hp.set_g_clim(minTemp, maxTemp)
    hp.set_g_clim(lensedMin, lensedMax)
    # plt.savefig('/scratch2/r/rbond/phamloui/pics/jul17_lensed_fwhm_0.0035/lensed_zmin%.1f_zmax%.1f.png' % (zMinCib, zMaxCib), dpi=220)
    # print 'plot diff'
    # hp.mollzoom(diffUnseen, xsize=2048, cmap=diffCmap, title="Diff @ %.2f<zmax<%.2f" % (zMinCib,zMaxCib), min=-diffMax, max=diffMax)#, norm=diffScheme)
    # hp.set_g_clim(-diffMax, diffMax)
    # plt.savefig('/scratch2/r/rbond/phamloui/pics/jul17_diff_fwhm_0.0035/diff_zmin%.1f_zmax%.1f.png' % (zMinCib, zMaxCib), dpi=220)
    # print 'plot kappa'
    # hp.mollzoom(kappaMap, xsize=2048, cmap=lensingCmap, title="Kappa @ 0.20<zmax<%.2f" % zMax)
    # print 'plot phi'
    # hp.mollview(phiMap, xsize=2048, cmap=lensingCmap, title="Phi @ 0.20<zmax<%.2f" % zMax, min=phiMin, max=phiMax)
    # hp.set_g_clim(phiMin, phiMax)

    i += 1
    plt.show()

# print "diffAbsMax =", diffAbsMax
# print "cibAbsMin =", cibAbsMin
# print "cibAbsMax =", cibAbsMax
