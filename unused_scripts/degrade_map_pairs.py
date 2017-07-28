# loops over all cib maps and downgrades them to desired nside
import numpy as np
import healpy as hp
# from kappa2phi import kap2phi                                                                       

def appendToFilename(filename, newStr):
    return "{0}_{2}.{1}".format(*filename.rsplit('.', 1) + [newStr])

zInit=0.0
zFinal=4.6
n=22
dz=0.2
i=0
# unlensedPrimary = hp.read_alm("/scratch2/r/rbond/phamloui/lenspix_files/cib_n4096_unlensed/cib_fullsky_ns4096_zmin0.20_zmax0.40_nu217_ns4096_tot_alm.fits")
while i<n:
    print "i:", i
    zMin = i*dz
    zMax = (i+1)*dz
    zMinCib = (i+1)*dz
    zMaxCib = (i+2)*dz

    #primary map                                                                                               
    mapFilename = "/scratch2/r/rbond/phamloui/lenspix_files/cib_v2_unlensed/cib_fullsky_ns4096_zmin%.2f_zmax%.2f_nu217_ns4096_tot.fits" % (zMinCib, zMaxCib)                                                                       
    downgradeFile = "/scratch2/r/rbond/phamloui/lenspix_files/cib_v2_unlensed/cib_fullsky_ns2048_zmin%.2f_zmax%.2f_nu217_ns2048_tot.fits" % (zMinCib, zMaxCib)
    print downgradeFile
    curmap = hp.read_map(mapFilename, verbose=False)
    downgraded = hp.ud_grade(curmap, 2048)
    hp.write_map(downgradeFile, downgraded)

    #kappa map
    if zMax%1 == 0: #phi files have weird formatting                                                        
        kappaMapFilename = "/scratch2/r/rbond/phamloui/lenspix_files/cib_v2_kappa/8Gpc_n4096_nb18_nt16_kap_sis_2_ns4096_zmin0.0_zmax%d_hp.fits" % (int(zMax))
        downgradeKapFile = "/scratch2/r/rbond/phamloui/lenspix_files/cib_v2_kappa/8Gpc_n2048_nb18_nt16_kap_sis_2_ns2048_zmin0.0_zmax%d_hp.fits" % (int(zMax))
    else:
        kappaMapFilename = "/scratch2/r/rbond/phamloui/lenspix_files/cib_v2_kappa/8Gpc_n4096_nb18_nt16_kap_sis_2_ns4096_zmin0.0_zmax%.1f_hp.fits" % (zMax)
        downgradeKapFile = "/scratch2/r/rbond/phamloui/lenspix_files/cib_v2_kappa/8Gpc_n2048_nb18_nt16_kap_sis_2_ns2048_zmin0.0_zmax%.1f_hp.fits" % (zMax)
    print downgradeKapFile
    curkapmap = hp.read_map(kappaMapFilename, verbose=False)
    downgradedKap = hp.ud_grade(curkapmap, 2048)
    hp.write_map(downgradeKapFile, downgradedKap)

    i += 1
