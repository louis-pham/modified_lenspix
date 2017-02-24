import matplotlib.pyplot as plt
import numpy as np
import healpy as hp

def kap2phi(writeMap=False):
    nside = 4096
    lmax = 5120
    comb_kappa_file = "kappa_maps/8Gpc_n4096_nb23_nt18_kap_comb.fits"
    phi_alm_file = "kappa_maps/comb_z_4.6_n" + str(nside) + "_lmax" + str(lmax) + "_phi_alm.fits"
    phi_map_file = "kappa_maps/comb_z_4.6_n" + str(nside) + "_lmax" + str(lmax) + "_phi_map.fits"
    
    print "===============kappa2phi==============="
    field_kappa_file = "kappa_maps/8Gpc_n4096_nb23_nt18_kap_field.fits"
    halo_kappa_file = "kappa_maps/8Gpc_n4096_nb23_nt18_kap_halo.fits"
    
    print "----Loading maps..."
    field_kappa_map = hp.read_map(field_kappa_file)
    halo_kappa_map = hp.read_map(halo_kappa_file) 
    print "----Done."

    #combine
    print "----Combining maps..." #convert to alm before combining?
    kappa_map = field_kappa_map + halo_kappa_map
    print "----Done."
    print "----Writing new combined map to file..."
    hp.write_map(comb_kappa_file, kappa_map)
    print "----Done."
    
    #no field kappa right now, so just override
    #kappa_map = halo_kappa_map
    
    #convert to alm
    print "----Converting kappa map to alm..."
    kappa_lm = hp.map2alm(kappa_map, lmax=lmax)
    print "----Done."
    
    #convert to phi (grav potential)
    print "----Converting kappa to phi..."
    l,m = hp.Alm.getlm(lmax)
    phi_lm = kappa_lm * (2.0 / (l*(l+1.0)))
    phi_lm[l==0] = 0
    
    print "----Writing phi alm to file..."
    hp.write_alm(phi_alm_file, phi_lm)
    print "----Done."

    #convert to map
    print "----Converting phi to map..."
    phi_map = hp.alm2map(phi_lm, nside, lmax=lmax)
    print "----Done."
    
    if writeMap:
        print "----Writing phi map to file..."
        hp.write_map(phi_map_file, phi_map)
        print "----Done."

    print "=============kappa2phi end============="
    return

if __name__ == "__main__":
    #todo: file selection
    kap2phi(True);
