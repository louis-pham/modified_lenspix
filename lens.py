from kappa2phi import kap2phi
import subprocess
import argparse

parser = argparse.ArgumentParser(description='Lens CMB primary.')
parser.add_argument('-ol', '--output_lmax', help='new lmax for the lensed map')
args = parser.parse_args()

generic_params = 'generic_params.ini'
specific_params = 'specific_params.ini'
out_file_root = "test_cib_lens"

#input filenames
field_kappa_file = "kappa_maps/8Gpc_n2048_nb23_nt18_kap_field.fits"
halo_kappa_file = "kappa_maps/8Gpc_n2048_nb23_nt18_kap_halo.fits"
unlensed_primary_file = "kappa_maps/cib_fullsky_ns2048_zmin0.0_zmax1.245_nu217_13579_normalized_alm.fits"

#phi output filenames
comb_kappa_file = "kappa_maps/8Gpc_n2048_nb23_nt18_kap_comb.fits"
phi_alm_file = "kappa_maps/" + out_file_root + "_n2048_lmax5120_phi_alm.fits"
phi_map_file = "kappa_maps/" + out_file_root + "_n2048_lmax5120_phi_map.fits"

nside, lmax = kap2phi(field_kappa_file, halo_kappa_file, unlensed_primary_file, comb_kappa_file, phi_alm_file, writeMap=True, phi_map_file=phi_map_file)

print 'Obtained parameters NSIDE:', nside, 'and LMAX:', lmax

#todo: use nside and lmax from above
#comb_kappa_file = "kappa_maps/8Gpc_" + nside + "_nb23_nt18_kap_comb.fits"
#phi_alm_file = "kappa_maps/" + out_file_root + "_" + nside + "_lmax" + lmax + "_phi_alm.fits"
#phi_map_file = "kappa_maps/" + out_file_root + "_" + nside + "_lmax" + lmax + "_phi_map.fits"

print 'Creating params file...'

if args.output_lmax:
    output_lmax = args.output_lmax
    print 'Output lmax set to', str(output_lmax) + '.'
else:
    output_lmax = lmax
    print 'Output lmax defaulted to input lmax', str(output_lmax) + '.'

subprocess.call(['cp', generic_params, specific_params])
#commas are used as delimiters for sed commands to avoid escaping (back)slashes -- this assumes there's no commas in the filenames
subprocess.call(['sed', '-i', 's,__NSIDEREPLACE__,' + str(nside) + ',g', specific_params])
subprocess.call(['sed', '-i', 's,__LMAXREPLACE__,' + str(lmax) + ',g', specific_params])
subprocess.call(['sed', '-i', 's,__OUTFILEROOTREPLACE__,' + out_file_root + ',g', specific_params])
subprocess.call(['sed', '-i', 's,__OUTPUTLMAXREPLACE__,' + str(output_lmax) + ',g', specific_params])
subprocess.call(['sed', '-i', 's,__PRIMARYFILEREPLACE__,' + unlensed_primary_file + ',g', specific_params])
subprocess.call(['sed', '-i', 's,__PHIFILEREPLACE__,' + phi_alm_file + ',g', specific_params])
print 'Params file created.'
print 'Running simlens...'
subprocess.call(['mpirun', '-np', '24', './simlens', specific_params])
print 'Simlens complete.'
