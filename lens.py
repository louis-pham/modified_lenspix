from kappa2phi import kap2phi
import subprocess
import argparse

parser = argparse.ArgumentParser(description='Lens CMB primary.')
parser.add_argument('-ol', '--output_lmax', help='new lmax for the lensed map')
args = parser.parse_args()

generic_params = 'generic_params.ini'
specific_params = 'specific_params.ini'
out_file_root = "test_file_stem"

#input filenames
field_kappa_file = "kappa_maps/8Gpc_n4096_nb23_nt18_kap_field.fits"
halo_kappa_file = "kappa_maps/8Gpc_n4096_nb23_nt18_kap_halo.fits"
unlensed_primary_file = "FromNERSC/ffp10_unlensed_scl_cmb_000_alm.fits"
#output filenames
comb_kappa_file = "kappa_maps/8Gpc_n4096_nb23_nt18_kap_comb.fits"
phi_alm_file = "kappa_maps/comb_z_4.6_n4096_lmax5120_phi_alm.fits"
phi_map_file = "kappa_maps/comb_z_4.6_n4096_lmax5120_phi_map.fits"

nside, lmax = kap2phi(field_kappa_file, halo_kappa_file, unlensed_primary_file, comb_kappa_file, phi_alm_file, writeMap=True, phi_map_file=phi_map_file)

print 'Obtained parameters NSIDE:', nside, 'and LMAX:', lmax
print 'Creating params file...'

if args.output_lmax:
    output_lmax = args.output_lmax
    print 'Output lmax set to', output_lmax + '.'
else:
    output_lmax = lmax
    print 'Output lmax defaulted to input lmax', output_lmax + '.'

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
