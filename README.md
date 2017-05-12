# modified_lenspix
---very much a WIP---
original code found at: http://cosmologist.info/lenspix/

kappa2phi.py - reads in kappa map and converts to phi map

showmaps.py - displays unlensed/lensed maps, and optionally kappa,phi,diff maps

Simlens.f90 - modified to read in external unlensed map and a phi map

params.ini - nside/lmax should match the nside/lmax of the input primary and phi map

## What does it do?
  Running lens.py will take the kappa (convergence) map and convert it to phi (gravitational potential). It will then read the nside and lmax of the primary map to be used in the lensing step. The script will then generate a "specific_params.ini" (created from the "generic_params.ini") based on the inputs provided. Lensing is simulated, and the result is with the provided filename. This will also yield an unlensed and lensed power spectrum.

## Usage
### Required files:
1. Unlensed TQU alm/map FITS file
2. Kappa map (combination of halo and field)

Run __python lens.py <kappa_map_in> <primary_in> <phi_alm_out> <lensed_map_out>__, where __<kappa_map_in>__ is the filename of the kappa map you want to lens with, __<primary_in>__ is the filename of the primary you want lensed, __<phi_alm_out>__ is the filename of the outputted phi alm, and __<lensed_map_out>__ is the filename for the lensed result.

Also, in pycamb_scripts there is the __power_compare.py__ script that allows one to compare the simulated lensed result with the theoretical lensed calculations based on CAMB. __** This is not functioning properly yet **__