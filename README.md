# modified_lenspix
---very much a WIP---
original code found at: http://cosmologist.info/lenspix/

## What does it do?
  Running lens.py will take the kappa (convergence) map and convert it to phi (gravitational potential) to be used as a lensing map. It will then read the nside and lmax of the primary map to be used in the lensing step. The script will then generate a "specific_params.ini" (created from the "generic_params.ini") based on the inputs provided. Lensing is simulated, and the lensed map is written to the provided filename. This will also yield a gradient phi map, and an unlensed and lensed power spectrum.

Much, if not all of the following assumes being used in a SciNet environment.

## Usage
### Required files:
1. Unlensed TQU alm/map FITS file
2. Kappa map (combination of halo and field)

Run "__python lens.py <kappa_map_in> <primary_in> <phi_alm_out> <lensed_map_out> -np <no. of processes>__", where __<kappa_map_in>__ is the filename of the kappa map you want to lens with, __<primary_in>__ is the filename of the primary alm you want lensed, __<phi_alm_out>__ is the filename of the outputted phi alm, and __<lensed_map_out>__ is the filename for the lensed result. __<no. of processes>__ is the same as the __mpirun -np__ argument, since the code simply calls this command at the end, as a result one does not need to call this script with mpirun, but it must done in a job environment. If -np is not provided it will default to 1. 

An additional __make_maps.sh__ script has been added for lensing of multiple maps, which we used for CIB lensing simulations. If using this, make sure the __-np__ arguments agrees with the script's arguments (the comments at the beginning of the file). 

## Helper scripts
There are several helper scripts included. The main ones would be:
* kappa2phi.py
* theory_sim_compare.py (in pycamb_scripts folder)
* makekappa.py
* makegif.py
* makegif.py
* plotgrad.py
* shell_z_integrand.py
* shellcls.py

__kappa2phi.py__ is the script that does the conversion of the kappa map to phi alm.

__theory_sim_compare.py__ is a script that allows one to compare the simulated lensed result with the theoretical lensed calculations based on CAMB. To use, provide the following in the script:
* primary alm filename
* phi alm filename
* lensed map filename

It will also need the theoretical power spectrum of the same maps used in Lenspix, which are:
* kappa power spectrum
* primary power spectrum

Finally, give a filename for the output from CAMB's convolution (the theoretical lensed power spectrum). 

__makekappa.py__ was used in conjunction with theory_sim_compare.py to take the correct power spectrum data from scalCls.dat or lenspotentialCls.dat and combine them with our own simulated kappa maps.

__makegif.py__ will take a set of CIB plots and combine them into a single animation.

__plotgrad.py__ will create a plot of the phi map with the gradients overplotted.

__shell_z_integrand.py__ takes the cl value of each CIB shell at l=500 and plots them.

__shellcls.py__ will plot the power spectra of the shells in one plot.

## Known issues
* lensing a map with many point sources will yield NaNs

## Todos
* different output lmax not yet implemented