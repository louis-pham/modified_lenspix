# modified_lenspix
---very much a WIP---
original code found at: http://cosmologist.info/lenspix/

kappa2phi.py - reads in kappa map and converts to phi map

showmaps.py - displays unlensed/lensed maps, and optionally kappa,phi,diff maps

Simlens.f90 - modified to read in external unlensed map and a phi map

params.ini - nside/lmax should match the nside/lmax of the input primary and phi map

## Usage
### Required files:
1. Unlensed TQU alm
2. Kappa alms (one halo and one field alm)

In __lens.py__, set the input and output filenames, the output file root, and desired parameter filenames, then run it. The user can also run the script with the '-ol' flag and provide a different lmax for the lensed map (not yet reflected in output).

## What does it do?
Running the lens.py script will take the kappa maps and convert them to phi maps. It will then read the nside and lmax of the primary CMB map to be used in the lensing step. The script will then generate a "specific_params.ini" (created from the "generic_params.ini") that contains the filenames and other parameters for lensing. 
Lensing is executed and saved in files with the _output_file_root_ specified in the script.
