# modified_lenspix
---very much a WIP---
original code found at: http://cosmologist.info/lenspix/

kappa2phi.py - reads in kappa map and converts to phi map

showmaps.py - displays unlensed/lensed maps, and optionally kappa,phi,diff maps

Simlens.f90 - modified to read in external unlensed map and a phi map

params.ini - nside/lmax should match the nside/lmax of the input primary and phi map

## Usage
### Input files:
1. Unlensed TQU alm
2. Kappa alm

* in __kappa2phi.py__, set the input and output filenames, then run it to get a phi alm file
* in __params.ini__, set the _nside_ and _lmax_ to match that of your two input files, and then set _primary_file_ and _phi_file_ to the filenames of the CMB alm and phi map. Also set the _out_file_root_ to your liking, this will be the file stem of all the output from __simlens__
* run the __simlens__ program -- the outputted files will be in the same directory
