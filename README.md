# ROMS-ICEPLUME Offline Package
Offline version of the ROMS-ICEPLUME module.

This repository contains the Buoyant Plume Thoery (BPT) computation code, and a few python/notebook scripts to extract required data from ROMS history/average files and easily call the code to compute the plume states. The Fotran source code in **src/** is copied from the compiled code of ROMS-ICEPLUME, plus some modifications to read in the velocity and tracer profiles from a text file.

The directory **python/** contains python scripts. The primary script is in **python/run_iceplume_ROMS_nc.py**, which serves two purposes - (1) Read in the ROMS grid, river, and output files, process them and save one profile to to several text files in **inputs/**; (2) Build the executable, and run the executable repeatedly to compute the plume states and save the outputs.

The outputed variables include:

* **z_r**, **z_w**, vertical coordinates; 

* **f**, volume flux of the upwelling plume;

* **w**, vertical velocity of the upwelling plume;

* **t**, **s**, **rho**, active tracer concentration;

* **a**, integrated area of the plume/ice contact area;

* **mInt**, integrated submarine melt rate;

* **ent**, **det**, entrainment and detrainment rate at each layer;

* **detI**, identifier of plume detrainment layers;

* **tAm**, **sAm**, **rhoAm**, ambient active tracer concentration;

* **m**, submarine melt rate of each layer;

* **dye**, passive tracer concentration of the plume.

For more information please refer to my PhD thesis.

I wrote this package because at this moment ROMS-ICEPLUME cannot properly write the plume states to the standard output files. This functionality is still being developed, and will likely be released in a future version.
