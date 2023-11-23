# Rho validation

The files included in this repository were used to generate the plots shown during the Run3 (2023) data-taking period. The RhoAnalyser was used to extract the values of rho. The forests were created using the script validation.py, which was being called by the crab_many_files.py submission script. The files were then processed by the script Rho_new_file.cpp, and the macro validation_multifiles.cpp was used for plotting.

Rho_new_file.cpp accepts two arguments, an input and an output file. I used that script for the validation, but it contains some unneeded sections - I tried to reduce the script - the Reduced_rho.cpp, which on first glance gives same results (at least for run XXX22). 

NB: HLT_HIMinimumBiasHF1AND_v1 is replaced by v2 for Run ...45 onwards.

validation_multifiles.cpp - in the body of the main script, the combined root files from the Rho_new_file/Reduced_rho.cpp are provided (with the names appearing in the histograms input by hand), as well as the destination and name of the output pdfs - one containing all fits to distributions, and the other containing the summary plots.

The resulting plots were used in the slides shown in meetings, found also on https://cernbox.cern.ch/s/W4QRGDBqYo7NED8

