# Assimilation code for using FESOM-REcoM with PDAF

# Directories:
* config/
  Namelist files and job script for Albedo

* fesom/
  Modified routines of FESOM and REcoM
  These files should be copied into the src/ directory of FESOM

* pdaf/
  FESOM-PDAF user code
  These routines should be linked or copied into src/pdaf of FESOM

The file CMakelists.txt in the main dirctory should replace the corresponding
file in the main dirceory of FESOM.


# Notes for compilation
* 
* The PDAF library should be compiled separately. The compilation has to be
   consistent with the compile configuration of FESOM.
   (e.g. same compiler and MPI library)
* The code was tested with PDAF V3.0
* src/CMakelists.txt is configured to look for the PDAF library
   in the directory PDAF/ in the main FESOM directory. On can
   compile PDAF elsewhere and link the PDAf main directory to PDAF/
