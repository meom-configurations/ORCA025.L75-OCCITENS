# ORCA025.L75-OCCITENS
OCCIPUT ensemble configuration with ORCA025.L75

# Reference code
## NEMO
  rev. 4521    of  http://forge.ipsl.jussieu.fr/nemo/svn/trunk 

  svn co -r 4521 http://forge.ipsl.jussieu.fr/nemo/svn/trunk NEMO_3.5beta

## XIOS
  rev. 565 of http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-1.0

  svn co -r 565 http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-1.0

# Description
 This repository hold the CONFIG directory to be used for reproducing this ensemble run.
 
 In MY_SRC there are the fortran modules differing or added to the standard NEMO reference.
 
 cpp_EORCA12.L75-MJMgd16.fcm shows the used CPP keys for this config, and arch_xxx.fcm were used for compilation on OCCIGEN (CINES) and CURIE( TGCC) super-computers.
 
 In EXP00, namelists and xml files are given.  Notice that most of the xml files are templates which are concatenated at run time to form the final iodef.xml file ( XIOS_1.0, for this configuration). This is because to each member of the ensemble  corresponds a NEMO context (NEMO.001, NEMO.002 etc... ) All the context are identical ( except for member 1 which has some extra inter-member diagnostics output ).


 
# References

Development of a probabilistic ocean modelling system based on NEMO 3.5: application at eddying resolution

Laurent Bessières , Stéphanie Leroux , Jean-Michel Brankart , Jean-Marc Molines , Marie-Pierre Moine ,
Pierre-Antoine Bouttier , Thierry Penduff , Laurent Terray , Bernard Barnier , and Guillaume Sérazin

Geosci. Model Dev., 10, 1091-1106, 2017
www.geosci-model-dev.net/10/1091/2017/
doi:10.5194/gmd-10-1091-2017


