###########################################################################
# 
#                                NTT Pipeline version 2.4.0
#
#                                INSTALLATION
###########################################################################

THIS README IS OBSOLETE!!!

NTT is written in python and requires the following package:

- Python 2.5 or Python 2.6 or Python 2.7
   these modules have to be installed:
            - numpy
	    - pyraf
	    - matplotlib
	    - pyfits

- sextractor
- swarp
- wget
- Iraf 
	    
##############################################################################
extract the files from the tarball
> tar -xvf NTT.tar

> cd NTT
> python setup.py install  (--record files.txt) 

##########################################################################
To uninstall a previus version 

- delete the ntt directory in your site-package path
- delete the ntt****.egg-info from the same directory
- delete all the PESSTO executable: 
PESSTO,PESSTOEFOSC2dSPEC,PESSTOFASTSPEC,PESSTOSOFI2dSPEC,PESSTOWISE         
PESSTOEFOSC1dSPEC,PESSTOEFOSCPHOT,PESSTOSOFI1dSPEC,PESSTOSOFIPHOT

or if during installation  you used the option: --record files.txt
you can run the following command in the terminal:
> cat files.txt | xargs sudo rm -rf
