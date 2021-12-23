#!/usr/bin/env python

import os,sys,string,re,glob
import ntt
from pyfits import open as popen
from ntt.util import readkey3, readhdr, readspectrum, delete, correctcard
import datetime
import time

keyword={}
keyword['efosc']={}
keyword['sofi']={}
keyword['efosc']['image']={'ABMAGLIM':'R','ABMAGSAT':'R','PSF_FWHM':'R','ELLIPTIC':'R',\
                             'CRVAL1':'R','CRPIX1':'R','CTYPE1':'S','CUNIT1':'S','CDELT1':'R','CSYER1':'R','CRDER1':'R',\
                             'CD1_1':'R','CD2_2':'R','CD1_2':'R','CD2_1':'R','CRVAL2':'R','CRPIX2':'R','CTYPE2':'S',\
                             'CUNIT2':'S','CDELT2':'R','CSYER2':'R','CRDER2':'R',\
                               'PHOTZP':'R','PHOTSYS':'S'}#,'PHOTZPER':'R'

keyword['sofi']['image']={'ABMAGLIM':'R','ABMAGSAT':'R','PSF_FWHM':'R','ELLIPTIC':'R',\
                             'CRVAL1':'R','CRPIX1':'R','CTYPE1':'S','CUNIT1':'S','CDELT1':'R','CSYER1':'R','CRDER1':'R',\
                             'CD1_1':'R','CD2_2':'R','CD1_2':'R','CD2_1':'R','CRVAL2':'R','CRPIX2':'R','CTYPE2':'S',\
                             'CUNIT2':'S','CDELT2':'R','CSYER2':'R','CRDER2':'R',\
                              'PHOTZP':'R','PHOTSYS':'S'}#,'PHOTZPER':'R'
keyword['efosc']['spectrum']={'SPECSYS':'S','EXT_OBJ':'L','CONTNORM':'L','TOT_FLUX':'L','FLUXERR':'R','DISPELEM':'S',\
                                 'WAVELMIN':'R','WAVELMAX':'R','LAMRMS':'R','LAMNLIN':'R','SPEC_BIN':'R','SPEC_ERR':'R',\
                                 'SPEC_SYE':'R','SPEC_RES':'R'}
keyword['sofi']['spectrum']={'SPECSYS':'S','EXT_OBJ':'L','CONTNORM':'L','TOT_FLUX':'L','FLUXERR':'R','DISPELEM':'S',\
                                 'WAVELMIN':'R','WAVELMAX':'R','LAMRMS':'R','LAMNLIN':'R','SPEC_BIN':'R','SPEC_ERR':'R',\
                                 'SPEC_SYE':'R','SPEC_RES':'R'}
keyword['sofi']['all']={'ORIGIN':'S','TELESCOP':'S','INSTRUME':'S','FILTER':'S','OBJECT':'S','RA':'R','DEC':'R','EQUINOX':'R',\
                    'RADECSYS':'S','EXPTIME':'R','TEXPTIME':'R','MJD-OBS':'R','MJD-END':'R','PROG_ID':'S','OBID':'I',\
                    'M_EPOCH':'L','SINGLEXP':'L','NCOMBINE':'I','OBSTECH':'S','PROCSOFT':'S','DATAMIN':'R','DATAMAX':'R',\
                    'GAIN':'R','PRODCATG':'S','BUNIT':'S','REFERENC':'S','FLUXCAL':'S',\
                            'DIT':'R','NDIT':'I','NJITTER':'I','NOFFSETS':'I','NUSTEP':'I'}
keyword['efosc']['all']={'ORIGIN':'S','TELESCOP':'S','INSTRUME':'S','FILTER':'S','OBJECT':'S','RA':'R','DEC':'R','EQUINOX':'R',\
                    'RADECSYS':'S','EXPTIME':'R','TEXPTIME':'R','MJD-OBS':'R','MJD-END':'R','PROG_ID':'S','OBID':'I',\
                    'M_EPOCH':'L','SINGLEXP':'L','NCOMBINE':'I','OBSTECH':'S','PROCSOFT':'S','DATAMIN':'R','DATAMAX':'R',\
                    'GAIN':'R','PRODCATG':'S','BUNIT':'S','REFERENC':'S','FLUXCAL':'S'}

keyword['efosc']['1D']={'VOCLASS':'S','VOPUB':'S','TITLE':'S','APERTURE':'R','TELAPSE':'R','TMID':'R','SPEC_VAL':'R','SPEC_BW':'R','SNR':'R'}
keyword['sofi']['1D']={'VOCLASS':'S','VOPUB':'S','TITLE':'S','APERTURE':'R','TELAPSE':'R','TMID':'R','SPEC_VAL':'R','SPEC_BW':'R','SNR':'R'}
keyword['efosc']['2D']={}
keyword['sofi']['2D']={}


typeobj={'R':float,'S':str,'I':int,'L':bool}

from optparse import OptionParser
description="> check header "
usage= "%prog \t listframes [option] "

if __name__ == "__main__":
    parser = OptionParser(usage=usage,description=description,version="%prog 1.0")
    parser.add_option("-v", "--verbose",dest="verbose",action="store_true")
    option,args = parser.parse_args()
    if len(args)<1:
        sys.argv.append('--help')
    option,args = parser.parse_args()
    lista=ntt.util.readlist(args[0])
    lista=ntt.sofiphotredudef.sortbyJD(lista)
    for img in lista:
        correctcard(img)
        hdr=readhdr(img)
        _instrume=readkey3(hdr,'instrume')
        _type=readkey3(hdr,'tech').lower()
        print(img,_type,_instrume)
        for key in keyword[_instrume]['all'].keys():
            if readkey3(hdr,key)==None:
                print(key,readkey3(hdr,key))
            else:
                if type(readkey3(hdr,key)) is not typeobj[keyword[_instrume]['all'][key]]:
                    print(key,readkey3(hdr,key))
                else:
                    if option.verbose:
                        print(key,readkey3(hdr,key))
        for key in keyword[_instrume][_type].keys():
            if readkey3(hdr,key)==None:
                print(key,readkey3(hdr,key))
            else:
                if type(readkey3(hdr,key)) is not typeobj[keyword[_instrume][_type][key]]:
                    print(key,readkey3(hdr,key))
                else:
                    if option.verbose:
                        print(key,readkey3(hdr,key))
        naxis=readkey3(hdr,'naxis')
        if naxis==2: aa='2D'
        else:       aa='1D'
        print(aa)
        for key in keyword[_instrume][aa].keys():
            if readkey3(hdr,key)==None:
                print(key,readkey3(hdr,key))
            else:
                if type(readkey3(hdr,key)) is not typeobj[keyword[_instrume][aa][key]]:
                    print(key,readkey3(hdr,key))
                else:
                    if option.verbose:
                        print (key,readkey3(hdr,key))
