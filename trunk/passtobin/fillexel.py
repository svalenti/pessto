#!/usr/bin/env python

import os,sys,string,re,glob
import ntt
from pyfits import open as popen
from ntt.util import readkey3, readhdr, readspectrum, delete
import datetime
import time

from optparse import OptionParser
description="> filling exelfile " 
usage= "%prog \t listframes [option] "

if __name__ == "__main__":
    parser = OptionParser(usage=usage,description=description,version="%prog 1.0")
    parser.add_option("-o", "--output",dest="output",default='',type="str", help='name -l list \t [%default]')
    option,args = parser.parse_args()
    if len(args)<1: # or  option.input=='listfiles': 
        sys.argv.append('--help')
    option,args = parser.parse_args()

    def julday(year,month,day,hour,min):
          import math
          univTime = float(hour)+(float(min)/60.)
          if ((100*float(year))+float(month)-190002.5) >= 0:
               sign = 1
          else:
               sign = -1
          part1 = 367 * float(year)
          part2 = math.floor((7*(float(year)+math.floor((float(month)+9)/12)))/4)
          part3 = float(day)+math.floor((275*float(month))/9)
          part4 = 1721013.5+(univTime/24)
          part5 = 0.5*sign
          jd = part1-part2+part3+part4-part5+0.5
          return jd

    parent_dir = os.getcwd()+'/'
    now=datetime.datetime.now()
    JD=julday(now.year,now.month,now.day,now.hour,now.minute)
    user=os.environ['USER']
    Y,M,D,H,m,s,x,x,x=time.gmtime()
    lista=ntt.util.readlist(args[0])
    if option.output:  _output=option.output
    else:     _output=user+'_'+str(JD)+'_exel.txt'

    lista=ntt.sofiphotredudef.sortbyJD(lista)
    for img in lista:
        hdr=readhdr(img)
        print readkey3(hdr,'ORIGFILE')
        print readkey3(hdr,'ARCFILE')
        print readkey3(hdr,'OBJECT')
        print readkey3(hdr,'date-obs')
        OBID=readkey3(hdr,'HIERARCH ESO OBS ID')
        _object=readkey3(hdr,'object')
        _start=readkey3(hdr,'HIERARCH ESO OBS START')
        _end=str(Y)+'-'+str(M)+'-'+str(D)+'T'+str(H)+':'+str(m)+':'+str(s)
        print 'A \t  OB done fully within constraints \nB \t OB mostly with in constrain\nC \t not usefull for science\n'
        _QC=raw_input('QC Grade (A,B,C) [A] ?')
        if not _QC:_QC='A'
        print 'seeing,clouds,wind (1.0,clear,low)'
        _observcond=raw_input('observing conditions [1.5,clear,10km/s]?')
        if not _observcond: _observcond='1.5,clear,10km/s'
        _RA=readkey3(hdr,'RA')
        _DEC=readkey3(hdr,'DEC')
        _origfile=readkey3(hdr,'ARCFILE')
        _comments=raw_input('comments ? ')
        line=str(OBID)+'\t'+str(_object)+'\t'+str(_start)+'\t'+str(_end)+'\t'+str(_QC)+'\t'+str(_observcond)+\
            '\t'+str(_RA)+'\t'+str(_DEC)+'\t'+str(_origfile)+'\t'+str(_comments)+'\n'
        print line
        if os.path.isfile(_output): f=open(_output,'a')
        else: f=open(_output,'w')
        f.write(line)
        f.close()
        print '\n info stored in  '+str(_output)
