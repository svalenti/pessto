try:      from astropy.io import fits as pyfits
except:   import pyfits

def ReadAscii2(ascifile):
    # print "LOGX:: Entering `ReadAscii2` method/function in %(__file__)s" %
    # globals()
    import string

    f = open(ascifile, 'r')
    ss = f.readlines()
    f.close()
    vec1, vec2 = [], []
    for line in ss:
        if line[0] != '#':
            vec1.append(float(string.split(line)[0]))
            vec2.append(float(string.split(line)[1]))
    return vec1, vec2


# ########################################################################
def readlist(listfile):

    from ntt.util import correctcard
    import string
    import sys
    import re
    import glob

    if '*' in listfile:
        imglist = glob.glob(listfile)
    elif ',' in listfile:
        imglist = string.split(listfile, sep=',')
    else:
        try:
            hdulist = pyfits.open(listfile)
        except:
            hdulist = []
        if hdulist:
            imglist = [listfile]
        else:
            try:
                ff = open(listfile, 'r')
                files = ff.readlines()
                ff.close()
                imglist = []
                for ff in files:
                    ff = re.sub(' ', '', ff)
                    if not ff == '\n' and ff[0] != '#':
                        ff = re.sub('\n', '', ff)
                        try:
                            hdulist = pyfits.open(ff)
                            imglist.append(ff)
                        except:
                            try:
                                correctcard(ff)
                                hdulist = pyfits.open(ff)
                                imglist.append(ff)
                            except:
                                pass
            except:
                sys.exit('\n##### Error ###\n file ' +
                         str(listfile) + ' do not  exist\n')
    if len(imglist) == 0:
        sys.exit('\n##### Error ###\nIf "' + str(listfile)
                 + '" is an image, it is corrupted \n or is not a list of image\n')
    return imglist


##############################################################################
def delete(listfile):

    import os
    import string
    import re
    import glob

    if listfile[0] == '@':
        ff = open(listfile[1:])
        files = ff.readlines()
        imglist = []
        for ff in files:
            ff = re.sub(' ', '', ff)
            if not ff == '\n' and ff[0] != '#':
                ff = re.sub('\n', '', ff)
                imglist.append(ff)
    elif ',' in listfile:
        imglist = string.split(listfile, sep=',')
    else:
        imglist = [listfile]
    lista = []
    for _file in imglist:
        lista = lista + glob.glob(_file)
    if lista:
        for _file in lista:
            try:
                os.system('rm ' + _file)
            except:
                pass


###############################################################
def readhdr(img):
    try:
        hdr = pyfits.open(img)[0].header
    except:
        import ntt
        try:
            ntt.util.correctcard(img)
        except:
            import sys

            sys.exit('image ' + str(img) +
                     ' is corrupted, delete it and start again')
        hdr = pyfits.open(img)[0].header
    return hdr


def readkey3(hdr, keyword):
    import sys
    import re
    import string
    try:
        if int( re.sub('\.', '', str(pyfits.__version__))[:2] ) <= 30:
            aa = 'HIERARCH '
        else:
            aa = ''
    except:
        aa = ''

    try:
        _instrume = hdr.get('INSTRUME').lower()
    except:
        _instrume = 'none'
    if _instrume == 'efosc':
        useful_keys = {'object': 'OBJECT',
                       'date-obs': 'DATE-OBS',
                       'ut': 'DATE-OBS',
                       'RA': 'RA',
                       'DEC': 'DEC',
                       'datamin': -100,
                       'datamax': 60000,
                       'observer': 'OBSERVER',
                       'exptime': 'EXPTIME',
                       'instrume': 'INSTRUME',
                       'JD': 'MJD-OBS',
                       'lamp': 'LMP_ID',
                       'esoprog': aa + 'ESO OBS PROG ID',
                       'filter': aa + 'ESO INS FILT1 NAME',
                       'grism': aa + 'ESO INS GRIS1 NAME',
                       'catg': aa + 'ESO DPR CATG',
                       'tech': aa + 'ESO DPR TECH',
                       'type': aa + 'ESO DPR TYPE',
                       'gain': aa + 'ESO DET OUT1 GAIN',
                       'ron': aa + 'ESO DET OUT1 RON',
                       'esoid': aa + 'ESO OBS ID',
                       'binx': aa + 'ESO DET WIN1 BINX',
                       'speed': aa + 'ESO DET READ SPEED',
                       'posang': aa + 'ESO ADA POSANG',
                       'airmass': aa + 'ESO TEL AIRM START',
                       'airmass1': aa + 'ESO TEL AIRM END',
                       'slit': aa + 'ESO INS SLIT1 NAME',
                       'obsmode': aa + 'ESO DPR CATG',
                       'telescop': 'TELESCOP'}
    elif _instrume == 'sofi':
        useful_keys = {'object': 'OBJECT',
                       'date-obs': 'DATE-OBS',
                       'ut': 'DATE-OBS',
                       'RA': 'RA',
                       'DEC': 'DEC',
                       'exptime': 'EXPTIME',
                       'observer': 'OBSERVER',
                       'gain': 5.4,
                       'ron': 2.1,
                       'instrume': 'INSTRUME',
                       'JD': 'MJD-OBS',
                       'lamp': 'LMP_ID',
                       'posang': aa + 'ESO ADA POSANG',
                       'esoprog': aa + 'ESO OBS PROG ID',
                       'esotplid': aa + 'ESO TPL ID',
                       'filter': aa + 'ESO INS FILT1 NAME',
                       'grism': aa + 'ESO INS OPTI2 NAME',
                       'catg': aa + 'ESO DPR CATG',
                       'tech': aa + 'ESO DPR TECH',
                       'type': aa + 'ESO DPR TYPE',
                       'dit': aa + 'ESO DET DIT',
                       'ndit': aa + 'ESO DET NDIT',
                       'nexp': aa + 'ESO TPL NEXP',
                       'esoid': aa + 'ESO OBS ID',
                       'airmass': aa + 'ESO TEL AIRM START',
                       'lamp1': aa + 'ESO INS LAMP1 NAME',
                       'lamp3': aa + 'ESO INS LAMP3 NAME',
                       'slit': aa + 'ESO INS OPTI1 NAME',
                       'xcum': aa + 'ESO SEQ CUMOFFSETX',
                       'ycum': aa + 'ESO SEQ CUMOFFSETY',
                       'obsmode': aa + 'ESO DPR CATG',
                       'telescop': 'TELESCOP'}
    else:
        useful_keys = {'object': 'OBJECT',
                       'date-obs': 'DATE-OBS'}
    if keyword in useful_keys:
        if type(useful_keys[keyword]) == float:
            value = useful_keys[keyword]
        else:
            value = hdr.get(useful_keys[keyword])
            if keyword == 'date-obs':
                import string
                import re

                try:
                    value = re.sub('-', '', string.split(value, 'T')[0])
                except:
                    pass
            elif keyword == 'ut':
                import string
                import re

                try:
                    value = string.split(value, 'T')[1]
                except:
                    pass
            elif keyword == 'JD':
                value = value + 0.5
            elif keyword == 'instrume':
                value = value.lower()
        if type(value) == str:
            value = re.sub('\#', '', value)
    else:
        if keyword == 'date-night':
            import datetime

            _date = readkey3(hdr, 'DATE-OBS')
            a = (datetime.datetime.strptime(string.split(_date, '.')[0], "20%y-%m-%dT%H:%M:%S") - datetime.timedelta(
                .5)).isoformat()
            value = re.sub('-', '', string.split(a, 'T')[0])
        else:
            try:
                value = hdr.get(keyword)
            except:
                sys.exit('Warning: keyword not valid')
        #    if type(value) == str:    value=re.sub('\#','',value)
    return value


#######################################################
def writeinthelog(text, logfile):
    # print "LOGX:: Entering `writeinthelog` method/function in %(__file__)s"
    # % globals()
    f = open(logfile, 'a')
    f.write(text)
    f.close()


################################################
def correctcard(img):
    import numpy as np
    # print "LOGX:: Entering `correctcard` method/function in %(__file__)s" %
    # globals()
    from numpy import asarray
    import re
    import os
    hdulist = pyfits.open(img)
    a = hdulist[0]._verify('exception')
    _header = hdulist[0].header
    hdulist.close()

    ######   change 20161003 
    #print a
    #for i in range(len(a)):
    #    if not a[i]:
    #        a[i] = ['']
    #ww = asarray([i for i in range(len(a)) if (re.sub(' ', '', a[i][0]) != '')])

    ww = np.where(a)[0]
    try:
        if int(re.sub('\.', '', str(pyfits.__version__))[:2]) <= 30:
            aa = 'HIERARCH '
        else:
            aa = ''
    except:
        aa = ''

    print aa
    print ww
    if len(ww) > 0:
        newheader = []
        headername = []
        for j in _header.items():
            headername.append(j[0])
            newheader.append(j[1])
        for i in ww:
            if len(headername[i]) > 8:
                data, hdr = pyfits.getdata(img, 0, header=True)
                comm = hdr.comments[aa + headername[i]]
                hdr.pop(aa + headername[i])
                cc = re.sub('HIERARCH ', '', headername[i])
                hdr.append(
                    ('HIERARCH ' + cc, newheader[i], comm), False, False)
                os.system('rm ' + img)
                pyfits.writeto(img, data, hdr)
            else:
                try:
                    imm = pyfits.open(img, mode='update')
                    imm.close(output_verify='silentfix', verbose=False)
                except:
                    imm = pyfits.open(img, mode='update')
                    _header = imm[0].header
                    for i in ww:
                        if headername[i]:
                            if len(headername[i]) > 8 and 'HIERARCH' not in headername[i]:
                                bb = 'HIERARCH '
                            else:
                                bb = ''
                            try:
                                _header.update(
                                    bb + headername[i], newheader[i])
                            except:
                                _header.update(bb + headername[i], 'xxxx')
                    imm.flush()
                    imm.close()
    else:
        pass


#    data, hdr = pyfits.getdata(img, 0, header=True)
#       print 'correction not needed'
##########################################################################

def updateheader(image, dimension, headerdict):
    # added to cut long header
    while len(max([str(headerdict[i][0]) for i in headerdict], key=len)) > 68:
        key = [i for i in headerdict]
        valori, commenti = zip(*[headerdict[i] for i in headerdict])
        num = valori.index(max([str(headerdict[i][0])
                                for i in headerdict], key=len))
        headerdict[key[num]] = [valori[num][0:68], commenti[num]]
        print 'warning: header to long, ', str(key[num]), str(valori[num][0:68]), str(commenti[num])
        #   keytochange=hdr.keys()[hdr.values().index(max([str(i) for i in hdr.values()],key=len))]
        #   hdr[keytochange]=[str(hdr[keytochange])[0:68]]

    if 'version' in dir(pyfits):
        try:
            imm = pyfits.open(image, mode='update')
            _header = imm[dimension].header
            for i in headerdict.keys():
                _header.update(i, headerdict[i][0], headerdict[i][1])
            imm.flush()
            imm.close()
        except:
            from ntt.util import correctcard
            print 'warning: problem to update header, try to correct header format ....'
            correctcard(image)
            try:
                print headerdict
                imm = pyfits.open(image, mode='update')
                _header = imm[dimension].header
                for i in headerdict.keys():
                    _header.update(i, headerdict[i][0], headerdict[i][1])
                imm.flush()
                imm.close()
            except:
                print 'error: not possible update header'
    else:
        #
        #
        # astropy.io.fits requre a tuple to update header
        #
        #
        imm = pyfits.open(image, mode='update')
        _header = imm[dimension].header
        for i in headerdict.keys():
            _header.update( { i : (headerdict[i][0], headerdict[i][1]) } )
        imm.flush()
        imm.close()
        
##########################################################################


def display_image(img, frame, _z1, _z2, scale, _xcen=0.5, _ycen=0.5, _xsize=1, _ysize=1, _erase='yes'):
    # print "LOGX:: Entering `display_image` method/function in %(__file__)s"
    # % globals()
    goon = 'True'
    import glob
    import subprocess
    import time
    import os

    ds9 = subprocess.Popen("ps -U" + str(os.getuid()) + "|grep -v grep | grep ds9", shell=True,
                           stdout=subprocess.PIPE).stdout.readlines()
    if len(ds9) == 0:
        subproc = subprocess.Popen('ds9', shell=True)
        time.sleep(3)

    if glob.glob(img):
        from pyraf import iraf

        iraf.images(_doprint=0)
        iraf.tv(_doprint=0)
        import string
        import os

        if _z2:
            try:
                sss = iraf.display(img, frame, xcen=_xcen, ycen=_ycen, xsize=_xsize, ysize=_ysize, erase=_erase,
                                   fill='yes', zscale='no', zrange='no', z1=_z1, z2=_z2, Stdout=1)
            except:
                print ''
                print '### ERROR: PROBLEM OPENING DS9'
                print ''
                goon = 'False'
        else:
            try:
                sss = iraf.display(img, frame, xcen=_xcen, ycen=_ycen, xsize=_xsize, ysize=_ysize, erase=_erase,
                                   fill='yes', Stdout=1)
            except:
                print ''
                print '### ERROR: PROBLEM OPENING DS9'
                print ''
                goon = False

        if scale and goon:
            answ0 = raw_input('>>> Cuts OK ? [y/n] ? [y] ')
            if not answ0:
                answ0 = 'y'
            elif answ0 == 'no' or answ0 == 'NO':
                answ0 = 'n'

            while answ0 == 'n':
                _z11 = float(string.split(string.split(sss[0])[0], '=')[1])
                _z22 = float(string.split(string.split(sss[0])[1], '=')[1])
                z11 = raw_input('>>> z1 = ? [' + str(_z11) + '] ? ')
                z22 = raw_input('>>> z2 = ? [' + str(_z22) + '] ? ')
                if not z11:
                    z11 = _z11
                else:
                    z11 = float(z11)
                if not z22:
                    z22 = _z22
                else:
                    z22 = float(z22)
                print z11, z22
                sss = iraf.display(img, frame, fill='yes', xcen=_xcen, ycen=_ycen, xsize=_xsize, ysize=_ysize,
                                   erase=_erase, zrange='no', zscale='no', z1=z11, z2=z22, Stdout=1)
                answ0 = raw_input('>>> Cuts OK ? [y/n] ? [y] ')
                if not answ0:
                    answ0 = 'y'
                elif answ0 == 'no' or answ0 == 'NO':
                    answ0 = 'n'
        if goon:
            _z1, _z2 = string.split(string.split(sss[0])[0], '=')[
                1], string.split(string.split(sss[0])[1], '=')[1]
    else:
        print 'Warning: image ' + str(img) + ' not found in the directory '
    return _z1, _z2, goon


###########################################################################
def searcharc(img, listarc):
    # print "LOGX:: Entering `searcharc` method/function in %(__file__)s" %
    # globals()
    import ntt
    import glob
    import numpy as np

    hdr = ntt.util.readhdr(img)
    JD = ntt.util.readkey3(hdr, 'JD')
    _instrume = ntt.util.readkey3(hdr, 'instrume')
    grism0 = ntt.util.readkey3(hdr, 'grism')
    filter0 = ntt.util.readkey3(hdr, 'filter')
    slit0 = ntt.util.readkey3(hdr, 'slit')
    if slit0 == 'slit5.0':
        slit0 = 'slit1.0'
    if not listarc:
        directory = ntt.__path__[
            0] + '/archive/' + str(_instrume) + '/arc/' + grism0 + '/' + filter0 + '/' + slit0
        listarc = glob.glob(directory + '/*fits')
    else:
        directory = ''
    if listarc:
        arcfile = ''
        distance = []
        goodlist = []
        for arc in listarc:
            hdra = ntt.util.readhdr(arc)
            JDarc = ntt.util.readkey3(hdra, 'JD')
            filter1 = ntt.util.readkey3(hdra, 'filter')
            grism1 = ntt.util.readkey3(hdra, 'grism')
            slit1 = ntt.util.readkey3(hdra, 'slit')
            if filter0 == filter1 and slit0 == slit1 and grism0 == grism1:
                goodlist.append(arc)
                distance.append(np.abs(JD - JDarc))
        if len(distance) >= 1:
            arcfile = goodlist[np.argmin(distance)]
        else:
            arcfile = ''
    else:
        arcfile = ''
    return arcfile, directory


###########################################################################
def searchsens(img, listsens):
    # print "LOGX:: Entering `searchsens` method/function in %(__file__)s" %
    # globals()
    import ntt
    import glob
    import numpy as np

    hdr = ntt.util.readhdr(img)
    JD = ntt.util.readkey3(hdr, 'JD')
    _instrume = ntt.util.readkey3(hdr, 'instrume')
    grism0 = ntt.util.readkey3(hdr, 'grism')
    filter0 = ntt.util.readkey3(hdr, 'filter')
    if not listsens:
        directory = ntt.__path__[0] + '/archive/' + \
            str(_instrume) + '/sens/' + grism0 + '/' + filter0
        listsens = glob.glob(directory + '/*fits')
    else:
        directory = ''
    if listsens:
        sensfile = ''
        distance = []
        goodlist = []
        for sens in listsens:
            hdrs = ntt.util.readhdr(sens)
            JDsens = ntt.util.readkey3(hdrs, 'JD')
            filter1 = ntt.util.readkey3(hdrs, 'filter') \
                if ntt.util.readkey3(hdrs, 'filter') else ntt.util.readkey3(hdrs, 'FILTER')
            grism1 = ntt.util.readkey3(hdrs, 'grism') \
                if ntt.util.readkey3(hdrs, 'grism') else ntt.util.readkey3(hdrs, 'GRISM')
            if filter0 == filter1 and grism0 == grism1:
                goodlist.append(sens)
                distance.append(np.abs(JD - JDsens))
        if len(distance) >= 1:
            sensfile = goodlist[np.argmin(distance)]
        else:
            sensfile = ''
    else:
        sensfile = ''
    return sensfile, directory


###########################################################################
def searchflat(img, listflat):
    # print "LOGX:: Entering `searchflat` method/function in %(__file__)s" %
    # globals()
    import ntt
    import glob
    import numpy as np

    hdr = ntt.util.readhdr(img)
    JD = ntt.util.readkey3(hdr, 'JD')
    _instrume = ntt.util.readkey3(hdr, 'instrume')
    filter0 = ntt.util.readkey3(hdr, 'filter')
    if not listflat:
        directory = ntt.__path__[0] + '/archive/' + \
            str(_instrume) + '/flat/' + filter0
        listflat = glob.glob(directory + '/*fits')
    else:
        directory = ''
    if listflat:
        faltfile = ''
        distance = []
        goodlist = []
        for flat in listflat:
            hdrf = ntt.util.readhdr(flat)
            JDflat = ntt.util.readkey3(hdrf, 'JD')
            filter1 = ntt.util.readkey3(hdrf, 'filter')
            if filter0 == filter1:
                goodlist.append(flat)
                distance.append(np.abs(JD - JDflat))
        if len(distance) >= 1:
            flatfile = goodlist[np.argmin(distance)]
        else:
            flatfile = ''
    else:
        flatfile = ''
    return flatfile, directory


###########################################################################
def choseclosest(img0, listimg):
    # print "LOGX:: Entering `choseclosest` method/function in %(__file__)s" %
    # globals()
    import ntt
    import numpy as np

    hdr0 = ntt.util.readhdr(img0)
    JD0 = ntt.util.readkey3(hdr0, 'JD')
    if listimg:
        distance = []
        for img in listimg:
            hdrf = ntt.util.readhdr(img)
            JDimg = ntt.util.readkey3(hdrf, 'JD')
            distance.append(np.abs(JD0 - JDimg))
        if len(distance) >= 1:
            image = listimg[np.argmin(distance)]
        else:
            image = ''
    else:
        image = ''
    return image

###########################################################################


def readstandard(standardfile):
    # print "LOGX:: Entering `readstandard` method/function in %(__file__)s" %
    # globals()
    import ntt
    import numpy as np
    import string
    import os

    if os.path.isfile(standardfile):
        listastandard = standardfile
    elif standardfile[0] == '/':
        listastandard = standardfile
    else:
        listastandard = ntt.__path__[0] + '/standard/stdlist/' + standardfile
    f = open(listastandard, 'r')
    liststd = f.readlines()
    f.close()
    star, ra, dec = [], [], []
    magnitude = []
    for i in liststd:
        if i[0] != '#':
            star.append(string.split(i)[0])
            _ra = string.split(string.split(i)[1], ':')
            _dec = string.split(string.split(i)[2], ':')
            ra.append(
                (float(_ra[0]) + ((float(_ra[1]) + (float(_ra[2]) / 60.)) / 60.)) * 15)
            if '-' in str(_dec[0]):
                dec.append(
                    (-1) * (np.abs(float(_dec[0])) + ((float(_dec[1]) + (float(_dec[2]) / 60.)) / 60.)))
            else:
                dec.append(
                    float(_dec[0]) + ((float(_dec[1]) + (float(_dec[2]) / 60.)) / 60.))
            try:
                magnitude.append(string.split(i)[3])
            except:
                magnitude.append(999)
    return np.array(star), np.array(ra), np.array(dec), np.array(magnitude)

##########################################################################


def readspectrum(img):
    # print "LOGX:: Entering `readspectrum` method/function in %(__file__)s" %
    # globals()
    from numpy import array
    import string

    fl = ''
    lam = ''
    graf = 1
    spec = pyfits.open(img)
    head = spec[0].header
    try:
        if spec[0].data.ndim == 1:
            fl = spec[0].data
        elif spec[0].data.ndim == 2:
            fl = spec[0].data[:, 0]
        elif spec[0].data.ndim == 3:
            fl = spec[0].data[0, 0, :]
    except:
        if spec[0].data.rank == 1:
            fl = spec[0].data
        elif spec[0].data.rank == 2:
            fl = spec[0].data[:, 0]
        elif spec[0].data.rank == 3:
            fl = spec[0].data[0, 0, :]
    naxis1 = head['naxis1']
    try:
        crpix1 = head['crpix1']
        crval1 = head['crval1']
        try:
            cdelt1 = head['cdelt1']
        except:
            cdelt1 = head['cd1_1']
        pix = array(range(1, naxis1 + 1, 1))
        pix = array(range(1, len(fl) + 1, 1))
        lam = (pix - crpix1) * cdelt1 + crval1
    except:
        try:
            WAT = head['WAT2_001']
            pix = array(range(1, naxis1 + 1, 1))
            crpix1 = string.split(string.split(WAT, '"')[1])[0]
            crval1 = string.split(string.split(WAT, '"')[1])[3]
            cdelt1 = string.split(string.split(WAT, '"')[1])[4]
            lam = (pix - float(crpix1)) * float(cdelt1) + float(crval1)
        except:
            graf = 0
    return lam, fl


###########################################################################
def pval(_xx, p):
    # print "LOGX:: Entering `pval` method/function in %(__file__)s" %
    # globals()
    _y = +p[0] + p[1] * _xx
    return _y


def residual(p, y, x):
    # print "LOGX:: Entering `residual` method/function in %(__file__)s" %
    # globals()
    for i in range(len(p)):
        err = (y - p[i] * x ** i)
    return err


#########################################################################
def defsex(namefile):
    # print "LOGX:: Entering `defsex` method/function in %(__file__)s" %
    # globals()
    import ntt
    import string

    sexfile = ntt.__path__[0] + '/standard/sex/default.sex'
    f = open(sexfile, 'r')
    ss = f.readlines()
    f.close()
    ff = open(namefile, 'w')
    for i in ss:
        if string.count(i, 'PARAMETERS_NAME') == 1:
            ff.write('PARAMETERS_NAME  "' +
                     ntt.__path__[0] + '/standard/sex/default.param"\n')
        elif string.count(i, 'FILTER_NAME') == 1:
            ff.write('FILTER_NAME  "' +
                     ntt.__path__[0] + '/standard/sex/default.conv"\n')
        elif string.count(i, 'STARNNW_NAME') == 1:
            ff.write('STARNNW_NAME "' +
                     ntt.__path__[0] + '/standard/sex/default.nnw"\n')
        else:
            ff.write(i)
    ff.close()
    return namefile


############################################################

def defswarp(namefile, imgname, _combine, gain=''):
    # print "LOGX:: Entering `defswarp` method/function in %(__file__)s" %
    # globals()
    import ntt
    import string
    import re

    if _combine.lower() in ['median']:
        _combine = 'MEDIAN'
    elif _combine.lower() in ['average']:
        _combine = 'AVERAGE'
    elif _combine.lower() in ['sum']:
        _combine = 'SUM'
    swarpfile = ntt.__path__[0] + '/standard/sex/default.swarp'
    f = open(swarpfile, 'r')
    ss = f.readlines()
    f.close()
    ff = open(namefile, 'w')
    for i in ss:
        if string.count(i, 'IMAGEOUT_NAME') == 1:
            ff.write('IMAGEOUT_NAME    ' + str(imgname) +
                     '  # Output filename \n')
        elif string.count(i, 'WEIGHTOUT_NAME') == 1:
            ff.write('WEIGHTOUT_NAME   ' + str(
                re.sub('.fits', '.weight.fits', imgname)) + '  # Output weight-map filename  \n')
        elif string.count(i, 'COMBINE_TYPE') == 1:
            ff.write('COMBINE_TYPE    ' + str(_combine) +
                     '  # MEDIAN,AVERAGE,MIN,MAX,WEIGHTED,CHI2 \n')
        elif string.count(i, 'GAIN_DEFAULT') == 1:
            if gain:
                ff.write('GAIN_DEFAULT    ' + str(gain) +
                         '  # Default gain if no FITS keyword found \n')
            else:
                ff.write(i)
        else:
            ff.write(i)
    ff.close()
    return namefile


##########################################################################
def archivefile(img, overwrite=True):
    # print "LOGX:: Entering `archivefile` method/function in %(__file__)s" %
    # globals()
    import os
    import ntt

    outputfile = ntt.util.readkey3(ntt.util.readhdr(img), 'ARCFILE')
    print outputfile
    if not overwrite and os.path.isfile(outputfile):
        answ = raw_input('overwrite file ' + outputfile + ' [[y]/n]? ')
        if not answ:
            answ = 'y'
    else:
        answ = 'y'
    if answ.lower() in ['yes', 'y']:
        os.system('cp ' + img + ' ' + outputfile)
    else:
        outputfile = ''
    return outputfile


##########################################################################
def airmass(img, overwrite=True, _observatory='lasilla'):
    # print "LOGX:: Entering `airmass` method/function in %(__file__)s" %
    # globals()
    import ntt
    from pyraf import iraf

    iraf.astutil(_doprint=0)
    hdr = ntt.util.readhdr(img)
    if readkey3(hdr, 'UTC'):
        _UT = (ntt.util.readkey3(hdr, 'UTC') +
               (ntt.util.readkey3(hdr, 'exptime') / 2)) / 3600
        _date = ntt.util.readkey3(hdr, 'date-obs')
        _date = _date[0:4] + '-' + _date[4:6] + '-' + _date[6:8]
        _RA = ntt.util.readkey3(hdr, 'RA') / 15
        _DEC = ntt.util.readkey3(hdr, 'DEC')
        f = file('airmass.txt', 'w')
        f.write('mst = mst ("' + str(_date) + '",' + str(_UT) +
                ', obsdb ("' + str(_observatory) + '", "longitude"))\n')
        f.write(
            'air = airmass (' + str(_RA) + ',' + str(_DEC) + ',mst, obsdb ("' + str(_observatory) + '", "latitude"))\n')
        f.write('print(air)\n')
        f.close()
        _air = iraf.astcalc(image=img, command="airmass.txt", Stdout=1)[0]
        try:
            _air = float(_air)
        except:
            _air = 999
        ntt.util.delete('airmass.txt')
        if overwrite and _air < 99.:
            ntt.util.updateheader(
                img, 0, {'AIRMASS': [_air, 'mean airmass computed with astcalc']})
    else:
        _air = ''
    return _air


##########################################################################
def dvex():
    # print "LOGX:: Entering `dvex` method/function in %(__file__)s" %
    # globals()
    dv = {}
    dv['line'] = {'Gr16': 300, 'Gr11': 430, 'Gr13': 200, 'GR': 150, 'GB': 430, 'Gr18': 430, 'Gr20': 430}
    dv['std'] = {'_t_order': 6, '_t_niter': 50, '_t_sample': '*', '_t_nlost': 20, '_width': 10, '_radius': 10,
                 '_weights': 'variance',
                 '_nsum': 30, '_t_step': 10, '_t_nsum': 10, '_lower': -10, '_upper': 10, '_b_sample': '-40:-20,20:40',
                 '_resize': 'no'}
    dv['obj'] = {'_t_order': 4, '_t_niter': 50, '_t_sample': '*', '_t_nlost': 20, '_width': 10, '_radius': 10,
                 '_weights': 'variance',
                 '_nsum': 40, '_t_step': 10, '_t_nsum': 10, '_lower': -5, '_upper': 5, '_b_sample': '-25:-15,15:25',
                 '_resize': 'yes'}
    return dv


##########################################################################

def phase3header(img):
    # print "LOGX:: Entering `phase3header` method/function in %(__file__)s" %
    # globals()
    import ntt
    import numpy as np

    img_data = pyfits.open(img)[0].data
    hdr = ntt.util.readhdr(img)

    hedvec0 = {'DATAMIN': [float(min(img_data[np.isfinite(img_data)])), 'Minimal pixel value'],
               'DATAMAX': [float(max(img_data[np.isfinite(img_data)])), 'Maximum pixel value'],
               'REFERENC': ['Smartt et al 2014', 'Bibliographic reference'],
               'ORIGIN': ['ESO', 'European Southern Observatory'],
               'PI-COI': ['Smartt', 'PI-COI name'],
               'PROCSOFT': ['ntt_' + str(ntt.__version__), 'pipeline version']}

    if ntt.util.readkey3(hdr, 'filter'):
        hedvec0['FILTER'] = [ntt.util.readkey3(hdr, 'filter'), 'Filter name']
    if ntt.util.readkey3(hdr, 'gain'):
        hedvec0['GAIN'] = [ntt.util.readkey3(
            hdr, 'gain'), 'Conversion from electrons to ADU']
    if ntt.util.readkey3(hdr, 'esoid'):
        hedvec0['OBID1'] = [
            int(str(ntt.util.readkey3(hdr, 'esoid'))), 'Observation block ID']
    if ntt.util.readkey3(hdr, 'esoprog'):
        hedvec0['PROG_ID'] = [ntt.util.readkey3(
            hdr, 'esoprog'), 'ESO program identification']
    if ntt.util.readkey3(hdr, 'tech'):
        hedvec0['OBSTECH'] = [ntt.util.readkey3(
            hdr, 'tech'), 'Observation technique']

    # added for DR2
    # check longest value in the dictionary header
    # hedvec0['LONGSTRN']=[True,'TRUE if header values longer than 68 character']
    # while len(max([str(i) for i in hdr.values()],key=len)) > 68:
    #   keytochange=hdr.keys()[hdr.values().index(max([str(i) for i in hdr.values()],key=len))]
    #   hdr[keytochange]=[str(hdr[keytochange])[0:68]]
    ntt.util.updateheader(img, 0, hedvec0)


##########################################################################
def rangedata(lista):
    # print "LOGX:: Entering `rangedata` method/function in %(__file__)s" %
    # globals()
    import ntt
    import numpy as np

    if len(lista) >= 1:
        JD, date = [], []
        _instrume = ntt.util.readkey3(ntt.util.readhdr(lista[0]), 'instrume')
        for img in lista:
            try:
                hdr = ntt.util.readhdr(img)
                JD.append(ntt.util.readkey3(hdr, 'JD'))
                date.append(ntt.util.readkey3(hdr, 'date-night'))
            except:
                pass
        datemin = date[np.argmin(np.array(JD))]
        month = {'01': 'jan', '02': 'feb', '03': 'mar', '04': 'apr', '05': 'may', '06': 'jun', '07': 'jul', '08': 'aug',
                 '09': 'sep', '10': 'oct', '11': 'nov', '12': 'dec'}
        datemax = date[np.argmax(np.array(JD))]
        stringa = month[datemin[4:6]] + datemin[2:4] + '_d' + str(datemin[-2:]) + 'to' + str(
            datemax[-2:]) + '_' + _instrume
    else:
        stringa = ''
    return stringa


##########################################################################
def name_duplicate(img, nome, ext):
    # print "LOGX:: Entering `name_duplicate` method/function in %(__file__)s"
    # % globals()
    import glob
    import ntt

    dimg = ntt.util.readkey3(ntt.util.readhdr(img), 'DATE-OBS')
    listafile = glob.glob(nome + '_?' + ext + '.fits') + \
        glob.glob(nome + '_??' + ext + '.fits')
    if len(listafile) == 0:
        nome = nome + "_1" + ext + '.fits'
    else:
        date = []
        for l in listafile:
            date.append(ntt.util.readkey3(ntt.util.readhdr(l), 'DATE-OBS'))
        if dimg in date:
            nome = listafile[date.index(dimg)]
        #         if overwrite:
        #            delete(nome)
        else:
            n = 1
            while nome + '_' + str(n) + str(ext) + '.fits' in listafile:
                n = n + 1
            nome = nome + '_' + str(n) + str(ext) + '.fits'
    return nome


###############################################################################
def correctobject(img, coordinatefile):
    # print "LOGX:: Entering `correctobject` method/function in %(__file__)s"
    # % globals()
    import re
    import ntt
    import numpy as np

    scal = np.pi / 180.
    std, rastd, decstd, magstd = ntt.util.readstandard(coordinatefile)
    img = re.sub('\n', '', img)
    ntt.util.correctcard(img)
    hdr = ntt.util.readhdr(img)
    _ra = ntt.util.readkey3(hdr, 'RA')
    _dec = ntt.util.readkey3(hdr, 'DEC')
    dd = np.arccos(np.sin(_dec * scal) * np.sin(decstd * scal) + np.cos(_dec * scal) * np.cos(decstd * scal) *
                   np.cos((_ra - rastd) * scal)) * ((180 / np.pi) * 3600)
    if min(dd) < 200:
        ntt.util.updateheader(
            img, 0, {'OBJECT': [std[np.argmin(dd)], 'Original target.']})
        aa, bb, cc = rastd[np.argmin(dd)], decstd[np.argmin(dd)], std[
            np.argmin(dd)]
    else:
        aa, bb, cc = '', '', ''
    return aa, bb, cc


##########################################################################
def archivingtar(outputlist, rawfile):
    # print "LOGX:: Entering `archivingtar` method/function in %(__file__)s" %
    # globals()
    import os
    import re
    import ntt

    print '\n### making a tar with pre-reduced frames ........ please wait'
    stringa = ''
    for img in outputlist:
        stringa = stringa + img + ' '
    stringa = stringa + rawfile
    ntt.util.delete(re.sub('raw.list', 'tar.gz', rawfile))
    os.system('tar -zcf ' + re.sub('raw.list',
                                   'tar.gz', rawfile) + ' ' + stringa)
    print '\n### tar file: ' + re.sub('raw.list', 'tar.gz', rawfile)


#################################################################
def repstringinfile(filein, fileout, string1, string2):
    # print "LOGX:: Entering `repstringinfile` method/function in
    # %(__file__)s" % globals()
    import re

    f = open(filein, 'r')
    ss = f.readlines()
    f.close()
    f = open(fileout, 'w')
    for n in range(len(ss)):
        if string1 in ss[n]:
            f.write(re.sub(string1, string2, ss[n]))
        else:
            f.write(ss[n])
    f.close()


###################################################

def StoN(img, ran=50):
    # print "LOGX:: Entering `StoN` method/function in %(__file__)s" %
    # globals()
    import ntt
    import numpy as np

    xx, yy = ntt.util.readspectrum(img)
    sntot = []
    xxmed = []
    for j in range(2, int((xx[-1] - xx[0]) / ran) - 2):
        #      aa,bb=xx[0]+j*ran-25,xx[0]+j*ran+25
        aa, bb = xx[0] + j * ran - \
            int(ran / 2.), xx[0] + j * ran + int(ran / 2.)
        ww = np.asarray([i for i in range(len(xx))
                         if ((xx[i] >= aa) & (xx[i] < bb))])
        snr = np.average(
            yy[ww]) / np.sqrt((sum(((yy[ww] - np.average(yy[ww]))) ** 2)) / (len(ww) - 1))
        sntot.append(snr)
        xxmed.append(xx[0] + j * ran)
    #   from pylab import ion,plot,show
    #   ion()
    #   plot(xxmed,sntot,'-r')
    #   raw_input('dd')
    #   show()
    #   plot(xx,yy,'-b')
    #   show()
    return np.mean(sntot)


def StoN2(img, show=False):
    # print "LOGX:: Entering `StoN2` method/function in %(__file__)s" %
    # globals()
    import numpy as np
    data, hdr0 = pyfits.getdata(img, header=True)
    yy1 = data[0][0]
    #      yy3=data[2][0]
    #      yy2=data[1][0]
    yy4 = data[3][0]
    xx = np.arange(len(data[0][0]))
    xxmed = hdr0['CRVAL1'] + xx * hdr0['CD1_1']
    sntot = yy1 / yy4
    if show:
        import pylab as pl

        pl.ion()
        pl.plot(xxmed, yy1, '-r', label='Spectrum')
        pl.plot(xxmed, sntot, '-b', label='StoN ' + str(np.mean(sntot)) + ' ')
        pl.legend(numpoints=1, markerscale=1.5, loc=1, ncol=1)
        pl.xlabel('Wavlength')
        pl.ylabel('StoN')
    return np.median(sntot)


################################################
def spectraresolution(img):
    # print "LOGX:: Entering `spectraresolution` method/function in
    # %(__file__)s" % globals()
    import ntt

    hdr = ntt.util.readhdr(img)
    _instrume = ntt.util.readkey3(hdr, 'instrume')
    _slit = ntt.util.readkey3(hdr, 'slit')
    _grism = ntt.util.readkey3(hdr, 'grism')
    _filter = ntt.util.readkey3(hdr, 'filter')
    risoluzioni = {}
    risoluzioni['efosc'] = {}
    risoluzioni['sofi'] = {}
    risoluzioni['efosc']['GG495', 'Gr13', 'slit1.5'] = 29.
    risoluzioni['efosc']['GG495', 'Gr13', 'slit1.0'] = 19.
    risoluzioni['efosc']['Free', 'Gr13', 'slit1.5'] = 29.
    risoluzioni['efosc']['Free', 'Gr13', 'slit1.0'] = 19.
    risoluzioni['efosc']['Free', 'Gr11', 'slit1.5'] = 22.
    risoluzioni['efosc']['Free', 'Gr11', 'slit1.0'] = 14.

    risoluzioni['efosc']['Free', 'Gr18', 'slit1.0'] = 8.19  # taken from Table 5 of the Manual.
    risoluzioni['efosc']['Free', 'Gr20', 'slit1.0'] = 2.0  # taken from Table 5 of the Manual, although this is actually for the 0.5" slit.

    risoluzioni['efosc']['OG530', 'Gr16', 'slit1.5'] = 22.
    risoluzioni['efosc']['OG530', 'Gr16', 'slit1.0'] = 16.
    risoluzioni['sofi']['GBF', 'GB', 'long_slit_1'] = 27.
    risoluzioni['sofi']['GBF', 'GB', 'long_slit_2'] = 35.
    risoluzioni['sofi']['GBR', 'GR', 'long_slit_1'] = 30.
    risoluzioni['sofi']['GBR', 'GR', 'long_slit_2'] = 38.
    if _instrume in risoluzioni.keys():
        if (_filter, _grism, _slit) in risoluzioni[_instrume].keys():
            return risoluzioni[_instrume][_filter, _grism, _slit]
        else:
            return ''
    else:
        return ''

################################################


def spectraresolution2(img0, ww=25):
    # print "LOGX:: Entering `spectraresolution2` method/function in
    # %(__file__)s" % globals()
    import string
    import re
    import ntt
    import numpy as np

    id = 'database/id' + re.sub('.fits', '', img0)
    img = 't' + img0
    data, hdr = pyfits.getdata(img, 0, header=True)
    crvals = ntt.util.readkey3(hdr, 'CRVAL2')
    cds = ntt.util.readkey3(hdr, 'CD2_2')
    ny = ntt.util.readkey3(hdr, 'NAXIS2')
    yy = data[:, int(ny / 2):int(ny / 2) + 10].mean(1)
    xx = np.arange(len(yy))
    zz = np.mean(yy)
    aa = crvals + (xx) * cds
    ntt.delete('new3.fits')
    hdu = pyfits.PrimaryHDU(yy)
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto('new3.fits')
########################################
#   use updateheader definition to make astropy compatible
#
#    hdulist[0].header.update('CRVAL1', crvals)
#    hdulist[0].header.update('CD1_1', cds)
#
    ntt.util.updateheader('new3.fits', 0, {'CRVAL1':[crvals,''],'CD1_1':[cds,'']})
    hdulist.close()
    #   read identified lines from id file
    f = open(id, 'r')
    ss = f.readlines()
    f.close()
    indices = [i for i, x in enumerate(ss) if "begin" in x]
    dd = ss[indices[1]:indices[2]]
    start = [i for i, x in enumerate(dd) if "features" in x][0] + 1
    stop = [i for i, x in enumerate(dd) if "function" in x][0]
    ff = dd[start:stop]
    lines = []
    for i in ff:
        lines.append(float(string.split(i)[2]))
    lines = np.compress((aa[0] < np.array(lines)) & (
        np.array(lines) < aa[-1]), np.array(lines))
    cursor = ''
    yym = np.interp(lines - ww, aa, yy)
    yyp = np.interp(lines + ww, aa, yy)
    for i in range(0, len(lines)):
        cursor = cursor + str(lines[i] - ww) + '  ' + str(yym[i]) + '  1   k\n'
        cursor = cursor + str(lines[i] + ww) + '  ' + str(yyp[i]) + '  1   k\n'
    cursor = cursor + str(lines[i] + ww) + '  ' + str(yyp[i]) + '  1   q\n'
    ff = open('_cursor', 'w')
    ff.write(cursor)
    ff.close()
    from pyraf import iraf

    aaa = iraf.noao.onedspec.bplot(
        'new3.fits', cursor='_cursor', spec2='', new_ima='', overwri='yes', Stdout=1)
    fw = []
    for i in aaa[1:]:
        fw.append(float(string.split(string.split(i, '=')[-1], 'k')[0]))
    ntt.delete('new3.fits,_cursor')
    # return mean(fw)
    return (aa[0] + ((aa[-1] - aa[0]) / 2)) / np.mean(fw)


##################################################

def limmag(img):
    # print "LOGX:: Entering `limmag` method/function in %(__file__)s" %
    # globals()
    import ntt
    import math
    import os
    import numpy as np

    hdr = ntt.util.readhdr(img)
    _ZP = ntt.util.readkey3(hdr, 'PHOTZP')
    _gain = ntt.util.readkey3(hdr, 'gain')
    _exptime = ntt.util.readkey3(hdr, 'exptime')
    _fwhm = ntt.util.readkey3(hdr, 'PSF_FWHM')
    _mbkg = ntt.util.readkey3(hdr, 'MBKG')  # background from sextractor
    _instrume = ntt.util.readkey3(hdr, 'instrume')

    if 'NCOMBINE' in hdr:
        _ncombine = ntt.util.readkey3(hdr, 'NCOMBINE')
    else:
        _ncombine = 1

    # CHANGE for DR2, EFRONN is defined after lim mag, need to compute again
    # here
    if _instrume == 'sofi':
        EFFRON = 12. * (math.sqrt(float(_ncombine)) /
                        math.sqrt(float(ntt.util.readkey3(hdr, 'ndit'))))
    else:
        if 'FLATCOR' in hdr and os.path.isfile(hdr['FLATCOR']):
            hdrn = ntt.util.readhdr(hdr['FLATCOR'])
            if 'NCOMBINE' in hdrn:
                nflat = hdrn['NCOMBINE']
            else:
                nflat = 1
        else:
            nflat = 1
        if 'ZEROCOR' in hdr and os.path.isfile(hdr['ZEROCOR']):
            hdrb = ntt.util.readhdr(hdr['ZEROCOR'])
            if 'NCOMBINE' in hdrb:
                nbias = hdrb['NCOMBINE']
            else:
                nbias = 1
        else:
            nbias = 1
        EFFRON = float(ntt.util.readkey3(hdr, 'ron')) * \
            math.sqrt(1. + 1. / float(nflat) + 1. / float(nbias))

    check = 1
    if not _ZP:
        check = 0
    elif _ZP == 9999:
        check = 0
    if not _gain:
        check = 0
    if not _fwhm:
        check = 0
    if not _mbkg:
        check = 0
    else:
        if float(_mbkg) <= 0:
            _mbkg = 1
    if check == 1:
        if _instrume == 'efosc':
            ps = ntt.util.readkey3(hdr, 'binx') * .12
        else:
            ps = 0.288
        npix = np.pi * (_fwhm / ps)
        _sn = 5  # signal to noise
        maglim = _ZP - 2.5 * np.log10(
            _sn / (_gain * _exptime) * (npix * (_mbkg * _gain) + (EFFRON ** 2) * npix) ** (1. / 2.))
        return maglim
    else:
        return ''

    # formula from McLean 1997)
#      n=pi*((_fwhm/ps)**2)
#      maglim=_ZP -2.5 * log10(sn * (1/_gain) * ((n*_mbkg/_exptime)**(.5)) )

####################################################################


def extractspectrum(img, dv, _ext_trace, _dispersionline, _interactive, _type, automaticex=False):
    # print "LOGX:: Entering `extractspectrum` method/function in
    # %(__file__)s" % globals()
    import glob
    import os
    import string
    import sys
    import re
    import ntt
    import datetime
    import numpy as np

    MJDtoday = 55927 + (datetime.date.today() - datetime.date(2012, 01, 01)).days
    from pyraf import iraf

    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.specred(_doprint=0)
    toforget = ['specred.apall', 'specred.transform']
    for t in toforget:
        iraf.unlearn(t)
    dv = ntt.dvex()
    hdr = ntt.util.readhdr(img)
    _gain = ntt.util.readkey3(hdr, 'gain')
    _rdnoise = ntt.util.readkey3(hdr, 'ron')
    _grism = ntt.util.readkey3(hdr, 'grism')
    iraf.specred.dispaxi = 2
    imgex = re.sub('.fits', '_ex.fits', img)
    imgfast = re.sub(string.split(img, '_')[-2] + '_', '', img)
    # imgfast=re.sub(str(MJDtoday)+'_','',img)
    if not os.path.isfile(imgex) and not os.path.isfile(
            'database/ap' + re.sub('.fits', '', img)) and not os.path.isfile(
            'database/ap' + re.sub('.fits', '', imgfast)):
        _new = 'yes'
        _extract = 'yes'
    else:
        if automaticex:
            if _interactive in ['Yes', 'yes', 'YES', 'y', 'Y']:
                answ = 'x'
                while answ not in ['o', 'n', 's']:
                    answ = raw_input(
                        '\n### New extraction [n], extraction with old parameters [o], skip extraction [s] ? [o]')
                    if not answ:
                        answ = 'o'
                if answ == 'o':
                    _new, _extract = 'no', 'yes'
                elif answ == 'n':
                    _new, _extract = 'yes', 'yes'
                else:
                    _new, _extract = 'yes', 'no'
            else:
                _new, _extract = 'no', 'yes'
        else:
            if _interactive in ['Yes', 'yes', 'YES', 'y', 'Y']:
                answ = 'x'
                while answ not in ['y', 'n']:
                    answ = raw_input(
                        '\n### do you want to extract again [[y]/n] ? ')
                    if not answ:
                        answ = 'y'
                if answ == 'y':
                    _new, _extract = 'yes', 'yes'
                else:
                    _new, _extract = 'yes', 'no'
            else:
                _new, _extract = 'yes', 'yes'
    if _extract == 'yes':
        ntt.util.delete(imgex)
        if _dispersionline:
            question = 'yes'
            while question == 'yes':
                _z1, _z2, goon = ntt.util.display_image(img, 1, '', '', False)
                dist = raw_input(
                    '\n### At which line do you want to extract the spectrum [' + str(dv['line'][_grism]) + '] ? ')
                if not dist:
                    dist = 400
                try:
                    dist = int(dist)
                    question = 'no'
                except:
                    print '\n### input not valid, try again:'
        else:
            dist = dv['line'][_grism]
        if _ext_trace in ['yes', 'Yes', 'YES', True]:
            lista = glob.glob('*ex.fits')
            if lista:
                for ii in lista:
                    print ii
                _reference = raw_input(
                    '\### which object do you want to use for the trace [' + str(lista[0]) + '] ? ')
                if not _reference:
                    _reference = lista[0]
                _reference = re.sub('_ex', '', _reference)
                _fittrac = 'no'
                _trace = 'no'
            else:
                sys.exit('\n### error: no extracted spectra in the directory')
        else:
            _reference = ''
            _fittrac = 'yes'
            _trace = 'yes'
        if _new == 'no':
            if not os.path.isfile('database/ap' + re.sub('.fits', '', img)):
                ntt.util.repstringinfile('database/ap' + re.sub('.fits', '', imgfast),
                                         'database/ap' +
                                         re.sub('.fits', '', img), re.sub(
                                             '.fits', '', imgfast),
                                         re.sub('.fits', '', img))
            _find = 'no'
            _recenter = 'no'
            _edit = 'no'
            _trace = 'no'
            _fittrac = 'no'
            _mode = 'h'
            _resize = 'no'
            _review = 'no'
            iraf.specred.mode = 'h'
            _interactive = 'no'
        else:
            iraf.specred.mode = 'q'
            _mode = 'q'
            _find = 'yes'
            _recenter = 'yes'
            _edit = 'yes'
            _review = 'yes'
            _resize = dv[_type]['_resize']
        iraf.specred.apall(img, output=imgex, referen=_reference, trace=_trace, fittrac=_fittrac, find=_find,
                           recenter=_recenter, edit=_edit,
                           nfind=1, extract='yes', backgro='fit', gain=_gain, readnoi=_rdnoise, lsigma=4, usigma=4,
                           format='multispec',
                           b_function='legendre', b_sample=dv[_type]['_b_sample'], clean='yes', pfit='fit1d',
                           lower=dv[_type]['_lower'], upper=dv[_type][
                               '_upper'], t_niter=dv[_type]['_t_niter'],
                           width=dv[_type]['_width'],
                           radius=dv[_type]['_radius'], line=dist, nsum=dv[
                               _type]['_nsum'], t_step=dv[_type]['_t_step'],
                           t_nsum=dv[_type]['_t_nsum'],
                           t_nlost=dv[_type]['_t_nlost'], t_sample=dv[
                               _type]['_t_sample'], resize=_resize,
                           t_order=dv[_type]['_t_order'],
                           weights=dv[_type]['_weights'], interactive=_interactive, review=_review, mode=_mode)
        ntt.util.repstringinfile('database/ap' + re.sub('.fits', '', img), 'database/ap' + re.sub('.fits', '', imgfast),
                                 re.sub('.fits', '', img), re.sub('.fits', '', imgfast))

        data, hdr = pyfits.getdata(imgex, 0, header=True)
        xxex = np.arange(len(data[0][0]))
        aaex = ntt.util.readkey3(hdr, 'CRVAL1') + \
            (xxex) * ntt.util.readkey3(hdr, 'CD1_1')
        # add sky from original image  for sofi .....probably better to move
        # out of this module
        _original = ntt.util.readkey3(hdr, 'ORIGFILE')
        _archive = ntt.util.readkey3(hdr, 'ARCFILE')
        _arc = ntt.util.readkey3(hdr, 'ARC')
        _instrume = ntt.util.readkey3(hdr, 'instrume')
        if _arc and _instrume != 'efosc':
            if os.path.isfile(_arc):
                if os.path.isfile(_archive):
                    imgstart = _archive
                elif os.path.isfile(_original):
                    imgstart = _original
                else:
                    imgstart = ''
                if imgstart:
                    ntt.util.delete('_tmp.fits')
                    iraf.specred.transform(input=imgstart, output='_tmp.fits', minput='',
                                           fitnames=re.sub('.fits', '', _arc), databas='database',
                                           x1='INDEF', x2='INDEF', y1='INDEF', y2='INDEF', flux='yes', mode='h',
                                           logfile='logfile')
                    yysky = pyfits.open('_tmp.fits')[0].data[:, 10]
                    ntt.util.delete('_tmp.fits')
                    data[2][0] = yysky
                    ntt.util.delete(imgex)
                    pyfits.writeto(imgex, np.float32(data), hdr)
                else:
                    print '\n### warning raw image for sky information not found'
                ##########################################
            #        ntt.util.updateheader(imgex,0,{'XMIN':[aaex[0],'min wavelength [Angstrom]'],'XMAX':[aaex[-1],'max wavelength [Angstrom]']})
    else:
        print '\n### skipping new extraction'
    return imgex

##########################################################################
