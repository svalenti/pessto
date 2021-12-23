import numpy as np

def xpa(arg):
    # print "LOGX:: Entering `xpa` method/function in %(__file__)s" % globals()
    import subprocess

    subproc = subprocess.Popen('xpaset -p ds9 ' + arg, shell=True)
    subproc.communicate()


def vizq(_ra, _dec, catalogue, radius):
    # print "LOGX:: Entering `vizq` method/function in %(__file__)s" %
    # globals()
    import os
    import string
    import re

    _site = 'vizier.u-strasbg.fr'
    # _site='vizier.cfa.harvard.edu'
    cat = {'usnoa2': ['I/252/out', 'USNO-A2.0', 'Rmag'], '2mass': ['II/246/out', '2MASS', 'Jmag'],
           'usnob1': ['I/284/out', 'USNO-B1.0', 'R2mag']}

    vizcommand ='vizquery -mime=tsv  ' + '-site=' + _site + ' -source=' + cat[catalogue][0] +\
                ' -c.ra=' + str(_ra) + ' -c.dec=' + str(_dec) + ' -c.eq=J2000 -c.rm=' + str(radius) +\
                ' -c.geom=b -oc.form=h -sort=_RA*-c.eq -out.add=_RAJ2000,_DEJ2000 -out.max=10000 -out=' +\
                cat[catalogue][1] + ' -out=' + cat[catalogue][2] + ''

    print(vizcommand)
    a = os.popen( vizcommand ).read()

    aa = a.split('\n')
    bb = []
    for i in aa:
        if i and i[0] != '#':
            bb.append(i)
    _ra, _dec, _name, _mag = [], [], [], []
    for ii in bb[3:]:
        aa = ii.split('\t')
        _ra.append(re.sub(' ', ':', aa[0]))
        _dec.append(re.sub(' ', ':', aa[1]))
        _name.append(aa[2])
        try:
            _mag.append(float(aa[3]))
        except:
            _mag.append(float(9999))
    dictionary = {'ra': _ra, 'dec': _dec, 'id': _name, 'mag': _mag}
    return dictionary

###########################################################################################

def vizq2(_ra, _dec, catalogue, radius):
    ''' Query vizquery '''
    import os, string, re
    #_site = 'vizier.cfa.harvard.edu'
    _site='vizier.u-strasbg.fr'
    cat = {'usnoa2': ['I/252/out', 'USNO-A2.0', 'Rmag'],
           'usnob1': ['I/284/out', 'USNO-B1.0', 'R2mag'],
           '2mass': ['II/246/out', '2MASS', 'Jmag,Hmag,Kmag'],
           'landolt': ['II/183A/table2', '', 'Vmag,B-V,U-B,V-R,R-I,Star,e_Vmag'],
           'apass': ['I/322A/out', '', 'Bmag,Vmag,gmag,rmag,imag,e_Vmag,e_Bmag,e_gmag,e_rmag,e_imag,UCAC4'],

           'sdss9': ['V/139/sdss9', '', 'objID,umag,gmag,rmag,imag,zmag,e_umag,e_gmag,e_rmag,e_imag,e_zmag,gc'],
           'sdss7': ['II/294/sdss7', '', 'objID,umag,gmag,rmag,imag,zmag,e_umag,e_gmag,e_rmag,e_imag,e_zmag,gc'],
           'sdss8': ['II/306/sdss8', '', 'objID,umag,gmag,rmag,imag,zmag,e_umag,e_gmag,e_rmag,e_imag,e_zmag,gc']}

    a = os.popen('vizquery -mime=tsv  -site=' + _site + ' -source=' + cat[catalogue][0] + \
                 ' -c.ra=' + str(_ra) + ' -c.dec=' + str(_dec) + ' -c.eq=J2000 -c.rm=' + str(radius) + \
                 ' -c.geom=b -oc.form=h -sort=_RA*-c.eq -out.add=_RAJ2000,_DEJ2000 -out.max=10000 -out=' + \
                 cat[catalogue][1] + ' -out=' + cat[catalogue][2] + '').read()
    print('vizquery -mime=tsv  -site=' + _site + ' -source=' + cat[catalogue][0] + \
          ' -c.ra=' + str(_ra) + ' -c.dec=' + str(_dec) + ' -c.eq=J2000 -c.rm=' + str(radius) + \
          ' -c.geom=b -oc.form=h -sort=_RA*-c.eq -out.add=_RAJ2000,_DEJ2000 -out.max=10000 -out=' + \
          cat[catalogue][1] + ' -out=' + cat[catalogue][2] + '')

    aa = a.split('\n')

    bb = []
    for i in aa:
        if i and i[0] != '#':
            bb.append(i)
    _ra, _dec, _name, _mag = [], [], [], []
    for ii in bb[3:]:
      if ii:
        aa = ii.split('\t')
        rr, dd = deg2HMS(ra=re.sub(' ', ':', aa[0]), dec=re.sub(' ', ':', aa[1]), round=False)
        _ra.append(rr)
        _dec.append(dd)
        _name.append(aa[2])

    dictionary = {'ra': _ra, 'dec': _dec, 'id': _name}
    sss = string.split(cat[catalogue][2], ',')
    for ii in sss:
        dictionary[ii] = []
    for ii in bb[3:]:
      if ii:
        aa = ii.split('\t')
        for gg in range(0, len(sss)):
            if sss[gg] not in ['UCAC4', 'id']:
                try:
                    dictionary[sss[gg]].append(float(aa[3 + gg]))
                except Exception as e:
                    print(e)
                    dictionary[sss[gg]].append(float(9999))
            else:
                dictionary[sss[gg]].append(str(aa[3 + gg]))

    if catalogue in ['sdss7', 'sdss9', 'sdss8']:
        dictionary['u'] = dictionary['umag']
        dictionary['g'] = dictionary['gmag']
        dictionary['r'] = dictionary['rmag']
        dictionary['i'] = dictionary['imag']
        dictionary['z'] = dictionary['zmag']
        dictionary['uerr'] = dictionary['e_umag']
        dictionary['gerr'] = dictionary['e_gmag']
        dictionary['rerr'] = dictionary['e_rmag']
        dictionary['ierr'] = dictionary['e_imag']
        dictionary['zerr'] = dictionary['e_zmag']
        for key in dictionary.keys():
            if key != 'r':
                dictionary[key] = np.compress((np.array(dictionary['r']) < 19) & (np.array(dictionary['r'] > 10)),
                                              dictionary[key])
        dictionary['r'] = np.compress((np.array(dictionary['r']) < 19) & (np.array(dictionary['r'] > 10)),
                                      dictionary['r'])

    elif catalogue == 'landolt':
        dictionary['B'] = np.array(dictionary['Vmag']) + np.array(dictionary['B-V'])
        dictionary['U'] = np.array(dictionary['B']) + np.array(dictionary['U-B'])
        dictionary['V'] = np.array(dictionary['Vmag'])
        dictionary['Verr'] = np.array(dictionary['e_Vmag'])
        dictionary['R'] = np.array(dictionary['Vmag']) - np.array(dictionary['V-R'])
        dictionary['I'] = np.array(dictionary['R']) - np.array(dictionary['R-I'])
        dictionary['id'] = np.array(dictionary['Star'])
    elif catalogue == 'apass':
        dictionary['B'] = np.array(dictionary['Bmag'])
        dictionary['V'] = np.array(dictionary['Vmag'])
        dictionary['g'] = np.array(dictionary['gmag'])
        dictionary['r'] = np.array(dictionary['rmag'])
        dictionary['i'] = np.array(dictionary['imag'])
        dictionary['Berr'] = np.array(dictionary['e_Bmag'], float) / 100.
        dictionary['Verr'] = np.array(dictionary['e_Vmag'], float) / 100.
        dictionary['gerr'] = np.array(dictionary['e_gmag'], float) / 100.
        dictionary['rerr'] = np.array(dictionary['e_rmag'], float) / 100.
        dictionary['ierr'] = np.array(dictionary['e_imag'], float) / 100.
        dictionary['id'] = np.array(dictionary['UCAC4'], str)
        for key in dictionary.keys():
            if key != 'r':
                dictionary[key] = np.compress((np.array(dictionary['r']) < 22) & (np.array(dictionary['r'] > 10.5)),
                                              dictionary[key])
        dictionary['r'] = np.compress((np.array(dictionary['r']) < 22) & (np.array(dictionary['r'] > 10.5)),
                                      dictionary['r'])
    elif catalogue == '2mass':
        dictionary['J'] = np.array(dictionary['Jmag'])
        dictionary['H'] = np.array(dictionary['Hmag'])
        dictionary['K'] = np.array(dictionary['Kmag'])
        dictionary['mag'] = dictionary['Jmag']
        dictionary['mag1'] = dictionary['Jmag']
        dictionary['mag2'] = dictionary['Hmag']
        dictionary['mag3'] = dictionary['Kmag']
    elif catalogue == 'usnoa2':
        dictionary['mag'] = dictionary['Rmag']
    elif catalogue == 'usnob1':
        dictionary['mag'] = dictionary['R2mag']
    return dictionary

########################################################################################

def wcsstart(img, CRPIX1='', CRPIX2=''):
    # print "LOGX:: Entering `wcsstart` method/function in %(__file__)s" %
    # globals()
    from numpy import pi, sin, cos
    from ntt.util import updateheader, readhdr, readkey3
    import ntt

    hdr = readhdr(img)
    _instrume = readkey3(hdr, 'instrume')
    _RA = readkey3(hdr, 'RA')
    _DEC = readkey3(hdr, 'DEC')
    _xdimen = readkey3(hdr, 'NAXIS1')
    _ydimen = readkey3(hdr, 'NAXIS2')
    angle = readkey3(hdr, 'posang')
    if _instrume == 'efosc':
        theta = (angle * pi / 180.)
        CDELT0 = 6.6888889999999995e-05
        if not CRPIX1:
            CRPIX1 = 506.651
        else:
            CRPIX1 = 506.651 + float(CRPIX1)
        if not CRPIX2:
            CRPIX2 = 488.909
        else:
            CRPIX2 = 488.909 + float(CRPIX2)
            # CDELT1=2.
        #        CDELT2=2.
        CD1_1 = CDELT0 * cos(theta)
        CD2_2 = CDELT0 * cos(theta)
        CD1_2 = abs(CDELT0) * (abs(CDELT0) / CDELT0) * sin(theta)
        CD2_1 = (-1) * abs(CDELT0) * (abs(CDELT0) / CDELT0) * sin(theta)
    elif _instrume == 'sofi':
        theta = ((angle + 270.) * pi / 180.)
        CDELT0 = 8.0000000000000E-5  # 7.988960e-05#8.00062491264e-05
        if not CRPIX1:
            _CRPIX1 = readkey3(hdr, 'NTCRPIX1')
            if _CRPIX1:
                CRPIX1 = float(_CRPIX1)
            else:
                CRPIX1 = 512.
        else:
            CRPIX1 = float(CRPIX1) + 512

        if not CRPIX2:
            _CRPIX2 = readkey3(hdr, 'NTCRPIX2')  # +10
            if _CRPIX2:
                CRPIX2 = float(_CRPIX2)
            else:
                CRPIX2 = 512.
        else:
            CRPIX2 = float(CRPIX2) + 521
        # CDELT1=1.
        #        CDELT2=1.
        CD1_1 = CDELT0 * cos(theta)
        CD2_2 = (-1) * CDELT0 * cos(theta)
        CD1_2 = abs(CDELT0) * (abs(CDELT0) / CDELT0) * sin(theta)
        CD2_1 = abs(CDELT0) * (abs(CDELT0) / CDELT0) * sin(theta)

    CTYPE1 = 'RA---TAN'
    CTYPE2 = 'DEC--TAN'
    CRVAL1 = _RA
    CRVAL2 = _DEC
    WCSDIM = 2
    LTM1_1 = 1.
    LTM2_2 = 1.
    WAT0_001 = 'system=image'
    WAT1_001 = 'wtype=tan axtype=ra'
    WAT2_001 = 'wtype=tan axtype=dec'
    ntt.util.updateheader(img, 0, {'CTYPE1': [CTYPE1, 'pixel coordinate system'],
                                   'CTYPE2': [CTYPE2, 'pixel coordinate system'],
                                   'CRVAL1': [CRVAL1, 'RA at ref pixel'],
                                   'CRVAL2': [CRVAL2, 'DEC at ref pixel'],
                                   'CRPIX1': [CRPIX1, 'Ref pixel in X'],
                                   'CRPIX2': [CRPIX2, 'Ref pixel in Y'],
                                   'CD1_1': [CD1_1, 'Transformation matrix element'],
                                   'CD2_2': [CD2_2, 'Transformation matrix element'],
                                   'CD1_2': [CD1_2, 'Transformation matrix element'],
                                   'CD2_1': [CD2_1, 'Transformation matrix element'],
                                   'WCSDIM': [WCSDIM, '']})


def efoscastroloop(imglist, catalogue, _interactive, number1, number2, number3, _fitgeo, _tollerance1,
                   _tollerance2, sexvec='', _guess=False, _numin=4, method='iraf', _CRPIX1='', _CRPIX2=''):
    # print "LOGX:: Entering `efoscastroloop` method/function in %(__file__)s"
    # % globals()
    import ntt
    from ntt.util import delete, readkey3, readhdr
    from ntt.efoscastrodef import efoscastrometry2
    from numpy import median, array
    import math
    import datetime
    import time
    # from datetime import datetime
    _imex = False
    for img in imglist:
        hdr = readhdr(img)
        _instrume = readkey3(hdr, 'instrume')
        if catalogue == 'inst' and _instrume == 'efosc':
            catalogue = 'usnob1'
        if catalogue == 'inst' and _instrume == 'sofi':
            catalogue = '2mass'
        if not sexvec:
            sexvec = ntt.efoscastrodef.sextractor(img)
        ###################
        if _guess:
            print('guess astrometry before starting ')
            ntt.efoscastrodef.wcsstart(img, _CRPIX1, _CRPIX2)
        ss = datetime.datetime.now()
        time.sleep(1)
        catvec = ntt.efoscastrodef.querycatalogue(catalogue, img, method)
        rmsx1, rmsy1, num1, fwhm1, ell1, ccc, bkg1, rasys1, decsys1 = efoscastrometry2(
            [img], catalogue, _interactive, number1, sexvec, catvec, guess=False, fitgeo=_fitgeo,
            tollerance1=_tollerance1, tollerance2=_tollerance2, _update='yes', imex=_imex, nummin=_numin)
        if rmsx1 > 1 or rmsy1 > 1:
            catvec = ntt.efoscastrodef.querycatalogue(catalogue, img, method)
            rmsx2, rmsy2, num2, fwhm2, ell2, ccc, bkg2, rasys2, decsys2 = efoscastrometry2(
                [img], catalogue, _interactive, number2, sexvec, catvec, guess=False, fitgeo=_fitgeo,
                tollerance1=_tollerance1, tollerance2=_tollerance2, _update='yes', imex=_imex, nummin=_numin)
            if rmsx2 > 1 or rmsy2 > 1:
                catvec = ntt.efoscastrodef.querycatalogue(
                    catalogue, img, method)
                rmsx3, rmsy3, num3, fwhm3, ell3, ccc, bkg3, rasys3, decsys3 = efoscastrometry2(
                    [img], catalogue, _interactive, number3, sexvec, catvec, guess=False, fitgeo=_fitgeo,
                    tollerance1=_tollerance1, tollerance2=_tollerance2, _update='yes', imex=_imex, nummin=_numin)
            else:
                rmsx3, rmsy3, num3, fwhm3, ell3, ccc, bkg3, rasys3, decsys3 = rmsx2, rmsy2, num2, fwhm2, ell2, ccc, bkg2, rasys2, decsys2
        else:
            rmsx3, rmsy3, num3, fwhm3, ell3, ccc, bkg3, rasys3, decsys3 = rmsx1, rmsy1, num1, fwhm1, ell1, ccc, bkg1, rasys1, decsys1
        ########################################
        if rmsx3 < 10 and rmsy3 < 10:
            if _instrume == 'efosc':
                fwhmgess3 = median(array(fwhm3)) * .68 * \
                    2.35 * readkey3(hdr, 'binx') * .12
                if _imex:
                    fwhmgessime = median(array(ccc)) * \
                        readkey3(hdr, 'binx') * .12
                else:
                    fwhmgessime = 9999
            elif _instrume == 'sofi':
                fwhmgess3 = median(array(fwhm3)) * .68 * 2.35 * 0.288
                if _imex:
                    fwhmgessime = median(array(ccc)) * 0.288
                else:
                    fwhmgessime = 9999
            ellgess3 = median(array(ell3))
        else:
            fwhmgess3 = 9999
            fwhmgessime = 9999
            ellgess3 = 9999
        if _instrume == 'efosc':
            mbkg3 = median(bkg3)
            ntt.util.updateheader(
                img, 0, {'MBKG': [mbkg3, 'background level']})
        else:
            mbkg3 = readkey3(hdr, 'MBKG')
        if fwhmgess3:
            print(fwhmgess3)
            ###################################################################
            # change saturation mag    2014-05-18
            ###################################################################
            if _instrume == 'efosc':
                magsat = -2.5 * math.log10(
                    (math.pi / (4 * math.log(2.))) * (60000. - float(mbkg3)) * ((float(fwhmgess3) / 0.24) ** 2))
            else:
                magsat = -2.5 * math.log10(
                    (math.pi / (4 * math.log(2.))) * (32000. - float(mbkg3)) * ((float(fwhmgess3) / 0.288) ** 2))
        else:
            magsat = 9999
    print(rmsx3, rmsy3, num3, fwhmgess3, ellgess3, fwhmgessime, rasys3, decsys3, magsat)
    if catalogue == '2mass':
        rasys3, decsys3 = 0.0000278, 0.0000278
    elif catalogue in ['usnoa2', 'usnob1']:
        rasys3, decsys3 = 0.0000556, 0.0000556
    return rmsx3, rmsy3, num3, fwhmgess3, ellgess3, fwhmgessime, rasys3, decsys3, magsat

    #if _instrume=='efosc':    V=(math.pi/(4*math.log(2)))*(60000-float(mbkg3))*(float(fwhmgess3)**2)
    # else:                     V=(math.pi/(4*math.log(2)))*(32000-float(mbkg3))*(float(fwhmgess3)**2)
    # magsat=-2.5*math.log10(V)


# ###########################################################################

def readtxt(ascifile):
    # print "LOGX:: Entering `readtxt` method/function in %(__file__)s" %
    # globals()
    import string

    f = open(ascifile, 'r')
    ss = f.readlines()
    f.close()
    columnname = []
    for i in range(0, len(ss)):
        if ss[i][0] == '#' and 'nfields' in ss[i]:
            num = int(string.split(ss[i])[-1])
            for g in range(1, num + 1):
                columnname.append(string.split(ss[i + g])[1])
            break
    ascidic = {}
    #    try:
    #        from numpy import genfromtxt
    #        data=genfromtxt(ascifile,str,unpack=True)
    #        for i in range(0,len(columnname)):
    #                 ascidic[columnname[i]]=data[i]
    #    except:
    for j in columnname:
        ascidic[j] = []
    for i in range(0, len(ss)):
        if ss[i][0] != '#':
            for j in range(0, len(columnname)):
                ascidic[columnname[j]].append(string.split(ss[i])[j])
    return ascidic


#########################################################################

def zeropoint(img, _field, method='iraf', verbose=False, _interactive=False):
    # print "LOGX:: Entering `zeropoint` method/function in %(__file__)s" %
    # globals()
    import string
    import os
    import re
    import sys
    from numpy import compress, array, median, zeros, std, mean, abs
    import math
    import ntt
    from ntt.util import delete, readhdr, readkey3
    from ntt import sqlcl
    from pyraf import iraf

    iraf.noao(_doprint=0, Stdout=0)
    iraf.digiphot(_doprint=0, Stdout=0)
    iraf.daophot(_doprint=0, Stdout=0)
    iraf.images(_doprint=0, Stdout=0)
    iraf.imcoords(_doprint=0, Stdout=0)
    iraf.proto(_doprint=0, Stdout=0)
    hdr = readhdr(img)
    _airmass = readkey3(hdr, 'AIRMASS')
    if not _airmass:
        print('\n### warning: airmass at starting exposure')
        _airmass = readkey3(hdr, 'airmass')
    _exptime = readkey3(hdr, 'exptime')
    _filter = readkey3(hdr, 'filter')
    _instrume = readkey3(hdr, 'instrume')
    _date = readkey3(hdr, 'date-night')
    _object = readkey3(hdr, 'object')

    #
    #  since the NDIT images are average for sofi
    #  the exptime to measure the object that should be used is DIT
    #
    if _instrume == 'sofi':
        _exptime = readkey3(hdr, 'DIT')
    else:
        _exptime = readkey3(hdr, 'exptime')

    kk = {'U': 0.46, 'u': 0.46, 'B': 0.27, 'g': 0.20, 'V': 0.12, 'r': 0.09, 'R': 0.09, 'i': 0.02, 'I': 0.02, 'z': 0.03,
          'J': 0.0, 'H': 0.0, 'K': 0.0}
    # zz = {'U640': 23.69, 'B639': 25.83, 'V641': 25.88, 'R642': 25.98, 'r784': 25.27, 'i705': 25.13, 'g782': 25.86,
    #      'z623': 24.41, 'J': 25, 'H': 25, 'Ks': 25}
    zz = {'U640': 23.655, 'B639': 25.755, 'V641': 25.830, 'R642': 25.967, 'r784': 25.673, 'i705': 25.112,
          'g782': 25.897, 'z623': 24.777, 'J': 25, 'H': 25, 'Ks': 25}

    if _filter in ['J', 'Js', 'H', 'Ks', 'K']:
        _field = '2mass'
        filters = {'J': 'J', 'Js': 'J', 'H': 'H', 'K': 'K', 'Ks': 'K'}
        colors = {'J': ['JH'], 'H': ['JH'], 'K': ['HK']}

        if verbose:
            print('Infrared image : J ,H, K or Ks')

        iraf.astcat(_doprint=0, Stdout=0)
        iraf.imcoords(_doprint=0, Stdout=0)
        iraf.noao.astcat.aregpars.rcrauni = ''
        iraf.noao.astcat.aregpars.rcdecuni = ''
        iraf.noao.astcat.catdb = ntt.__path__[
            0] + '/standard/cat/catalogue.dat'
        _ra = readkey3(hdr, 'RA')
        _dec = readkey3(hdr, 'DEC')
        stdcoo = querycatalogue('2mass', img, method)
        rastd, decstd = array(stdcoo['ra'], float), array(stdcoo['dec'], float)
        standardpix = {'ra': stdcoo['x'],
                       'dec': stdcoo['y'], 'id': stdcoo['id']}
        if verbose:
            print(stdcoo.keys())
            print(stdcoo['J'])
            print(stdcoo['H'])
            print(stdcoo['K'])
            print(stdcoo['mag1'])
            print(stdcoo['mag2'])
            print(stdcoo['mag3'])

    else:
        # check if it is landolt field
        stdcooL = ntt.efoscastrodef.readtxt(
            ntt.__path__[0] + '/standard/cat/landolt.cat')
        rastdL, decstdL = array(stdcooL['ra'], float), array(
            stdcooL['dec'], float)
        delete('tmp.stdL.pix')
        iraf.wcsctran(ntt.__path__[0] + '/standard/cat/landolt.cat', 'tmp.stdL.pix', img, inwcs='world',
                      units='degrees degrees', outwcs='logical', columns='1 2', formats='%10.1f %10.1f', verbose='no')
        standardpixL = ntt.efoscastrodef.readtxt('tmp.stdL.pix')
        xstdL = standardpixL['ra']
        ystdL = standardpixL['dec']
        idstdL = standardpixL['id']
        xstdL = compress(
            (array(xstdL, float) < readkey3(hdr, 'naxis1')) & (array(xstdL, float) > 0) & (array(ystdL, float) > 0) & (
                array(ystdL, float) < readkey3(hdr, 'naxis2')), xstdL)
        # check if it is sloan field
        _ra = readkey3(hdr, 'RA')
        _dec = readkey3(hdr, 'DEC')
        magsel0, magsel1 = 12, 18
        _ids = ntt.efoscastrodef.sloan2file(
            _ra, _dec, 4, float(magsel0), float(magsel1), '_tmpsloan.cat')

        ascifile = '_tmpsloan.cat'
        stdcooS = ntt.efoscastrodef.readtxt(ascifile)
        rastdS, decstdS = array(stdcooS['ra'], float), array(
            stdcooS['dec'], float)
        delete('tmp.stdS.pix')
        iraf.wcsctran(ascifile, 'tmp.stdS.pix', img, inwcs='world', units='degrees degrees', outwcs='logical',
                      columns='1 2', formats='%10.1f %10.1f', verbose='no')
        standardpixS = ntt.efoscastrodef.readtxt('tmp.stdS.pix')
        xstdS = standardpixS['ra']
        ystdS = standardpixS['dec']
        idstdS = standardpixS['id']
        xstdS = compress(
            (array(xstdS, float) < readkey3(hdr, 'naxis1')) & (array(xstdS, float) > 0) &
            (array(ystdS, float) > 0) & (array(ystdS, float) < readkey3(hdr, 'naxis2')), xstdS)

        ##############
        if _filter in ['U640', 'B639', 'V641', 'R642']:
            if _field == 'sloan':
                standardpix, stdcoo = {'ra': [9999], 'dec': [
                    9999], 'id': [1]}, {'ra': [9999], 'dec': [9999]}
            # standardpix={'ra':[9999],'dec':[9999],'id':[1]}
            else:
                _field = 'landolt'
                filters = {'U640': 'U', 'B639': 'B',
                           'V641': 'V', 'R642': 'R', 'i705': 'I'}
                colors = {'U': ['UB'], 'B': ['BV'], 'V': [
                    'BV', 'VR'], 'R': ['VR', 'RI'], 'I': ['RI']}
                if len(xstdL) >= 1:
                    standardpix = standardpixL
                    stdcoo = stdcooL
                elif len(xstdS) >= 1:
                    standardpix = standardpixS
                    stdcoo = stdcooS
                    stdcoo = ntt.efoscastrodef.transformsloanlandolt(stdcoo)
                    print('\n### transform sloan in landolt')
                else:
                    standardpix, stdcoo = {'ra': [9999], 'dec': [
                        9999], 'id': [1]}, {'ra': [9999], 'dec': [9999]}
        elif _filter in ['g782', 'r784', 'z623']:
            if _field == 'landolt':
                standardpix, stdcoo = {'ra': [9999], 'dec': [
                    9999], 'id': [1]}, {'ra': [9999], 'dec': [9999]}
            else:
                _field = 'sloan'
                filters = {'i705': 'i', 'g782': 'g', 'r784': 'r', 'z623': 'z'}
                colors = {'i': ['ri'], 'r': ['ri'], 'g': ['gr'], 'z': ['iz']}
                if len(xstdS) >= 1:
                    standardpix = standardpixS
                    stdcoo = stdcooS
                elif len(xstdL) >= 1:
                    standardpix = standardpixL
                    stdcoo = stdcooL
                    stdcoo = ntt.efoscastrodef.transformlandoltsloan(stdcoo)
                    print('\n### transform landolt to sloan')
                else:
                    standardpix, stdcoo = {'ra': [9999], 'dec': [
                        9999], 'id': [1]}, {'ra': [9999], 'dec': [9999]}
        elif _filter in ['i705']:
            if len(xstdL) >= 1 and _field != 'sloan':
                _field = 'landolt'
                filters = {'U640': 'U', 'B639': 'B',
                           'V641': 'V', 'R642': 'R', 'i705': 'I'}
                colors = {'U': ['UB'], 'B': ['BV'], 'V': [
                    'BV', 'VR'], 'R': ['VR', 'RI'], 'I': ['RI']}
                standardpix = standardpixL
                stdcoo = stdcooL
            elif len(xstdS) >= 1 and len(xstdL) == 0 and _field != 'sloan':
                _field = 'landolt'
                filters = {'U640': 'U', 'B639': 'B',
                           'V641': 'V', 'R642': 'R', 'i705': 'I'}
                colors = {'U': ['UB'], 'B': ['BV'], 'V': [
                    'BV', 'VR'], 'R': ['VR', 'RI'], 'I': ['RI']}
                standardpix = standardpixS
                stdcoo = stdcooS
                print('\n### transform sloan in landolt')
                stdcoo = ntt.efoscastrodef.transformsloanlandolt(stdcoo)
            elif len(xstdS) >= 1 and len(xstdL) == 0 and _field == 'sloan':
                _field = 'sloan'
                filters = {'i705': 'i', 'g782': 'g', 'r784': 'r', 'z623': 'z'}
                colors = {'i': ['ri'], 'r': ['ri'], 'g': ['gr'], 'z': ['iz']}
                standardpix = standardpixS
                stdcoo = stdcooS
            elif len(xstdL) >= 1 and len(xstdS) == 0 and _field == 'sloan':
                _field = 'sloan'
                filters = {'i705': 'i', 'g782': 'g', 'r784': 'r', 'z623': 'z'}
                colors = {'i': ['ri'], 'r': ['ri'], 'g': ['gr'], 'z': ['iz']}
                standardpix = standardpixL
                stdcoo = stdcooL
                print('\n### transform landolt in sloan')
                stdcoo = ntt.efoscastrodef.transformlandoltsloan(stdcoo)
            else:
                standardpix, stdcoo = {'ra': [9999], 'dec': [
                    9999], 'id': [1]}, {'ra': [9999], 'dec': [9999]}

    xstd = standardpix['ra']
    ystd = standardpix['dec']
    idstd = standardpix['id']
    rastd, decstd = array(stdcoo['ra'], float), array(stdcoo['dec'], float)
    xstd0 = compress(
        (array(xstd, float) < readkey3(hdr, 'naxis1')) & (array(xstd, float) > 0) & (array(ystd, float) > 0) & (
            array(ystd, float) < readkey3(hdr, 'naxis2')), xstd)
    if len(xstd0) > 1:  # go only if standard stars are in the field  ##########
        magstd0 = {}
        airmass0 = {}
        print('\n###  standard field: ' + str(_field))
        # sextractor on standard field
        namesex = ntt.util.defsex('default.sex')
        os.system('sex ' + img + ' -c ' + namesex + ' > _logsex')
        delete(namesex)
        delete('_logsex')
        xsex = iraf.proto.fields('detections.cat', fields='2', Stdout=1)
        ysex = iraf.proto.fields('detections.cat', fields='3', Stdout=1)
        fw = iraf.proto.fields('detections.cat', fields='8', Stdout=1)
        f = open('detection_sex.pix', 'w')
        for i in range(0, len(xsex)):
            f.write(str(xsex[i]) + ' ' + str(ysex[i]) + '\n')
        f.close()
        delete('detection_sex.coo')
        iraf.wcsctran('detection_sex.pix', 'detection_sex.coo', img, inwcs='logical', units='degrees degrees',
                      outwcs='world', columns='1 2', formats='%10.8f %10.8f')
        rasex = compress(array(iraf.proto.fields('detection_sex.coo', fields='1', Stdout=1)) != '',
                         array(iraf.proto.fields('detection_sex.coo', fields='1', Stdout=1)))
        decsex = compress(array(iraf.proto.fields('detection_sex.coo', fields='2', Stdout=1)) != '',
                          array(iraf.proto.fields('detection_sex.coo', fields='2', Stdout=1)))
        rasex = array(rasex, float)
        decsex = array(decsex, float)
        result = {}
        fileph = {}
        ystd0 = compress(
            (array(xstd, float) < readkey3(hdr, 'naxis1')) & (
                array(xstd, float) > 0) & (array(ystd, float) > 0)
            & (array(ystd, float) < readkey3(hdr, 'naxis2')), ystd)
        rastd0 = compress(
            (array(xstd, float) < readkey3(hdr, 'naxis1')) & (
                array(xstd, float) > 0) & (array(ystd, float) > 0)
            & (array(ystd, float) < readkey3(hdr, 'naxis2')), rastd)
        decstd0 = compress(
            (array(xstd, float) < readkey3(hdr, 'naxis1')) & (
                array(xstd, float) > 0) & (array(ystd, float) > 0)
            & (array(ystd, float) < readkey3(hdr, 'naxis2')), decstd)
        idstd0 = compress(
            (array(xstd, float) < readkey3(hdr, 'naxis1')) & (
                array(xstd, float) > 0) & (array(ystd, float) > 0)
            & (array(ystd, float) < readkey3(hdr, 'naxis2')), idstd)
        ###################
        colorvec = colors[filters[_filter]]
        if _field == 'landolt':
            print('\n###  landolt system')
            for _filtlandolt in 'UBVRI':
                if _filtlandolt == filters[_filter]:
                    airmass0[_filtlandolt] = _airmass
                else:
                    airmass0[_filtlandolt] = 1
                magstd0[_filtlandolt] = compress(
                    (array(xstd, float) < readkey3(hdr, 'naxis1')) & (array(xstd, float) > 0) &
                    (array(ystd, float) > 0) & (array(ystd, float) < readkey3(hdr, 'naxis2')), stdcoo[_filtlandolt])
            fileph['mU'] = zeros(len(rastd0)) + 999
            fileph['mB'] = zeros(len(rastd0)) + 999
            fileph['mV'] = zeros(len(rastd0)) + 999
            fileph['mR'] = zeros(len(rastd0)) + 999
            fileph['mI'] = zeros(len(rastd0)) + 999
            fileph['V'] = magstd0['V']
            fileph['BV'] = array(array(magstd0['B'], float) -
                                 array(magstd0['V'], float), str)
            fileph['UB'] = array(array(magstd0['U'], float) -
                                 array(magstd0['B'], float), str)
            fileph['VR'] = array(array(magstd0['V'], float) -
                                 array(magstd0['R'], float), str)
            fileph['RI'] = array(array(magstd0['R'], float) -
                                 array(magstd0['I'], float), str)
        elif _field == 'sloan':
            for _filtsloan in 'ugriz':
                if _filtsloan == filters[_filter]:
                    airmass0[_filtsloan] = _airmass
                else:
                    airmass0[_filtsloan] = 1
                magstd0[_filtsloan] = compress(
                    (array(xstd, float) < readkey3(hdr, 'naxis1')) & (array(xstd, float) > 0) &
                    (array(ystd, float) > 0) & (array(ystd, float) < readkey3(hdr, 'naxis2')), stdcoo[_filtsloan])
            fileph['mu'] = zeros(len(rastd0)) + 999
            fileph['mg'] = zeros(len(rastd0)) + 999
            fileph['mr'] = zeros(len(rastd0)) + 999
            fileph['mi'] = zeros(len(rastd0)) + 999
            fileph['mz'] = zeros(len(rastd0)) + 999
            fileph['r'] = magstd0['r']
            fileph['gr'] = array(array(magstd0['g'], float) -
                                 array(magstd0['r'], float), str)
            fileph['ri'] = array(array(magstd0['r'], float) -
                                 array(magstd0['i'], float), str)
            fileph['ug'] = array(array(magstd0['u'], float) -
                                 array(magstd0['g'], float), str)
            fileph['iz'] = array(array(magstd0['i'], float) -
                                 array(magstd0['z'], float), str)
        elif _field == '2mass':
            print('\n###  2mass system')
            for _filtlandolt in 'JHK':
                if _filtlandolt == filters[_filter]:
                    airmass0[_filtlandolt] = _airmass
                else:
                    airmass0[_filtlandolt] = 1
                magstd0[_filtlandolt] = compress(
                    (array(xstd, float) < readkey3(hdr, 'naxis1')) & (array(xstd, float) > 0) &
                    (array(ystd, float) > 0) & (array(ystd, float) < readkey3(hdr, 'naxis2')), stdcoo[_filtlandolt])

            fileph['mJ'] = zeros(len(magstd0['J'])) + 999
            fileph['mH'] = zeros(len(magstd0['H'])) + 999
            fileph['mK'] = zeros(len(magstd0['K'])) + 999
            fileph['J'] = magstd0['J']
            fileph['JH'] = array(array(magstd0['J'], float) -
                                 array(magstd0['H'], float), str)
            fileph['HK'] = array(array(magstd0['H'], float) -
                                 array(magstd0['K'], float), str)
        distvec, pos0, pos1 = ntt.efoscastrodef.crossmatch(array(rastd0), array(decstd0), array(rasex), array(decsex),
                                                           10)
        #    http://www.astromatic.net/forum/showthread.php?tid=516
        # if multiply FLUX_RADIUS by 1.9-2.0 you get a value of the fwhm
        #        fwhm0=(median(array(fw,float)[pos1]))*.68*2.35  #
        fwhm0 = (median(array(fw, float)[pos1])) * 1.9
        iraf.noao.digiphot.mode = 'h'
        iraf.noao.digiphot.daophot.photpars.zmag = 0
        iraf.noao.digiphot.daophot.daopars.psfrad = fwhm0 * 4
        iraf.noao.digiphot.daophot.daopars.fitrad = fwhm0
        iraf.noao.digiphot.daophot.daopars.sannulus = fwhm0 * 4
        iraf.noao.digiphot.daophot.fitskypars.annulus = fwhm0 * 4
        iraf.noao.digiphot.daophot.photpars.apertures = fwhm0 * 3
        iraf.noao.digiphot.daophot.datapars.readnoi = readkey3(hdr, 'ron')
        iraf.noao.digiphot.daophot.datapars.epadu = readkey3(hdr, 'gain')
        iraf.noao.digiphot.daophot.datapars.datamin = -1000
        iraf.noao.digiphot.daophot.datapars.datamax = 60000
        iraf.noao.digiphot.daophot.datapars.exposure = 'EXPTIME'
        iraf.noao.digiphot.daophot.datapars.airmass = 'AIRMASS'
        iraf.noao.digiphot.daophot.datapars.filter = 'FILTER'
        iraf.noao.digiphot.daophot.daopars.function = 'gauss'
        iraf.noao.digiphot.daophot.daopars.varord = 0
        iraf.noao.digiphot.daophot.daopars.fitsky = 'yes'
        iraf.noao.digiphot.daophot.daopars.recenter = 'yes'

        if _instrume == 'sofi':
            iraf.noao.digiphot.daophot.datapars.exposure = 'DIT'
        else:
            iraf.noao.digiphot.daophot.datapars.exposure = 'EXPTIME'

        zero = []
        fil = open(re.sub('.fits', '.ph', img), 'w')
        fil.write(str(_instrume) + ' ' + str(_date) + '\n')
        fil.write('*** ' + _object + ' ' + str(len(rastd0)) + '\n')
        if _field == 'landolt':
            fil.write('%6.6s\t%6.6s\t%6.6s\t%6.6s\t%6.6s\n' %
                      (str(1), str(1), str(1), str(1), str(1)))  # exptime
            fil.write('%6.6s\t%6.6s\t%6.6s\t%6.6s\t%6.6s\n' % (
                str(airmass0['U']), str(airmass0['B']), str(airmass0['V']), str(airmass0['R']), str(airmass0['I'])))
        elif _field == 'sloan':
            fil.write('%6.6s\t%6.6s\t%6.6s\t%6.6s\t%6.6s\n' %
                      (str(1), str(1), str(1), str(1), str(1)))  # exptime
            fil.write('%6.6s\t%6.6s\t%6.6s\t%6.6s\t%6.6s\n' % (
                str(airmass0['u']), str(airmass0['g']), str(airmass0['r']), str(airmass0['i']), str(airmass0['z'])))
        elif _field == '2mass':
            fil.write('%6.6s\t%6.6s\t%6.6s\n' %
                      (str(1), str(1), str(1)))  # exptime
            fil.write('%6.6s\t%6.6s\t%6.6s\n' % (
                str(airmass0['J']), str(airmass0['H']), str(airmass0['K'])))
        for i in range(0, len(pos1)):
            gg = open('tmp.one', 'w')
            gg.write(str(xsex[pos1[i]]) + ' ' + str(ysex[pos1[i]]) + '\n')
            gg.close()
            try:
                #############  magnitude  ##########
                #    daophot magnitude
                #    mag  =   mag(daophot) - K * airmass - log10(exptime)
                #    exptime for sofi is DIT
                #    exptime for efosc is EXPTIME
                #
                phot = iraf.noao.digiphot.daophot.phot(image=img, output='', coords='tmp.one', verify='no',
                                                       interactive='no', Stdout=1)
                mag0 = float(string.split(phot[0])[4])
                #                mag=mag0-kk[filters[_filter]]*float(_airmass)
                mag = mag0 - 2.5 * \
                    math.log10(float(_exptime)) - \
                    kk[filters[_filter]] * float(_airmass)
            except:
                mag0 = 999
                mag = 999

            jj = pos0[i]
            # in the ph file we stroe the insturmental magnitude not corrected
            # by K * airmass
            # instrumental mangitude of std in pos0[i]
            fileph['m' + filters[_filter]][jj] = mag
            if _field == 'landolt':
                stringastandard = '%12.12s\t%7.7s\t%7.7s\t%7.7s\t%7.7s\t%7.7s' % (
                    idstd0[jj], fileph['V'][jj], fileph['BV'][jj], fileph['UB'][jj], fileph['VR'][jj], fileph['RI'][jj])
                fil.write('%7.7s\t%7.7s\t%7.7s\t%7.7s\t%7.7s\t%60.60s\n'
                          % (str(fileph['mU'][jj]), str(fileph['mB'][jj]), str(fileph['mV'][jj]), str(fileph['mR'][jj]),
                             str(fileph['mI'][jj]), str(stringastandard)))
            elif _field == 'sloan':
                stringastandard = '%12.12s\t%7.7s\t%7.7s\t%7.7s\t%7.7s\t%7.7s' % (
                    idstd0[jj], fileph['r'][jj], fileph['gr'][jj], fileph['ug'][jj], fileph['ri'][jj], fileph['iz'][jj])
                fil.write('%7.7s\t%7.7s\t%7.7s\t%7.7s\t%7.7s\t%60.60s\n' % (
                    str(fileph['mu'][jj]), str(fileph['mg'][jj]), str(
                        fileph['mr'][jj]), str(fileph['mi'][jj]),
                    str(fileph['mz'][jj]), str(stringastandard)))
            elif _field == '2mass':
                stringastandard = '%12.12s\t%7.7s\t%7.7s\t%7.7s' % (
                    idstd0[jj], fileph['J'][jj], fileph['JH'][jj], fileph['HK'][jj])
                fil.write('%7.7s\t%7.7s\t%7.7s\t%60.60s\n' % (str(fileph['mJ'][jj]), str(fileph['mH'][jj]),
                                                              str(fileph['mK'][jj]), str(stringastandard)))
            zero.append(
                float(float(magstd0[filters[_filter]][jj])) - float(mag))

        fil.close()

        for col in colorvec:
            col0 = compress(
                (array(xstd, float) < readkey3(hdr, 'naxis1')) & (array(xstd, float) > 0) &
                (array(ystd, float) > 0) & (array(ystd, float) < readkey3(hdr, 'naxis2')), stdcoo[col[0]])
            col1 = compress(
                (array(xstd, float) < readkey3(hdr, 'naxis1')) & (array(xstd, float) > 0) &
                (array(ystd, float) > 0) & (array(ystd, float) < readkey3(hdr, 'naxis2')), stdcoo[col[1]])
            colstd0 = array(col0, float) - array(col1, float)
            ################## sex  ######################
            colore = []
            for i in range(0, len(pos1)):
                colore.append(colstd0[pos0[i]])
            colore = compress(abs(array(zero)) < 50, array(colore))
            zero = compress(abs(array(zero)) < 50, array(zero))
            if len(colore) == 0:
                b, a, RR = 9999, 9999, 9999
                result = ''
                print('no calibration, ' + _filter + ' ' + _field)
                if _filter not in zz:
                    ntt.util.updateheader(
                        img, 0, {'PHOTZP': [9999., 'MAG=-2.5*log(data)+PHOTZP']})
                    ntt.util.updateheader(
                        img, 0, {'PHOTZPER': [9999., 'error in PHOTZP']})
                    ntt.util.updateheader(
                        img, 0, {'FLUXCAL': ['UNCALIBRATED', 'Certifies the validity of PHOTZP']})
                else:
                    ntt.util.updateheader(
                        img, 0, {'PHOTZP': [zz[_filter], 'MAG=-2.5*log(data)+PHOTZP']})
                    #ntt.util.updateheader(img, 0, {'PHOTZPER': [2.0, 'error in PHOTZP']})
                    ntt.util.updateheader(
                        img, 0, {'PHOTZPER': [999, 'error in PHOTZP']})
                    ntt.util.updateheader(
                        img, 0, {'FLUXCAL': ['ABSOLUTE', 'Certifies the validity of PHOTZP']})
            elif len(colore) > 1:
                if not _interactive:
                    ntt.util.updateheader(
                        img, 0, {'PHOTZP': [mean(zero), 'MAG=-2.5*log(data)+PHOTZP']})
                    ntt.util.updateheader(
                        img, 0, {'PHOTZPER': [std(zero) / len(zero), 'error in PHOTZP']})
                    ntt.util.updateheader(
                        img, 0, {'FLUXCAL': ['ABSOLUTE', 'Certifies the validity of PHOTZP']})
                    if _field == 'landolt':
                        a, b, RR = ntt.efoscastrodef.linreg(colore, zero)
                    elif _field == '2mass':
                        b, a, RR = mean(zero), 0, 0
                    elif _field == 'sloan':
                        b, a, RR = mean(zero), 0, 0
                    xx = [min(array(colore)), max(array(colore))]
                    yy = ntt.efoscastrodef.pval(array(xx), [b, a])
                    result[filters[_filter] + col] = [a, b, RR]
                else:
                    print('do zeropoint interactively')
                    import pylab as plt

                    a, sa, b, sb = ntt.efoscastrodef.fitcol(
                        colore, zero, _filter, col, 0.0)
                    RR = 0  # 1 - residual/meanerror
                    ntt.util.updateheader(
                        img, 0, {'PHOTZP': [a, 'MAG=-2.5*log(data)+PHOTZP']})
                    ntt.util.updateheader(
                        img, 0, {'PHOTZPER': [sa, 'error in PHOTZP']})
                    ntt.util.updateheader(
                        img, 0, {'FLUXCAL': ['ABSOLUTE', 'Certifies the validity of PHOTZP']})
                    xx = [min(array(colore)), max(array(colore))]
                    yy = ntt.efoscastrodef.pval(array(xx), [b, a])
                    result[filters[_filter] + col] = [a, b, RR]
            else:
                ntt.util.updateheader(
                    img, 0, {'PHOTZP': [zero[0], 'MAG=-2.5*log(data)+PHOTZP ']})
                ntt.util.updateheader(
                    img, 0, {'PHOTZPER': [1.0, 'error in PHOTZP']})
                ntt.util.updateheader(
                    img, 0, {'FLUXCAL': ['ABSOLUTE', 'Certifies the validity of PHOTZP']})
                b, a, RR = zero[0], 0, 0
                xx = colore
                yy = zero
                result[filters[_filter] + col] = [0, zero[0], 0]

            if verbose and len(colore) > 0:
                from pylab import plot, show, ion, xlim, ylim, legend, setp, gca, draw, clf

                clf()
                ion()
                _label = '%2s   %2s   %5.5s  %5.5s' % (
                    filters[_filter], col, str(b), str(a))
                plot(xx, yy, '-', label='fit')
                plot(colore, zero, 'o', label=_label)
                if len(colore) > 1:
                    xlim(min(array(colore)) - .2, max(array(colore)) + .2)
                    ylim(min(array(zero)) - 1, max(array(zero) + 1))
                legend(numpoints=1, markerscale=1.5)
                show()
    else:
        if _filter not in zz:
            ntt.util.updateheader(
                img, 0, {'PHOTZP': [9999., 'MAG=-2.5*log(data)+PHOTZP']})
            ntt.util.updateheader(
                img, 0, {'PHOTZPER': [9999., 'error in PHOTZP']})
            ntt.util.updateheader(
                img, 0, {'FLUXCAL': ['UNCALIBRATED', 'Certifies the validity of PHOTZP']})
            print('no calibration, ' + _filter + ' ' + _field)
        else:
            ntt.util.updateheader(
                img, 0, {'PHOTZP': [zz[_filter], 'MAG=-2.5*log(data)+PHOTZP']})
            ntt.util.updateheader(
                img, 0, {'PHOTZPER': [2.0, 'error in PHOTZP']})
            ntt.util.updateheader(
                img, 0, {'FLUXCAL': ['ABSOLUTE', 'Certifies the validity of PHOTZP']})
            print('calibration from default zero point')
        result = ''

    if _field == 'filter':
        if _filter in ['J', 'Js', 'H', 'Ks', 'K']:
            _field = '2mass'
        elif _filter in ['U640', 'B639', 'V641', 'R642', 'i705']:
            _field = 'landolt'
        elif _filter in ['g782', 'r784', 'z623']:
            _field = 'sloan'
    if _field == 'landolt':
        ntt.util.updateheader(
            img, 0, {'PHOTSYS': ['VEGA', 'Photometric system VEGA or AB']})
    elif _field == '2mass':
        ntt.util.updateheader(
            img, 0, {'PHOTSYS': ['VEGA', 'Photometric system VEGA or AB']})
    elif _field == 'sloan':
        ntt.util.updateheader(
            img, 0, {'PHOTSYS': ['AB', 'Photometric system VEGA or AB']})
    else:
        ntt.util.updateheader(
            img, 0, {'PHOTSYS': ['NULL', 'Photometric system VEGA or AB']})
    delete('tmp.*,detection*')
    return result


#####################################################################
def pval(_xx, p):
    # print "LOGX:: Entering `pval` method/function in %(__file__)s" %
    # globals()
    _y = +p[0] + p[1] * _xx
    return _y


def crossmatch(_ra0, _dec0, _ra1, _dec1, tollerance):  # degree,degree,degree,degree,arcsec
    # print "LOGX:: Entering `crossmatch` method/function in %(__file__)s" %
    # globals()
    from numpy import pi, cos, sin, array, argmin, min, arccos

    scal = pi / 180.
    distvec = []
    pos0 = []
    pos1 = []
    for jj in range(0, len(_ra0)):
        distance = arccos(sin(array(_dec1) * scal) * sin(_dec0[jj] * scal) + cos(array(_dec1) * scal) *
                          cos(_dec0[jj] * scal) * cos((array(_ra1) - _ra0[jj]) * scal))
        if min(distance) <= tollerance * pi / (180 * 3600):
            distvec.append(min(distance))
            pos0.append(jj)
            pos1.append(argmin(distance))
    return distvec, pos0, pos1


###################################################################

def onkeypress(event):
    # print "LOGX:: Entering `onkeypress` method/function in %(__file__)s" %
    # globals()
    global idd, _col, _dmag, testo, lines, pol, sss, f, fixcol, sigmaa, sigmab
    import numpy as np
    import pylab as plt

    xdata, ydata = event.xdata, event.ydata
    dist = np.sqrt((xdata - _col) ** 2 + (ydata - _dmag) ** 2)
    ii = np.argmin(dist)
    if event.key == 'd':
        __col, __dmag = _col.tolist(), _dmag.tolist()
        plt.plot(_col[ii], _dmag[ii], 'xk', ms=25)
        del __col[ii], __dmag[ii]
        _col, _dmag = np.array(__col), np.array(__dmag)

    idd = range(len(_col))
    _fixcol = fixcol
    if len(_col[idd]) == 1:
        _fixcol = 0.0
    else:
        _fixcol = fixcol

    if _fixcol == '':
        pol = np.polyfit(_col, _dmag, 1, full=True)
        xx = [min(_col) - .1, max(_col) + .1]
        yy = np.polyval(pol[0], xx)
    else:
        pol = [[0, 0], []]
        pol[0][1] = np.mean(np.array(_dmag) - np.array(_col) * float(_fixcol))
        pol[0][0] = _fixcol
        xx = [min(_col) - .1, max(_col) + .1]
        yy = np.polyval(pol[0], xx)
        pol[1] = np.std(
            abs(pol[0][1] - (np.array(_dmag[idd]) - np.array(_col[idd]) * float(_fixcol))))
    #        pol[1]=std(array(_dmag)-array(_col)*float(_fixcol))

    try:
        sigmae = np.sqrt(pol[1] / (len(idd) - 2))
    except:
        sigmae = 0
    sigmaa = sigmae * np.sqrt(1. / len(idd) + (
        np.mean(_col[idd]) ** 2) / np.sum((_col[idd] - np.mean(_col[idd])) ** 2))
    sigmab = sigmae * \
        np.sqrt(1 / np.sum((_col[idd] - np.mean(_col[idd])) ** 2))
    if _fixcol != '':
        sigmab = 0.0

    lines.pop(0).remove()
    lines = plt.plot(xx, yy, 'r-')
    plt.ylim(min(_dmag) - .2, max(_dmag) + .2)
    plt.xlabel(sss)
    plt.title(f)

    try:
        plt.setp(testo, text='%5.3f + %s* %5.3f [%4.3f  %4.3f]' %
                             (pol[0][1], sss, pol[0][0], sigmaa, sigmab))
    except:
        plt.setp(testo, text='%5.3f + %s* %5.3f [%4.3f  %4.3f]' %
                             (pol[0][1], sss, pol[0][0], sigmaa, sigmab))


####################################################################

def onclick(event):
    # print "LOGX:: Entering `onclick` method/function in %(__file__)s" %
    # globals()
    global idd, _col, _dmag, testo, lines, pol, sss, f, fixcol, sigmaa, sigmab
    import numpy as np
    import pylab as plt

    xdata, ydata = event.xdata, event.ydata

    if xdata and ydata:
        dist = np.sqrt((xdata - _col) ** 2 + (ydata - _dmag) ** 2)
        ii = np.argmin(dist)
        if event.button == 2:
            if ii not in idd:
                idd.append(ii)
        if event.button == 1:
            if ii in idd:
                idd.remove(ii)

    nonincl = []
    for i in range(len(_col)):
        if i not in idd:
            nonincl.append(i)

    _fixcol = fixcol
    if len(_col[idd]) == 1:
        _fixcol = 0.0
    else:
        _fixcol = fixcol

    if _fixcol == '':
        pol = np.polyfit(_col[idd], _dmag[idd], 1, full=True)
        xx = [min(_col) - .1, max(_col) + .1]
        yy = np.polyval(pol[0], xx)
    else:
        pol = [[0, 0], [0]]
        pol[0][1] = np.mean(_dmag[idd] - _col[idd] * float(_fixcol))
        pol[0][0] = _fixcol
        xx = [min(_col) - .1, max(_col) + .1]
        yy = np.polyval(pol[0], xx)
        pol[1] = np.std(np.abs(
            pol[0][1] - (np.array(_dmag[idd]) - np.array(_col[idd]) * float(_fixcol))))

    try:
        sigmae = np.sqrt(pol[1] / (len(idd) - 2))
    except:
        sigmae = 0

    sigmaa = sigmae * np.sqrt(1. / len(idd) + (
        np.mean(_col[idd]) ** 2) / np.sum((_col[idd] - np.mean(_col[idd])) ** 2))
    sigmab = sigmae * \
        np.sqrt(1 / np.sum((_col[idd] - np.mean(_col[idd])) ** 2))
    if _fixcol != '':
        sigmab = 0.0

    plt.plot(_col, _dmag, 'ok')
    plt.plot(_col[nonincl], _dmag[nonincl], 'ow')
    lines.pop(0).remove()
    lines = plt.plot(xx, yy, 'r-')
    plt.ylim(min(_dmag) - .2, max(_dmag) + .2)
    plt.xlabel(sss)
    plt.title(f)
    try:
        plt.setp(testo, text='%5.3f + %s* %5.3f [%4.3f  %4.3f]' %
                             (pol[0][1], sss, pol[0][0], sigmaa, sigmab))
    except:
        plt.setp(testo, text='%5.3f + %s* %5.3f [%4.3f  %4.3f]' %
                             (pol[0][1], sss, pol[0][0], sigmaa, sigmab))


#################################################

def fitcol(col, dmag, band, color, fissa=''):
    # print "LOGX:: Entering `fitcol` method/function in %(__file__)s" %
    # globals()
    global idd, _col, _dmag, testo, lines, pol, sss, f, fixcol, sigmaa, sigmab
    import pylab as plt
    import numpy as np

    plt.ion()
    fig = plt.figure()
    _dmag = dmag[:]
    _col = col[:]
    sss = band
    f = color
    fixcol = fissa
    _col = np.array(_col)
    _dmag = np.array(_dmag)
    idd = range(len(_col))
    plt.plot(_col, _dmag, 'ok')

    _fixcol = fixcol
    if len(_col[idd]) == 1:
        _fixcol = 0.0
    else:
        _fixcol = fixcol

    if _fixcol == '':
        pol = np.polyfit(_col[idd], _dmag[idd], 1, full=True)
        xx = [min(_col) - .1, max(_col) + .1]
        yy = np.polyval(pol[0], xx)
    else:
        pol = [[0, 0], []]
        pol[0][1] = np.mean(np.array(_dmag)[idd] -
                            np.array(_col)[idd] * fixcol)
        pol[0][0] = fixcol
        xx = [min(_col) - .1, max(_col) + .1]
        yy = np.polyval(pol[0], xx)
        pol[1] = np.std(
            np.abs(pol[0][1] - (np.array(_dmag) - np.array(_col) * float(_fixcol))))
    if len(_col) <= 2:
        pol = [pol[0], [0]]

    lines = plt.plot(xx, yy, 'r-')
    plt.ylim(min(_dmag) - .2, max(_dmag) + .2)
    plt.xlim(min(xx), max(xx))
    plt.xlabel(sss)
    plt.title(f)

    try:
        sigmae = np.sqrt(pol[1] / (len(idd) - 2))
    except:
        sigmae = 0
    sigmaa = sigmae * np.sqrt(1. / len(idd) + (
        np.mean(_col[idd]) ** 2) / np.sum((_col[idd] - np.mean(_col[idd])) ** 2))
    sigmab = sigmae * np.sqrt(1 / sum((_col[idd] - np.mean(_col[idd])) ** 2))

    try:
        testo = plt.figtext(.2, .85, '%5.3f + %s* %5.3f [%4.3f  %4.3f]' %
                            (pol[0][1], sss, pol[0][0], sigmaa, sigmab))
    except:
        testo = plt.figtext(.2, .85, '%5.3f + %s* %5.3f [%4.3f  %4.3f]' %
                            (pol[0][1], sss, pol[0][0], sigmaa, sigmab))

    kid = fig.canvas.mpl_connect('key_press_event', onkeypress)
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.draw()
    raw_input(
        'left-click mark bad, right-click unmark, <d> remove. Return to exit ...')
    plt.close()
    if _fixcol == '':
        sigmaa, sigmab = sigmaa[0], sigmab[0]
    return pol[0][1], sigmaa, pol[0][0], sigmab


####################################################################
def linreg(X, Y):
    # print "LOGX:: Entering `linreg` method/function in %(__file__)s" %
    # globals()
    from math import sqrt

    """
    Summary
        Linear regression of y = ax + b
    Usage
        real, real, real = linreg(list, list)
    Returns coefficients to the regression line "y=ax+b" from x[] and y[], and R^2 Value
    """
    if len(X) != len(Y):
        raise ValueError('unequal length')
    N = len(X)
    Sx = Sy = Sxx = Syy = Sxy = 0.0
    for x, y in map(None, X, Y):
        Sx = Sx + x
        Sy = Sy + y
        Sxx = Sxx + x * x
        Syy = Syy + y * y
        Sxy = Sxy + x * y
    det = Sxx * N - Sx * Sx
    a, b = (Sxy * N - Sy * Sx) / det, (Sxx * Sy - Sx * Sxy) / det
    meanerror = residual = 0.0
    for x, y in map(None, X, Y):
        meanerror = meanerror + (y - Sy / N) ** 2
        residual = residual + (y - a * x - b) ** 2
    RR = 1 - residual / meanerror
    ss = residual / (N - 2)
    Var_a, Var_b = ss * N / det, ss * Sxx / det
    return a, b, RR


##########################################################################
##########################################
def querysloan(ra1, dec1, radius, mr1, mr2):
    # print "LOGX:: Entering `querysloan` method/function in %(__file__)s" %
    # globals()
    import ntt
    from ntt import sqlcl
    import os
    import sys
    import shutil
    import string

    righe = sqlcl.query('select P.objID, P.ra, P.dec , P.u, P.g, P.r, P.i, P.z, P.type ' +
                        'from PhotoPrimary as P , dbo.fGetNearbyObjEq(' + str(ra1) + ', ' + str(dec1) + ', ' +
                        str(radius) + ' ) N' + ' where P.objID = N.objID').readlines()
    _id, _ra, _dec, _u, _g, _r, _i, _z, _type = [], [], [], [], [], [], [], [], []
    for i in righe[1:]:
        if len(string.split(i, ',')) == 9:
            _id0, _ra0, _dec0, _u0, _g0, _r0, _i0, _z0, _type0 = string.split(
                i, ',')
            if mr1 and mr2:
                if mr1 <= float(_r0) <= mr2:
                    _id.append(_id0)
                    _ra.append(float(_ra0))
                    _dec.append(float(_dec0))
                    try:
                        _u.append(float(_u0))
                    except:
                        _u.append(float(9999.))
                    try:
                        _g.append(float(_g0))
                    except:
                        _g.append(float(9999.))
                    try:
                        _r.append(float(_r0))
                    except:
                        _r.append(float(9999.))
                    try:
                        _i.append(float(_i0))
                    except:
                        _i.append(float(9999.))
                    try:
                        _z.append(float(_z0))
                    except:
                        _z.append(float(9999.))
                    _type.append(_type0[:-1])
            else:
                _id.append(_id0)
                _ra.append(float(_ra0))
                _dec.append(float(_dec0))
                try:
                    _u.append(float(_u0))
                except:
                    _u.append(float(9999.))
                try:
                    _g.append(float(_g0))
                except:
                    _g.append(float(9999.))
                try:
                    _r.append(float(_r0))
                except:
                    _r.append(float(9999.))
                try:
                    _i.append(float(_i0))
                except:
                    _i.append(float(9999.))
                try:
                    _z.append(float(_z0))
                except:
                    _z.append(float(9999.))
                _type.append(_type0[:-1])
    return _id, _ra, _dec, _u, _g, _r, _i, _z, _type


#################################################################
def sloan2file(ra1, dec1, radius, magsel0, magsel1, _output):
    # print "LOGX:: Entering `sloan2file` method/function in %(__file__)s" %
    # globals()
    import os
    import re
    import string
    import ntt
    from numpy import array, compress

    ra2 = float(ra1) * 15  # degree
    _ids, _ras, _decs, _us, _gs, _rs, _is, _zs, _type = ntt.efoscastrodef.querysloan(ra1, dec1, float(radius), magsel0,
                                                                                     magsel1)
    headersloan = '# BEGIN CATALOG HEADER \n\
# catdb sqlcl.query \n# catname dss1@cadc \n# nquery 4 \n#     ra ' + str(ra2) + ' hours \n\
#     dec ' + str(dec1) + ' degrees \n#     radius ' + str(radius) + '  minutes \n#     qsystem J2000.0 INDEF \n# type btext \n\
# nheader 1 \n#     csystem J2000.0 \n# nfields 8 \n#       ra   1 0 d degrees %10.5f \n\
#     dec    2 0 d degrees %10.5f \n#     id   3 0 c NDEF %11s \n\
#     u    4 0 d degrees %6.3f \n#     g    5 0 d degrees %6.3f \n#     r   6 0 r INDEF %6.3f \n\
#     i   7 0 r INDEF %6.3f \n#     z   8 0 r INDEF %6.3f \n'
    ff = open(_output, 'w')
    ff.write(headersloan)
    ff.close()
    if _ids:
        for i in range(0, len(_type)):
            try:
                _type[i] = float(_type[i])
            except:
                _type[i] = 9999
        _ras = array(_ras)
        _decs = array(_decs)
        _type = array(_type)
        _ids = array(_ids)
        _rass = compress(_type == 6, _ras)
        _decss = compress(_type == 6, _decs)
        _idss = compress(_type == 6, _ids)
        _gss = compress(_type == 6, _gs)
        _rss = compress(_type == 6, _rs)
        _iss = compress(_type == 6, _is)
        _zss = compress(_type == 6, _zs)
        _uss = compress(_type == 6, _us)
        ff = open(_output, 'a')
        for i in range(0, len(_idss)):
            ff.write('%12.12s\t%12.12s\t%s\t%8.8s\t%8.8s\t%8.8s\t%8.8s\t%8.8s\n' %
                     (str(_rass[i]), str(_decss[i]), str(_idss[i][-6:]), str(_uss[i]), str(_gss[i]), str(_rss[i]),
                      str(_iss[i]), str(_zss[i])))
        ff.close()
    else:
        _idss = ''
    return _idss


#######################################################################

def transformsloanlandolt(stdcoo):
    # print "LOGX:: Entering `transformsloanlandolt` method/function in
    # %(__file__)s" % globals()
    from numpy import array
    # lupton 2005  ###### BVRI
    # Jester et al. (2005) ### U
    if 'u' and 'g' and 'r' in stdcoo:
        stdcoo['B'] = array(stdcoo['g'], float) + 0.3130 * (
            array(stdcoo['g'], float) - array(stdcoo['r'], float)) + 0.2271
        stdcoo['U'] = array(stdcoo['B'], float) + 0.78 * \
            (array(stdcoo['u'], float) - array(stdcoo['g'], float)) - 0.88
    if 'g' and 'r' in stdcoo:
        if 'B' not in stdcoo:
            stdcoo['B'] = array(stdcoo['g'], float) + 0.3130 * (
                array(stdcoo['g'], float) - array(stdcoo['r'], float)) + 0.2271
        if 'V' not in stdcoo:
            stdcoo['V'] = array(stdcoo['g'], float) - 0.5784 * (
                array(stdcoo['g'], float) - array(stdcoo['r'], float)) - 0.0038
    if 'r' and 'i' in stdcoo:
        if 'R' not in stdcoo:
            stdcoo['R'] = array(stdcoo['r'], float) - 0.2936 * (
                array(stdcoo['r'], float) - array(stdcoo['i'], float)) - 0.1439
        if 'I' not in stdcoo:
            stdcoo['I'] = array(stdcoo['r'], float) - 1.2444 * (
                array(stdcoo['r'], float) - array(stdcoo['i'], float)) - 0.3820
    if 'g' and 'r' in stdcoo:
        if 'R' not in stdcoo:
            stdcoo['R'] = array(stdcoo['r'], float) - 0.1837 * (
                array(stdcoo['g'], float) - array(stdcoo['r'], float)) - 0.0971
    if 'i' and 'z' in stdcoo:
        if 'I' not in stdcoo:
            stdcoo['I'] = array(stdcoo['i'], float) - 0.3780 * (
                array(stdcoo['i'], float) - array(stdcoo['z'], float)) - 0.3974
    for fil in 'UBVRI':
        if fil not in stdcoo:
            stdcoo[fil] = array(stdcoo['r'], float) - array(stdcoo['r'], float)
    return stdcoo


#####################################################

def transformlandoltsloan(stdcoo):
    # print "LOGX:: Entering `transformlandoltsloan` method/function in %(__file__)s" % globals()
    #  jordi et al 2006
    from numpy import array

    if 'B' and 'V' in stdcoo:
        if 'g' not in stdcoo:
            stdcoo['g'] = array(stdcoo['V'], float) + 0.630 * (
                array(stdcoo['B'], float) - array(stdcoo['V'], float)) - 0.124
    if 'V' and 'R' in stdcoo:
        if 'r' not in stdcoo:
            VR = (array(stdcoo['V'], float) - array(stdcoo['R'], float))
            a = array(stdcoo['R'], float) + 0.267 * VR + 0.088  # V-R < 0.93
            b = array(stdcoo['R'], float) + 0.77 * VR - 0.37  # V-R > 0.93
            stdcoo['r'] = [(a[i], b[i])[VR[i] <= 0.93]
                           for i in range(0, len(VR))]
    if 'R' and 'I' in stdcoo:
        if 'i' not in stdcoo:
            stdcoo['i'] = array(stdcoo['I'], float) - 0.247 * (
                array(stdcoo['R'], float) - array(stdcoo['I'], float)) + 0.329
    if 'V' and 'R' and 'I' in stdcoo:
        if 'r' not in stdcoo:
            VR = (array(stdcoo['V'], float) - array(stdcoo['R'], float))
            a = array(stdcoo['R'], float) + 0.267 * VR + 0.088  # V-R < 0.93
            b = array(stdcoo['R'], float) + 0.77 * VR - 0.37  # V-R > 0.93
            stdcoo['r'] = [(a[i], b[i])[VR[i] <= 0.93]
                           for i in range(0, len(VR))]
        if 'z' not in stdcoo:
            stdcoo['z'] = array(stdcoo['r'], float) - 1.584 * (
                array(stdcoo['R'], float) - array(stdcoo['I'], float)) + (0.386)
    for fil in 'ugriz':
        if fil not in stdcoo:
            stdcoo[fil] = array(stdcoo['V'], float) - array(stdcoo['V'], float)
    return stdcoo


####################################################
# MagLim = ZP - 2.5*log10(5 * sqrt(pi*rcore^2) * skynoise / exptime) - apcor -percorr
# apcor =0
# percorr =0
# MagLim = ZP - 2.5*log10(5 * sqrt(pi*FWHM^2) * skynoise / exptime)
# M_limiting = M + 2.5 * log10 (SNR/3).
# where M = Magnitude of star and SNR is the SNR of the star.
# This gives you a 3-sigma limiting magnitude.
##############    sextractor   ##################

def sextractor(img):
    # print "LOGX:: Entering `sextractor` method/function in %(__file__)s" %
    # globals()
    from ntt import defsex, delete
    import os
    from pyraf import iraf
    from iraf import proto
    from numpy import compress, array, asarray

    namesex = defsex('default.sex')
    os.system('sex ' + img + ' -c ' + namesex + ' > _logsex')
    delete(namesex)
    delete('_logsex')
    xpix = iraf.proto.fields('detections.cat', fields='2', Stdout=1)
    ypix = iraf.proto.fields('detections.cat', fields='3', Stdout=1)
    cm = iraf.proto.fields('detections.cat', fields='4', Stdout=1)
    cl = iraf.proto.fields('detections.cat', fields='7', Stdout=1)
    # flux radius in sextractor PHOT_FLUXFRAC    0.5
    fw = iraf.proto.fields('detections.cat', fields='8', Stdout=1)
    ell = iraf.proto.fields('detections.cat', fields='9', Stdout=1)
    bkg = iraf.proto.fields('detections.cat', fields='10', Stdout=1)

    cl = compress((array(xpix) != ''), array(cl, float))
    cm = compress((array(xpix) != ''), array(cm, float))
    fw = compress((array(xpix) != ''), array(fw, float))
    ell = compress((array(xpix) != ''), array(ell, float))
    bkg = compress((array(xpix) != ''), array(bkg, float))
    ypix = compress((array(xpix) != ''), array(ypix, float))
    xpix = compress((array(xpix) != ''), array(xpix, float))
    try:
        ww = asarray([i for i in range(len(xpix)) if (
            (xpix[i] < 960) or (ypix[i] < 980))])
        cl, cm, fw, ell, xpix, ypix, bkg = cl[ww], cm[
            ww], fw[ww], ell[ww], xpix[ww], ypix[ww], bkg[ww]

        ww = asarray([i for i in range(len(xpix)) if (
            (xpix[i] > 20) or (ypix[i] < 980))])
        cl, cm, fw, ell, xpix, ypix, bkg = cl[ww], cm[
            ww], fw[ww], ell[ww], xpix[ww], ypix[ww], bkg[ww]

        ww = asarray([i for i in range(len(xpix)) if (xpix[i] > 3)])
        cl, cm, fw, ell, xpix, ypix, bkg = cl[ww], cm[
            ww], fw[ww], ell[ww], xpix[ww], ypix[ww], bkg[ww]

        cl = compress((array(fw) <= 15) & (array(fw) >= -2), array(cl))
        cm = compress((array(fw) <= 15) & (array(fw) >= -2), array(cm))
        xpix = compress((array(fw) <= 15) & (array(fw) >= -2), array(xpix))
        ypix = compress((array(fw) <= 15) & (array(fw) >= -2), array(ypix))
        ell = compress((array(fw) <= 15) & (array(fw) >= -2), array(ell))
        bkg = compress((array(fw) <= 15) & (array(fw) >= -2), array(bkg))
        fw = compress((array(fw) <= 15) & (array(fw) >= -2), array(fw))
    except:
        xpix, ypix, fw, cl, cm, ell = [], [], [], [], [], []
        print('\n### ERROR Filtering the sextractor detections, please check that sextractor is working ......')
    delete('detections.cat')
    return xpix, ypix, fw, cl, cm, ell, bkg


##########################################################################
#####
def querycatalogue(catalogue, img, method='iraf'):

    # print "LOGX:: Entering `querycatalogue` method/function in %(__file__)s" % globals()
    # print "LOGX:: catalogue %(catalogue)s" % locals()
    # print "LOGX:: img %(img)s" % locals()
    # print "LOGX:: method %(method)s" % locals()
    from pyraf import iraf
    from numpy import array, compress
    import string
    import re
    import sys
    import os
    from ntt import delete
    from ntt.util import readhdr, readkey3, delete
    import ntt

    hdr = readhdr(img)
    _ra = readkey3(hdr, 'RA')
    _dec = readkey3(hdr, 'DEC')
    iraf.imcoords(_doprint=0, Stdout=0)
    iraf.astcat(_doprint=0, Stdout=0)
    toforget = ['imcoords', 'astcat', 'tv']
    for t in toforget:
        try:
            iraf.unlearn(t)
        except:
            pass

    iraf.noao.astcat.aregpars.rcrauni = ''
    iraf.noao.astcat.aregpars.rcdecuni = ''
    iraf.noao.astcat.catdb = ntt.__path__[0] + '/standard/cat/catalogue.dat'
    iraf.noao.astcat.aregpars.rcra = _ra / 15
    iraf.noao.astcat.aregpars.rcdec = _dec
    iraf.noao.astcat.aregpars.rrawidt = 15.
    iraf.noao.astcat.aregpars.rdecwid = 15.
    delete('tmp.catalogue')
    delete('tmp.catalogue.pix')

    if method == 'iraf':
        if catalogue == 'usnoa2':
            lll = iraf.noao.astcat.agetcat(
                'pars', 'STDOUT', catalog='usno2@noao', verbose='no', Stdout=1)
        elif catalogue == 'usnob1':
            lll = iraf.noao.astcat.agetcat(
                'pars', 'STDOUT', catalog='usnob1@noao', verbose='no', Stdout=1)
        elif catalogue == '2mass':
            lll = iraf.noao.astcat.agetcat(
                'pars', 'STDOUT', catalog='twomass@irsa', verbose='yes', Stdout=1)
        # elif catalogue=='gsc1':
        # lll=iraf.noao.astcat.agetcat('pars','STDOUT',catalog='gsc1@cadc',verbose='no',Stdout=1)
        else:
            if os.path.isfile(ntt.__path__[0] + '/standard/cat/' + catalogue):
                ff = open(ntt.__path__[0] + '/standard/cat/' + catalogue)
                lll = ff.readlines()
                ff.close()
                for ii in range(0, len(lll)):
                    lll[ii] = re.sub('\n', '', lll[ii])
                print('catalogue from user')
            else:
                sys.exit('Error: catalogue ' + str(catalogue) +
                         'not in the list [usnob1,usnoa2,2mass]')
            ########

        # REMOVE COMMENT LINES
        lllNew = []
        for ll in lll:
            if ll[:2] != "\ " and ll[:2] != "| ":
                lllNew.append(ll)
        lll = lllNew

        # for i, v in enumerate(lll):
        #     print i, v

        # FIND LINE WITH nfield value
        indfield = [i for i in range(0, len(lll)) if 'nfields' in lll[i]]

        fields = int(lll[indfield[0]].split()[-1])

        stdcoo = {}
        column = {}
        for j in range(indfield[0] + 1, indfield[0] + fields + 1):
            if lll[j].split()[1] not in column:
                column[lll[j].split()[1]] = int(lll[j].split()[2])
            if lll[j].split()[1] not in stdcoo:
                stdcoo[lll[j].split()[1]] = []

        startIndex = lll.index('# END CATALOG HEADER') + 2
        for i in lll[startIndex:]:
            for j in stdcoo.keys():
                val = i.split()[column[j] - 1]
                if j in ["ra", "dec"]:
                    val = val.replace("d", ":").replace(
                        "h", ":").replace("m", ":").replace("s", "")
                stdcoo[j].append(val)

        colonne3 = str(int(column['ra'])) + ' ' + str(int(column['dec']))
        if catalogue in ['usnoa2', 'usnob1', '2mass', 'gsc1']:
            colonne4 = {'usnoa2': 'mag1', 'usnob1': 'R2mag',
                        '2mass': 'mag1', 'gsc1': 'mag'}
        else:
            for jj in column.keys():
                if jj in ['U', 'B', 'V', 'R', 'I', 'g', 'r', 'i', 'z']:
                    colonne4 = {catalogue: jj}
                    break

    elif method == 'vizir':
#       replace vizq with vizq2 to be consistent
#

        stdcoo = ntt.efoscastrodef.vizq2(_ra, _dec, catalogue, 10)

        lll = ['# END CATALOG HEADER', '#']
        for ff in range(0, len(stdcoo['ra'])):
            lll.append(str(stdcoo['ra'][ff]) + '  ' +
                       str(stdcoo['dec'][ff]) + '  ' + str(stdcoo['mag'][ff]))
        colonne4 = {'usnoa2': 'mag', 'usnob1': 'mag',
                    '2mass': 'mag', 'gsc1': 'mag'}
        colonne3 = ' 1   2 '
        column = {'ra': 1, 'dec': 2, 'r': 3}

    ddd2 = iraf.wcsctran('STDIN', 'STDOUT', img, Stdin=lll, Stdout=1, inwcs='world', units='degree degrees',
                         outwcs='logical', columns=colonne3, formats='%10.1f %10.1f')

    xx, yy = [], []
    for i in ddd2[ddd2.index('# END CATALOG HEADER') + 2:]:
        if i:
            xx.append(float(i.split()[column['ra'] - 1]))
            yy.append(float(i.split()[column['dec'] - 1]))
    #        colonne4={'usnoa2':'mag1','usnob1':'R2mag','2mass':'mag1','gsc1':'mag'}
    #######
    acoo1 = []
    apixx1, apixy1, am1, apix1 = [], [], [], []
    for i in range(0, len(stdcoo['ra'])):
        acoo1.append(str(stdcoo['ra'][i]) + ' ' + str(stdcoo['dec'][i]))
        apix1.append(str(xx[i]) + ' ' + str(yy[i]))
        am1.append(stdcoo[colonne4[catalogue]][i])

    if catalogue == '2mass':
        for jj in range(0, len(am1)):
            try:
                am1[jj] = float(re.sub('L', '', str(am1[jj])))
            except:
                am1[jj] = 999

    for key in stdcoo.keys():
        stdcoo[key] = compress((array(xx) < int(int(hdr['NAXIS1']) + 100)) & (array(xx) > -100) & (
            array(yy) < int(int(hdr['NAXIS2']) + 100)) & (array(yy) > -100), array(stdcoo[key]))
    stdcoo['coo'] = compress((array(xx) < int(int(hdr['NAXIS1']) + 100)) & (array(xx) > -100) & (
        array(yy) < int(int(hdr['NAXIS2']) + 100)) & (array(yy) > -100), array(acoo1))
    stdcoo['pix'] = compress((array(xx) < int(int(hdr['NAXIS1']) + 100)) & (array(xx) > -100) & (
        array(yy) < int(int(hdr['NAXIS2']) + 100)) & (array(yy) > -100), array(apix1))
    stdcoo['mag'] = compress((array(xx) < int(int(hdr['NAXIS1']) + 100)) & (array(xx) > -100) & (
        array(yy) < int(int(hdr['NAXIS2']) + 100)) & (array(yy) > -100), array(am1, float))
    stdcoo['x'] = compress((array(xx) < int(int(hdr['NAXIS1']) + 100)) & (array(xx) > -100) & (
        array(yy) < int(int(hdr['NAXIS2']) + 100)) & (array(yy) > -100), array(xx, float))
    stdcoo['y'] = compress((array(xx) < int(int(hdr['NAXIS1']) + 100)) & (array(xx) > -100) & (
        array(yy) < int(int(hdr['NAXIS2']) + 100)) & (array(yy) > -100), array(yy, float))
    return stdcoo


##########################################################################
def efoscastrometry2(lista, catalogue, _interactive, number, sexvec, catvec, guess=False, fitgeo='xyscale',
                     tollerance1=100, tollerance2=30, _update='yes', imex=False, nummin=4):
    # print "LOGX:: entering efoscastrometry2\n"
    import os
    import string
    import re
    import sys
    import numpy
    import math
    from numpy import array, compress, argsort, sort, asarray
    from numpy import round, mean, std, sqrt, median
    from numpy import argmin, isnan, abs, genfromtxt
    import ntt
    import time
    import datetime
    from ntt.efoscastrodef import wcsstart
    from ntt.util import delete, readhdr, readkey3, display_image, defsex
    from pyraf import iraf

    xpix, ypix, fw, cl, cm, ell, bkg = sexvec
    acoo1, apix1, am1 = catvec['coo'], catvec['pix'], catvec['mag']

    # catalogue
    iraf.noao(_doprint=0, Stdout=0)
    iraf.imcoords(_doprint=0, Stdout=0)
    iraf.tv(_doprint=0, Stdout=0)
    iraf.tv.rimexam.backgrou = 'yes'
    iraf.astcat(_doprint=0, Stdout=0)
    toforget = ['imcoords', 'astcat', 'tv']
    for t in toforget:
        try:
            iraf.unlearn(t)
        except:
            pass
    verbose = False
    if _interactive:
        verbose = True
    img = lista[0]
    hdr = readhdr(img)
    _instrume = readkey3(hdr, 'instrume')
    if _instrume == 'efosc':
        magsel0 = 7.0
        magsel1 = 21.
    elif _instrume == 'sofi':
        magsel0 = 7.0
        magsel1 = 21.
    #    if guess:     ntt.efoscastrodef.wcsstart(img)#,500,500)
    _CRPIX1 = readkey3(hdr, 'CRPIX1')
    _CRPIX2 = readkey3(hdr, 'CRPIX2')
    if verbose:
        display_image(img, 1, '', '', False)
        iraf.tvmark(1, 'STDIN', Stdin=list(apix1), mark="circle", number='yes', label='no', radii=20, nxoffse=5,
                    nyoffse=5, color=205, txsize=4)
        raw_input('mark catalogue ' + str(len(apix1)))
    else:
        #        ss=datetime.datetime.now()
        time.sleep(.7)
    answ = 'yes'
    magsel11 = magsel1
    mlim = 0
    while answ == 'yes':
        amcut1 = compress((array(am1) > magsel0) &
                          (array(am1) < magsel11), am1)
        if len(amcut1) <= number:
            answ = 'no'
            magsel11 = magsel1 + mlim + .5
        else:
            mlim = mlim - .5
            magsel11 = magsel1 + mlim

    amcut = compress((array(am1) > magsel0) & (array(am1) < magsel11), am1)
    apixcut = compress((array(am1) > magsel0) & (
        array(am1) < magsel11), apix1)  # usno x y  cut_list
    acoocut = compress((array(am1) > magsel0) & (
        array(am1) < magsel11), acoo1)  # usno ra dec  cut_list

    rausno = compress((array(am1) > magsel0) & (
        array(am1) < magsel11), array(catvec['ra'], float))
    decusno = compress((array(am1) > magsel0) & (
        array(am1) < magsel11), array(catvec['dec'], float))
    xusno, yusno = [], []
    for i in apixcut:
        xusno.append(float(string.split(i)[0]))
        yusno.append(float(string.split(i)[1]))
    xusno, yusno = array(xusno), array(yusno)
    #################################################################
    if verbose:
        iraf.tvmark(1, 'STDIN', Stdin=list(apixcut), mark="circle", number='yes', label='no', radii=8, nxoffse=5,
                    nyoffse=5, color=204, txsize=2)
        raw_input('brightest ' + str(number) + ' objects')

    ##############    sextractor   ##################
    if len(xpix) >= number:
        cm = array(cm, float)
        xpix = xpix[argsort(cm)][0:number]
        ypix = ypix[argsort(cm)][0:number]
        fw = fw[argsort(cm)][0:number]
        ell = ell[argsort(cm)][0:number]
        cm = cm[argsort(cm)][0:number]
    if verbose:
        sexpix = []
        for i in range(0, len(xpix)):
            sexpix.append(str(xpix[i]) + ' ' + str(ypix[i]))
        iraf.tvmark(1, 'STDIN', Stdin=list(sexpix), mark="circle", number='yes', label='no', radii=8, nxoffse=5,
                    nyoffse=5, color=206, txsize=2)
        raw_input('print sex ' + str(len(sexpix)))

    xsex, ysex = array(xpix), array(ypix)
    fwsex = array(fw)
    ellsex = array(ell)
    #####################################################################
    max_sep = tollerance1
    xdist, ydist = [], []
    for i in range(len(xusno)):
        dist = sqrt((xusno[i] - xsex) ** 2 + (yusno[i] - ysex) ** 2)
        idist = argmin(dist)
        if dist[idist] < max_sep:
            xdist.append(xusno[i] - xsex[idist])
            ydist.append(yusno[i] - ysex[idist])
    if len(xdist) >= 2:
        xoff, xstd = round(median(xdist), 2), round(std(xdist), 2)
        yoff, ystd = round(median(ydist), 2), round(std(ydist), 2)
        _xdist, _ydist = array(xdist), array(ydist)
        __xdist = compress((abs(_xdist - xoff) < 3 * xstd) &
                           (abs(_ydist - yoff) < 3 * ystd), _xdist)
        __ydist = compress((abs(_xdist - xoff) < 3 * xstd) &
                           (abs(_ydist - yoff) < 3 * ystd), _ydist)
        if len(__xdist) >= 2:
            xoff, xstd = round(median(__xdist), 2), round(std(__xdist), 2)
            yoff, ystd = round(median(__ydist), 2), round(std(__ydist), 2)
        else:
            xoff, yoff = 0, 0
    else:
        xoff, yoff = 0, 0
    if isnan(xoff):
        xoff = 0
    if isnan(yoff):
        yoff = 0
    _CRPIX1 = readkey3(hdr, 'CRPIX1')
    _CRPIX2 = readkey3(hdr, 'CRPIX2')
    ntt.util.updateheader(img, 0, {'CRPIX1': [_CRPIX1 - xoff, '']})
    ntt.util.updateheader(img, 0, {'CRPIX2': [_CRPIX2 - yoff, '']})
    xusno2_new = xusno - xoff
    yusno2_new = yusno - yoff
    #####################################################################
    max_sep = tollerance2
    fwhm = []
    fwhm2 = []
    ell = []
    xref = []
    iraf.tv(_doprint=0, Stdout=0)
    iraf.tv.rimexam.backgrou = 'yes'
    vettoretran = []
    # print "LOGX:: xusno2_new: %(xusno2_new)s" % locals()
    for i in range(len(xusno2_new)):
        dist = sqrt((xusno2_new[i] - xsex) ** 2 + (yusno2_new[i] - ysex) ** 2)
        idist = argmin(dist)
        if dist[idist] < max_sep:
            xref.append(xsex[idist])
            vettoretran.append(
                str(rausno[i]) + ' ' + str(decusno[i]) + ' ' + str(xsex[idist]) + ' ' + str(ysex[idist]) + ' \n')
            fwhm.append(fwsex[idist])
            ell.append(ellsex[idist])
            if imex:
                gg = open('tmp.one', 'w')
                gg.write(str(xsex[idist]) + ' ' + str(ysex[idist]) + '\n')
                gg.close()
                ime = iraf.imexam(input=img, frame=1, logfile='', keeplog='yes', imagecur='tmp.one', wcs='logical',
                                  use_disp='no', Stdout=1)
                try:
                    _fwhm2 = median(compress(array(string.split(ime[3])[-3:], float) < 99,
                                             (array(string.split(ime[3])[-3:], float))))
                    fwhm2.append(_fwhm2)
                except:
                    pass
    # print "LOGX:: xref: %(xref)s" % locals()
    # print "LOGX:: nummin: %(nummin)s" % locals()
    if len(xref) >= nummin:
        _ccmap1 = iraf.ccmap('STDIN', 'STDOUT', images=img, Stdin=vettoretran, fitgeome=fitgeo, xcolum=3, xxorder=2,
                             yyorder=2, ycolum=4, lngcolum=1, latcolumn=2, lngunit='degrees', update='No',
                             interact='No', maxiter=3, Stdout=1)
        if 'rms' in _ccmap1[_ccmap1.index('Wcs mapping status') + 1]:
            try:
                rmsx, rmsy = array(
                    string.split(string.split(_ccmap1[_ccmap1.index('Wcs mapping status') + 1], ':')[-1])[0:2], float)
                if float(rmsx) < 2. and float(rmsy) < 2.:
                    _ccmap1 = iraf.ccmap('STDIN', 'STDOUT', images=img, Stdin=vettoretran, fitgeome=fitgeo, xcolum=3,
                                         xxorder=2, yyorder=2, ycolum=4, lngcolum=1, latcolumn=2, lngunit='degrees',
                                         update='Yes', interact='No', maxiter=3, Stdout=1)
                    xy = iraf.wcsctran('STDIN', output="STDOUT", Stdin=vettoretran, Stdout=1, image=img,
                                       inwcs='physical', outwcs='world', column="3 4", formats='%10.6f %10.6f',
                                       verbose='yes')[3:]
                    print('\n### update astrometry with non linear order')
            except:
                rmsx, rmsy = array(
                    string.split(string.split(_ccmap1[_ccmap1.index('Wcs mapping status') + 1], ':')[-1])[0:2])
            rasys = 0.0000278
            decsys = 0.0000278
        else:
            rmsx, rmsy, rasys, decsys = 999, 999, 999, 999
    else:
        rmsx, rmsy = 990, 999
        rasys, decsys = 999, 999
    return rmsx, rmsy, len(xref), fwhm, ell, fwhm2, bkg, rasys, decsys


##########################################################################

def querydss(RA, DEC, ImSize='10', output='tmpdss.fits'):
    # print "LOGX:: Entering `querydss` method/function in %(__file__)s" %
    # globals()
    import urllib
    import urllib2

    dss = 'DSS1'
    query_url = "http://archive.eso.org/dss/dss/image?ra=" + str(RA) + "&dec=" + str(DEC) + \
                "&x=" + str(ImSize) + "&y=" + str(ImSize) + \
        "&Sky-Survey=" + dss + "&mime-type=download-fits"
    U = urllib2.urlopen(query_url)
    R0 = U.read()
    U.close()
    f = open(output, 'w')
    f.write(R0)
    f.close()
    return output


###################################################################

def crossmatchxy(_xx0, _yy0, _xx1, _yy1, tollerance):  # pixel,pixel,pixel,pixel,pixel
    # print "LOGX:: Entering `crossmatchxy` method/function in %(__file__)s" %
    # globals()
    import numpy as np

    distvec = []
    pos0 = []
    pos1 = []
    for jj in range(0, len(_xx0)):
        distance = np.sqrt((_xx0[jj] - _xx1) ** 2 + (_yy0[jj] - _yy1) ** 2)
        if min(distance) <= tollerance:
            distvec.append(np.min(distance))
            pos0.append(jj)
            pos1.append(np.argmin(distance))
    return distvec, pos0, pos1

#######################################################

################################################################################

def deg2HMS(ra='', dec='', round=False):
      import string
      RA, DEC= '', ''
      if dec:
          if string.count(str(dec),':')==2:
              dec00=string.split(dec,':')
              dec0, dec1, dec2 = float(dec00[0]),float(dec00[1]),float(dec00[2])
              if '-' in str(dec0):
                 DEC=(-1)*((dec2/60.+dec1)/60.+((-1)*dec0))
              else:
                 DEC=(dec2/60.+dec1)/60.+dec0
          else:
             dec0 = abs(int(dec))
             dec1=int((abs(dec)-abs(dec0))*(60))
             dec2=((((abs(dec))-abs(dec0))*60)-abs(dec1))*60
             if str(dec)[0]=='-':
                DEC = '-'+'00'[len(str(dec0)):]+str(dec0)+':'+'00'[len(str(dec1)):]+str(dec1)+':'+'00'[len(str(int(dec2))):]+str(dec2)[:6]
             else:
                DEC = '+'+'00'[len(str(dec0)):]+str(dec0)+':'+'00'[len(str(dec1)):]+str(dec1)+':'+'00'[len(str(int(dec2))):]+str(dec2)[:6]
      if ra:
          if string.count(str(ra),':')==2:
              ra00=string.split(ra,':')
              ra0,ra1,ra2=float(ra00[0]),float(ra00[1]),float(ra00[2])
              RA=((ra2/60.+ra1)/60.+ra0)*15.
          else:
              ra0=int(ra/15.)
              ra1=int(((ra/15.)-ra0)*(60))
              ra2=((((ra/15.)-ra0)*60)-ra1)*60
              RA='00'[len(str(ra0)):]+str(ra0)+':'+'00'[len(str(ra1)):]+str(ra1)+':'+'00'[len(str(int(ra2))):]+str(ra2)[:6]
      if ra and dec:          return RA, DEC
      else:                   return RA or DEC

##############################################################################
