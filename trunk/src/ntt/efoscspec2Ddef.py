def aperture(img):

    from astropy.io import fits as pyfits
    import re
    import os

    hdr = pyfits.open(img)[0].header
    xmax = hdr['NAXIS1']
    center = float(xmax) / 2.
    xmin = -500
    img2 = re.sub('.fits', '', img)
    line = "# Sun 13:10:40 16-Jun-2013\nbegin	aperture " + img2 + " 1 " + str(center) + "  500.0\n" + \
           "	 image	" + img2 + "\n	aperture	1\n	beam	1\n	center	" + str(center) + "  500.0\n" + \
           "	 low	-450. " + str(xmin) + "\n	high	450. " + str(xmax) + "\n" \
        "	 background\n	 xmin -850.\n" + \
           "	 xmax 850.\n	 function chebyshev\n		order 1\n		sample *\n" + \
           "	 naverage -3\n	 niterate 0\n		low_reject 3.\n		high_reject 3.\n" + \
           "	 grow 0.\n	 axis	1\n	 curve	5\n		2.\n		1.\n" + \
           "	 1.\n		1020.\n		0\n"
    if not os.path.isdir('database/'):
        os.mkdir('database/')

    f = open('database/ap' + img2, 'w')
    f.write(line)
    f.close()

##########################################################################


def choseflat(obj, listflat, setup, _JD0, _interactive):
    # print "LOGX:: Entering `choseflat` method/function in %(__file__)s" %
    # globals()
    import ntt
    import numpy as np

    differences = []
    obidflat = []
    flatgood = []
    _OBID = ntt.util.readkey3(ntt.util.readhdr(obj), 'esoid')
    listflat = ntt.sofiphotredudef.sortbyJD(listflat)
    for flat in listflat:
        hdrf = ntt.util.readhdr(flat)
        _JDf = ntt.util.readkey3(hdrf, 'JD')
        differences.append(np.abs(float(_JD0) - float(_JDf)))
        obidflat.append(ntt.util.readkey3(hdrf, 'esoid'))
    if _interactive:
        for flat in listflat:
            _JDf = ntt.util.readkey3(ntt.util.readhdr(flat), 'JD')
            ntt.util.display_image(flat, 1, '', '', False)
            answ = raw_input(
                '### good/bad/stop(enough files, go on) [[g],b,s]')
            if not answ:
                answ = 'g'
            if answ in ['G', 'g', 'good', 'Good']:
                flatgood.append(flat)
            elif answ in ['stop', 'S', 'Stop', 's']:
                break
    else:
        if obidflat.count(_OBID) >= 3:
            flatgood = list(np.compress(np.array(obidflat) == _OBID, listflat))
            print('### Flat field in the same OB !!')
        else:
            print('### ', str(_OBID))
            inds = np.array(differences).argsort()
            for i in range(0, 3):
                flatgood.append(listflat[inds[i]])
    return flatgood


def continumsub(imagefile, _order1, _order2):
    # print "LOGX:: Entering `continumsub` method/function in %(__file__)s" %
    # globals()
    import ntt
    from pyraf import iraf
    iraf.noao(_doprint=0, Stdout=0)
    iraf.imred(_doprint=0, Stdout=0)
    iraf.specred(_doprint=0, Stdout=0)
    toforget = ['specred.continuum']
    for t in toforget:
        iraf.unlearn(t)
    ntt.util.delete('tsky.fits')
    iraf.specred.continuum(imagefile, output='tsky.fits', type='difference',
                           interact='no', function='legendre', niterat=300, low_rej=3, high_re=2, sample='*',
                           order=_order1, ask='YES')
    ntt.util.delete(imagefile)
    iraf.continuum('tsky.fits', output=imagefile, type='difference',
                   interact='no', function='spline1', overrid='yes', niterat=10, low_rej=3, high_re=1, sample='*',
                   order=_order2, ask='YES')
    ntt.util.delete('tsky.fits')
    return imagefile


def skyfrom2d(fitsfile, skyfile, interac=True):
    # print "LOGX:: Entering `skyfrom2d` method/function in %(__file__)s" %
    # globals()
    import ntt
    from astropy.io import fits as pyfits
    import numpy as np

    yy1 = pyfits.open(fitsfile)[0].data[:, :].mean(1)
    crval2 = pyfits.open(fitsfile)[0].header.get('CRVAL2')
    cd2 = pyfits.open(fitsfile)[0].header.get('CD2_2')

    ntt.util.delete('new3.fits')
    hdu = pyfits.PrimaryHDU(yy1)
    hdulist = pyfits.HDUList([hdu])
#    hdulist[0].header.update('CRVAL1', crval2)
#    hdulist[0].header.update('CD1_1', cd2)
    hdulist[0].header['CRVAL1']= crval2
    hdulist[0].header['CD1_1'] = cd2
    hdulist.writeto('new3.fits')
    hdulist.close()

    fitsfile = ntt.efoscspec2Ddef.continumsub('new3.fits', 6, 1)
    yy1 = pyfits.open(fitsfile)[0].data
    xx1 = np.arange(len(yy1))
    aa1 = crval2 + (xx1) * cd2
    ntt.util.delete('new3.fits')

    skyff = pyfits.open(skyfile)[0].data
    crval1 = pyfits.open(skyfile)[0].header.get('CRVAL1')
    cd1 = pyfits.open(skyfile)[0].header.get('CD1_1')
    skyxx = np.arange(len(skyff))
    skyaa = crval1 + (skyxx) * cd1
    shift = ntt.efoscspec2Ddef.checkwavelength_arc(
        aa1, yy1, skyaa, skyff, 5500, 6500, interac)
    return shift


def checkwavelength_arc(xx1, yy1, xx2, yy2, xmin, xmax, inter=True):

    import numpy as np

    minimo = max(min(xx1), min(xx2)) + 50
    massimo = min(max(xx1), max(xx2)) - 50
    yy1 = [0 if e < 0 else e for e in np.array(yy1)]
    yy2 = [0 if e < 0 else e for e in np.array(yy2)]
    _shift, integral = [], []
    for shift in range(-500, 500, 1):
        xxnew = xx1 + shift / 10.
        yy2interp = np.interp(xxnew, xx2, yy2)
        yy2timesyy = yy2interp * yy1
        xxcut = np.compress((np.array(xxnew) >= minimo) & (
            np.array(xxnew) <= massimo), np.array(xxnew))
        yycut = np.compress((np.array(xxnew) >= minimo) & (
            np.array(xxnew) <= massimo), np.array(yy2timesyy))
        integrale = np.trapz(yycut, xxcut)
        integral.append(integrale)
        _shift.append(shift / 10.)
    result = _shift[integral.index(max(integral))]
    if inter:
        # import matplotlib as mpl
        #   mpl.use("TKAgg")
        import pylab as pl
        pl.ion()
        pl.clf()
        ratio = np.trapz(yy1, xx1) / np.trapz(yy2, xx2)
        yy3 = np.array(yy2) * float(ratio)
        xx4 = xx1 + result
        pl.plot(xx1, yy1, label='spectrum')
        pl.plot(xx2, yy3, label='reference sky', lw=2.5)
        pl.plot(xx4, yy1, label='shifted spectrum')
        pl.legend(numpoints=1, markerscale=1.5)
        if xmin != '' and xmax != '':
            pl.xlim(xmin, xmax)
    return result


def imreplace_region(img):

    import ntt
    from pyraf import iraf

    _grism = ntt.util.readkey3(ntt.util.readhdr(img), 'grism')
    iraf.imutil(_doprint=0, Stdout=0)
    iraf.unlearn('imutil.imreplace')
    if _grism == 'Gr13':
        iraf.imutil.imreplace(
            img + '[*,1:200]', value=1, lower='INDEF', upper='INDEF')
        print('### replace pixel 1:200 with 1 (y axes)')
    elif _grism in ['Gr11']:
        iraf.imutil.imreplace(
            img + '[*,1:300]', value=1, lower='INDEF', upper='INDEF')
        print('### replace pixel 1:300 with 1 (y axes)')
    elif _grism in ['Gr18', 'Gr20']:
        iraf.imutil.imreplace(
            img + '[*,1:50]', value=1, lower='INDEF', upper='INDEF')
        print('### replace pixel 1:50 with 1 (y axes)')
    else:
        print('### no replace ')


def efoscspecreduction(files, _interactive, _dobias, _doflat, _listflat, _listbias, _listarc, _cosmic, _verbose=False):

    import ntt
    import string
    import re
    import os
    import glob
    import sys
    from astropy.io import fits as pyfits

    import numpy as np
    from pyraf import iraf
    iraf.noao(_doprint=0, Stdout=0)
    iraf.imred(_doprint=0, Stdout=0)
    iraf.ccdred(_doprint=0, Stdout=0)
    iraf.twodspec(_doprint=0, Stdout=0)
    iraf.longslit(_doprint=0, Stdout=0)
    iraf.specred(_doprint=0, Stdout=0)
    toforget = ['ccdred.flatcombine', 'ccdred.zerocombine', 'ccdproc', 'specred.apall', 'longslit.identify',
                'longslit.reidentify',
                'specred.standard', 'longslit.fitcoords', 'specred.transform', 'specred.response']
    for t in toforget:
        iraf.unlearn(t)
    iraf.longslit.dispaxi = 2
    iraf.longslit.mode = 'h'
    iraf.identify.fwidth = 7
    iraf.specred.dispaxi = 2
    iraf.specred.mode = 'h'
    iraf.ccdproc.darkcor = 'no'
    iraf.ccdproc.fixpix = 'no'
    iraf.ccdproc.trim = 'no'
    iraf.ccdproc.flatcor = 'no'
    iraf.ccdproc.overscan = 'no'
    iraf.ccdproc.zerocor = 'no'
    iraf.ccdproc.biassec = ''
    iraf.ccdproc.ccdtype = ''
    iraf.ccdred.instrument = "/dev/null"
    iraf.ccdred.instrument = "ccddb$kpno/camera.dat"

    iraf.set(direc=ntt.__path__[0] + '/')
    if _interactive:
        _inter = 'yes'
    else:
        _inter = 'no'
    if _verbose:
        iraf.ccdred.verbose = 'yes'
        iraf.specred.verbose = 'yes'
    else:
        iraf.specred.verbose = 'no'
        iraf.ccdred.verbose = 'no'
    import datetime

    now = datetime.datetime.now()
    datenow = now.strftime('20%y%m%d%H%M')
    MJDtoday = 55927 + (datetime.date.today() - datetime.date(2012, 01, 01)).days
    outputlist = []
    _gain = ntt.util.readkey3(ntt.util.readhdr(
        re.sub('\n', '', files[0])), 'gain')
    _rdnoise = ntt.util.readkey3(
        ntt.util.readhdr(re.sub('\n', '', files[0])), 'ron')
    _biassec0 = '[3:1010,1026:1029]'
    objectlist = {}
    biaslist = {}
    flatlist = {}
    flatlistd = {}
    arclist = {}
    for img in files:
        img = re.sub('\n', '', img)
        hdr = ntt.util.readhdr(img)
        _naxis1 = hdr.get('NAXIS1')
        _naxis2 = hdr.get('NAXIS2')
        if _naxis1 != 1030 or _naxis2 != 1030:
            print('### warning dimension of ' + str(img) + ' not good !!!!')
            ntt.util.writeinthelog(
                'image ' + str(img) + ' different dimension =\n', './logNTT.txt')
        else:
            _date = ntt.util.readkey3(hdr, 'date-night')
            _exptime = ntt.util.readkey3(hdr, 'exptime')
            _type = ''
            if float(_exptime) == 0.0:
                if _date not in biaslist:
                    biaslist[_date] = [img]
                else:
                    biaslist[_date].append(img)
                _type = 'bias'
            if not _type:
                _imagetype = ntt.util.readkey3(hdr, 'tech')
                if _imagetype != 'SPECTRUM':
                    _type = 'photometric data'
            if not _type:
                _object = ntt.util.readkey3(hdr, 'object')
                _grism = ntt.util.readkey3(hdr, 'grism')
                _filter = ntt.util.readkey3(hdr, 'filter')
                _slit = ntt.util.readkey3(hdr, 'slit')
                if _object.lower() == 'flat':
                    if (_grism, _filter, _slit) not in flatlist:
                        flatlist[_grism, _filter, _slit] = [img]
                    else:
                        flatlist[_grism, _filter, _slit].append(img)
                    _type = 'flat'
                elif _object.lower() == 'dome':
                    if (_grism, _filter, _slit) not in flatlistd:
                        flatlistd[_grism, _filter, _slit] = [img]
                    else:
                        flatlistd[_grism, _filter, _slit].append(img)
                    _type = 'flat'
                elif _object.lower() == 'wave':
                    if (_grism, _filter, _slit) not in arclist:
                        arclist[_grism, _filter, _slit] = [img]
                    else:
                        arclist[_grism, _filter, _slit].append(img)
                    _type = 'arc'
            if not _type:
                if (_grism, _filter, _slit) not in objectlist:
                    objectlist[_grism, _filter, _slit] = [img]
                else:
                    objectlist[_grism, _filter, _slit].append(img)
                _type = 'objects'
            if not _type:
                print('### ', _object, _filter, _grism, _slit, _imagetype)

    if _verbose:
        if _listarc:
            print('### arclist= \n', _listarc)
        else:
            print('### arclist= \n', arclist)
        print('### flatlist= \n', flatlist)
        print('### objectlist= \n', objectlist)
    # #############################################
    check1, check2 = [], []
    for ii in objectlist.keys():
        if ii not in arclist:
            check1.append(ii)
        if ii not in flatlist:
            check2.append(ii)

    if len(check2) > 0:
        for ii in check2:
            print('\n###Warning: flat with setup ' + str(ii) + ' are missing')
            answ = raw_input(
                '\n### skip this setup from reduction [s] or exit [e] or go on [g] ? [s] ')
            if not answ:
                answ = 's'
            if answ in ['s', 'S']:
                try:
                    objectlist.pop(ii)
                except:
                    pass
            elif answ in ['e']:
                sys.exit(
                    '\n### add to the following directory the missing files and try again.')
            else:
                pass

    if len(check1) > 0:
        for ii in check1:
            print('\n###Warning: arc with setup ' + str(ii) + ' are missing')
            answ = raw_input(
                '\n### skip this setup from reduction [s] or exit [e] or go on [g] ?  [s] ')
            if not answ:
                answ = 's'
            if answ in ['s', 'S']:
                try:
                    objectlist.pop(ii)
                except:
                    pass
            elif answ in ['e']:
                sys.exit(
                    '\n### download the missing calibrations from ESO archive and try again.')
            else:
                pass

    ###### masterbias  #################
    if _dobias:
        if _listbias:
            masterbiaslist = _listbias
        else:
            masterbiaslist = []
            if biaslist:
                for _date in biaslist:
                    biaslist[_date] = ntt.efoscphotredudef.rejectbias(
                        biaslist[_date], _interactive)
                    if len(biaslist[_date]) >= 3:
                        masterbiasfile = 'bias_' + str(_date) + '.fits'
                        f = open('biaslist', 'w')
                        for img in biaslist[_date]:
                            f.write(img + '\n')
                        f.close()
                        try:
                            ntt.util.delete(masterbiasfile)
                            iraf.ccdred.zerocombine('@biaslist', output=masterbiasfile, combine='median',
                                                    reject='ccdclip', ccdtype=' ', rdnoise=_rdnoise, gain=_gain,
                                                    process='no', Stdout=1)
                            ntt.util.correctcard(masterbiasfile)
                            masterbiaslist.append(masterbiasfile)
                            num = 0
                            for img in biaslist[_date]:
                                num = num + 1
                                ntt.util.updateheader(masterbiasfile, 0, {
                                    'PROV' + str(num): [ntt.util.readkey3(ntt.util.readhdr(img), 'ARCFILE'),
                                                        'Originating file']})
                                ntt.util.updateheader(masterbiasfile, 0, {
                                    'TRACE' + str(num): [ntt.util.readkey3(ntt.util.readhdr(img), 'ARCFILE'),
                                                         'Originating file']})
                            ntt.util.updateheader(masterbiasfile, 0,
                                                  {'SINGLEXP': [False, 'TRUE if resulting from single exposure']})
                            ntt.util.updateheader(masterbiasfile, 0,
                                                  {'M_EPOCH': [True, 'TRUE if resulting from multiple epochs']})
                        except:
                            ntt.util.writeinthelog(
                                'Warning ' +
                                str(biaslist[_date]) +
                                ' problem with this list of bias \n',
                                './logNTT.txt')
    ################################################################
    for setup in objectlist:
        listobject = objectlist[setup]
        _grism0 = setup[0]
        if _grism0 == 'Gr16':
            _trimsec0 = '[100:950,1:950]'
        elif _grism0 == 'Gr13':
            if setup[1] == 'Free':
                _trimsec0 = '[100:950,1:1015]'
            elif setup[1] == 'GG495':
                _trimsec0 = '[100:950,250:1015]'
            elif setup[1] == 'OG530':
                _trimsec0 = '[100:950,300:1015]'
        elif _grism0 in ['Gr11', 'Gr18', 'Gr20']:
            if setup[1] == 'GG495':
                _trimsec0 = '[100:710,40:1015]'
            else:
                _trimsec0 = '[100:950,5:1015]'
        else:
            _trimsec0 = '[100:950,5:1015]'

        if setup[0] == 'Gr16' and setup[1] == 'Free':
            _order = 20
            _sample = '*'
        elif setup[0] == 'Gr16' and setup[1] == 'GG495':
            _order = 60
            _sample = '*'
        elif setup[0] == 'Gr16' and setup[1] == 'OG530':
            _order = 70
            _sample = '*'
        elif setup[0] in ['Gr11', 'Gr18', 'Gr20'] and setup[1] == 'Free':
            _order = 35
            _sample = '*'
        elif setup[0] == 'Gr20' and setup[1] == 'GG495':
            _order = 35
            _sample = '*'
        elif setup[0] == 'Gr13' and setup[1] == 'Free':
            _order = 90
            _sample = '*'
        else:
            _order = 35
            _sample = '*'
        if setup in flatlist:
            listflat = flatlist[setup]
        elif setup in flatlistd:
            listflat = flatlistd[setup]
        else:
            listflat = ''
        if setup in arclist:
            listarc = arclist[setup]
        else:
            listarc = ''
        if _interactive:
            print('### ' + str(listobject))
            print('### ' + str(setup))
            answ = raw_input('### do you want to reduce this setup [[y],n] ? ')
            if not answ:
                answ = 'y'
        else:
            answ = 'y'
        if answ in ['YES', 'yes', 'y', 'Y', 'Yes']:
            ################# bias #########################################
            if _dobias:
                _zerocor = 'yes'
                masterbiaslist = masterbiaslist + \
                    glob.glob(ntt.__path__[0] + '/archive/efosc/bias/*')
                allframes = listobject + listflat + listarc
                biasneeded = []
                tmasterbiaslist = []
                for img0 in allframes:
                    bias0 = ntt.util.choseclosest(img0, masterbiaslist)
                    if bias0 not in biasneeded:
                        biasneeded.append(bias0)
                for masterbias in biasneeded:
                    if masterbias[0] == '/':
                        os.system('cp ' + masterbias + ' ' +
                                  string.split(masterbias, '/')[-1])
                        masterbias = string.split(masterbias, '/')[-1]

                    tmasterbias = re.sub('.fits', '_' + setup[0] + '_' + setup[1] + '_' + str(MJDtoday)
                                         + '.fits', masterbias)

                    tmasterbiaslist.append(tmasterbias)

                    if not ntt.util.readkey3(ntt.util.readhdr(masterbias), 'TRIM'):
                        ntt.util.delete(tmasterbias)
                        iraf.ccdproc(masterbias, output=tmasterbias, overscan="no", trim="yes", zerocor='no',
                                     flatcor='no', zero='', ccdtype='', fixpix='no', trimsec=_trimsec0, biassec='',
                                     readaxi='column', Stdout=1)
                    else:
                        os.system('cp ' + masterbias + ' ' + tmasterbias)

                    ntt.util.correctcard(tmasterbias)
                    headervecb = {'M_EPOCH': [True, 'TRUE if resulting from multiple epochs'],
                                  'SINGLEXP': [False, 'TRUE if resulting from single exposure'],
                                  'FILETYPE': [21201, 'bias']}
                    ntt.util.updateheader(tmasterbias, 0, headervecb)
                    try:
                        pyv = int(re.sub('\.', '', str(pyfits.__version__))[:2])
                    except:
                        pyv = 40 # astropy.pyfits do not have a version, set high number to avoud next if

                    if pyv <= 30:
                        ntt.util.updateheader(tmasterbias, 0, {'HIERARCH ESO INS FILT1 NAME': [setup[1], 'Filter name.']})
                        ntt.util.updateheader(tmasterbias, 0, {
                            'HIERARCH ESO INS GRIS1 NAME': ['Gr#' + re.sub('Gr', '', str(setup[0])), 'OPTIi name.']})
                    else:
                        # was this a bag ? (img instead of tmasterbias)
                        imm = pyfits.open(tmasterbias, mode='update')
                        header = imm[0].header
                        try:
                            header.pop('ESO INS FILT1 NAME')
                        except:
                            pass
                        try:
                            header.pop('ESO INS GRIS1 NAME')
                        except:
                            pass
                        header['ESO INS GRIS1 NAME'] = (
                            str(setup[1]), 'Filter name.')
                        header['ESO INS FILT1 NAME'] = (
                            str(setup[0]), 'OPTIi name.')
                        imm.flush()
                        imm.close()
                    if tmasterbias not in outputlist:
                        outputlist.append(tmasterbias)
                    if _verbose:
                        print(tmasterbias)
            else:
                _zerocor = 'no'

            print(tmasterbiaslist)

            ################   make all flats of this setup   #################
            if _doflat:
                _flatcor = 'yes'
                if _listflat:
                    flatgood = _listflat  # flat list from reducer
                elif listflat:  # flat in the  raw data
                    flatgood = []
                    IDflat = {}
                    for ff in listflat:
                        OBID = ntt.util.readkey3(ntt.util.readhdr(ff), 'esoid')
                        if OBID not in IDflat:
                            IDflat[OBID] = []
                        IDflat[OBID].append(ff)

                    for ID in IDflat.keys():
                        if len(IDflat[ID]) >= 3:
                            _date = ntt.util.readkey3(
                                ntt.util.readhdr(IDflat[ID][0]), 'date-night')
                            masterflat = 'flat_' + str(_date) + '_' + str(setup[0]) + '_' + str(setup[1]) + '_' + str(
                                setup[2]) + '_' + str(ID) + '_' + str(MJDtoday) + '.fits'
                            f = open('_flatlist', 'w')
                            o = open('_oflatlist', 'w')
                            for ff in IDflat[ID]:
                                f.write(ff + '\n')
                                o.write('o' + ff + '\n')
                                ntt.util.delete('o' + ff)
                            f.close()
                            o.close()
                            ntt.util.delete(masterflat)
                            if _dobias:
                                tmasterbias = ntt.util.choseclosest(
                                    IDflat[ID][0], tmasterbiaslist)
                            else:
                                tmasterbias = ''

                            iraf.ccdproc('@_flatlist', output='@_oflatlist', overscan="no", trim="yes", darkcor='no',
                                         fixpix='no', zerocor=_zerocor, flatcor="no", trimsec=_trimsec0,
                                         biassec='', zero=tmasterbias, readaxi='column', ccdtype='', Stdout=1)
                            iraf.ccdred.flatcombine('"@_oflatlist"', output=masterflat, combine='average',
                                                    reject='none', ccdtype=' ', rdnoise=_rdnoise, gain=_gain,
                                                    process='no', Stdout=1)
                            ntt.util.correctcard(masterflat)
                            ntt.util.delete('_flatlist')
                            ntt.util.delete('_oflatlist')
                            if masterflat not in outputlist:
                                outputlist.append(masterflat)
                            hedvec = {'M_EPOCH': [True, 'TRUE if resulting from multiple epochs'],
                                      'ZEROCOR': [tmasterbias, ''],
                                      'SINGLEXP': [False, 'TRUE if resulting from single exposure'],
                                      'FILETYPE': [21102, 'flat field']}
                            ntt.util.updateheader(masterflat, 0, hedvec)
                            ntt.util.delete('n' + masterflat)

                            ntt.efoscspec2Ddef.aperture(masterflat)

                            iraf.specred.apflatten(masterflat, output='n' + masterflat, interac=_inter, find='no',
                                                   recenter='no', resize='no', edit='no', trace='no', fittrac='no',
                                                   fitspec='yes', flatten='yes', aperture='', pfit='fit2d',
                                                   clean='no', function='spline3', order=_order, sample='*', mode='ql')
                            # print _order, _inter, masterflat
                            # print minpixel, maxpixel
                            #raw_input('test flat')

                            ###################################################
                            # iraf.specred.response(masterflat, normaliz=masterflat + '[' + str(minpixel) + ':' + str(
                            #    maxpixel) + ',*]', response='n' + masterflat, interac=_inter, thresho='INDEF',
                            #                      sample=_sample, naverage=2, function='spline3', low_rej=3,
                            #                      high_rej=3, order=_order, niterat=20, grow=0, graphic='stdgraph')
                            ###################################################

                            ntt.efoscspec2Ddef.imreplace_region(
                                'n' + masterflat)
                            nmasterflat = 'n' + masterflat
                            if nmasterflat not in outputlist:
                                outputlist.append(nmasterflat)
                            if nmasterflat not in flatgood:
                                flatgood.append(nmasterflat)
                            ntt.util.updateheader(nmasterflat, 0, {'FILETYPE': [
                                                  21203, 'normalized flat field']})
                            ntt.util.updateheader(
                                nmasterflat, 0, {'TRACE1': [masterflat, 'Originating file']})
                            num = 0
                            for img in IDflat[ID]:
                                ntt.util.delete('o' + img)
                                num = num + 1
                                ntt.util.updateheader(masterflat, 0, {
                                    'PROV' + str(num): [ntt.util.readkey3(ntt.util.readhdr(img), 'ARCFILE'),
                                                        'Originating file']})
                                ntt.util.updateheader(nmasterflat, 0, {
                                    'PROV' + str(num): [ntt.util.readkey3(ntt.util.readhdr(img), 'ARCFILE'),
                                                        'Originating file']})
                                ntt.util.updateheader(masterflat, 0, {
                                    'TRACE' + str(num): [ntt.util.readkey3(ntt.util.readhdr(img), 'ARCFILE'),
                                                         'Originating file']})

                else:
                    flatgood = []
            else:
                _flatcor = 'no'
                flatgood = []
            ###################################################################
            for obj in listobject:
                hdr0 = ntt.util.readhdr(obj)
                _object0 = ntt.util.readkey3(hdr0, 'object')
                _object0 = re.sub(' ', '', _object0)
                _object0 = re.sub('/', '_', _object0)
                _date0 = ntt.util.readkey3(hdr0, 'date-night')
                nameout0 = str(_object0) + '_' + str(_date0)
                for _set in setup:
                    nameout0 = nameout0 + '_' + _set
                nameout0 = nameout0 + '_' + str(MJDtoday)
                nameout0 = ntt.util.name_duplicate(obj, nameout0, '')

                if _dobias:
                    tmasterbias = ntt.util.choseclosest(obj, tmasterbiaslist)
                else:
                    tmasterbias = ''

                print('### ' + str(obj) + '  -> ', str(nameout0), '\n')
                ntt.util.display_image(obj, 1, '', '', False)
                if len(flatgood) == 1:
                    _flatcor = 'yes'
                    nmasterflat = flatgood[0]
                elif len(flatgood) >= 2:
                    _JD0 = ntt.util.readkey3(hdr0, 'JD')
                    OBID = ntt.util.readkey3(hdr0, 'esoid')
                    _ra0 = ntt.util.readkey3(hdr0, 'RA')
                    _dec0 = ntt.util.readkey3(hdr0, 'DEC')
                    nmasterflat = ''
                    from numpy import cos, sin, arccos, array, pi, argmin

                    scal = pi / 180.
                    distance = []
                    JDvec = []
                    for ff in flatgood:
                        hdrf = ntt.util.readhdr(ff)
                        if ntt.util.readkey3(hdrf, 'esoid') == OBID:
                            nmasterflat = ff
                            _flatcor = 'yes'
                            break
                        else:
                            _JD1 = ntt.util.readkey3(hdrf, 'JD')
                            _ra1 = ntt.util.readkey3(hdrf, 'RA')
                            _dec1 = ntt.util.readkey3(hdrf, 'DEC')
                            # some flat,bias do not have ra and dec keywords
                            if not _ra1:
                                print('Warning: missing keyword RA')
                                _ra1 = 0
                            if not _dec1:
                                print('Warning: missing keyword DEC')
                                _dec1 = 0
                            distance.append(arccos(
                                sin(_dec1 * scal) * sin(_dec0 * scal) + cos(_dec1 * scal) * cos(_dec0 * scal) * cos(
                                    _ra1 - _ra0) * scal))
                            JDvec.append(np.abs(_JD0 - _JD1))
                    if _verbose:
                        print(JDvec)
                        print(distance)
                        print(flatgood)
                    if not nmasterflat:
                        #       select closer RA and DEC
                        #nmasterflat = flatgood[argmin(distance)]
                        # select closer in time
                        nmasterflat = flatgood[argmin(JDvec)]
                        _flatcor = 'yes'
                else:
                    _flatcor = 'no'
                    nmasterflat = ''
                ##################################################
                # pre-reduce  image   and arc
                ntt.util.delete(nameout0)
                iraf.ccdproc(obj, output=nameout0, overscan="no", trim="yes", zerocor=_zerocor, flatcor=_flatcor,
                             zero=tmasterbias,
                             trimsec=_trimsec0, biassec='', flat=nmasterflat, readaxi='column', Stdout=1)
                ntt.util.correctcard(nameout0)
                hedvec = {'M_EPOCH': [False, 'TRUE if resulting from multiple epochs'],
                          'SINGLEXP': [True, 'TRUE if resulting from single exposure'],
                          'ZEROCOR': [tmasterbias, ''], 'FLATCOR': [nmasterflat, ''],
                          'FILETYPE': [22104, 'pre-reduced spectroscopic frame'],
                          'NCOMBINE': [1, 'Number of raw science data'],
                          'PROV1': [ntt.util.readkey3(ntt.util.readhdr(nameout0), 'ARCFILE'), 'Originating file'],
                          'TRACE1': [ntt.util.readkey3(ntt.util.readhdr(nameout0), 'ARCFILE'), 'Originating file']}
                ntt.util.updateheader(nameout0, 0, hedvec)
                if nameout0 not in outputlist:
                    outputlist.append(nameout0)

                arcfile = ''
                if _listarc:
                    arcfile = ntt.util.searcharc(obj, _listarc)[0]
                if not arcfile:
                    arcfile = ntt.util.searcharc(obj, listarc)[0]
                if not arcfile:
                    arcfile = ntt.util.searcharc(obj, '')[0]
                if arcfile:
                    if arcfile[0] == '/':
                        os.system('cp ' + arcfile + ' arc_' + nameout0)
                    else:
                        ntt.util.delete('arc_' + nameout0)
                        iraf.ccdproc(arcfile, output='arc_' + nameout0, overscan="no", trim="yes", zerocor=_zerocor,
                                     flatcor=_flatcor, zero=tmasterbias,
                                     trimsec=_trimsec0, biassec='', flat=nmasterflat, readaxi='column', Stdout=1)
                        if not os.path.isfile('arc_' + nameout0):
                            os.system('cp ' + arcfile + ' arc_' + nameout0)
                        if 'arc_' + nameout0 not in outputlist:
                            outputlist.append('arc_' + nameout0)
                        ntt.util.correctcard('arc_' + nameout0)
                        hedvec = {'M_EPOCH': [False, 'TRUE if resulting from multiple epochs'],
                                  'SINGLEXP': [True, 'TRUE if resulting from single exposure'],
                                  'ZEROCOR': [tmasterbias, ''], 'FLATCOR': [nmasterflat, ''],
                                  'FILETYPE': [22104, 'pre-reduced spectroscopic frame'],
                                  'PROV1': [ntt.util.readkey3(ntt.util.readhdr('arc_' + nameout0), 'ARCFILE'),
                                            'Originating file'],
                                  'TRACE1': [ntt.util.readkey3(ntt.util.readhdr('arc_' + nameout0), 'ARCFILE'),
                                             'Originating file']}
                        ntt.util.updateheader('arc_' + nameout0, 0, hedvec)
                    arcfile = 'arc_' + nameout0
                else:
                    sys.exit('Warning: arcfile not found')
                ##############
                if _cosmic:
                    # print cosmic rays rejection
                    ntt.cosmics.lacos(nameout0, output='', gain=_gain, readn=_rdnoise, xorder=9, yorder=9, sigclip=4.5,
                                      sigfrac=0.5, objlim=1, verbose=True, interactive=False)
                    print('\n### cosmic rays rejections ........ done ')
                    ntt.util.updateheader(nameout0, 0, {
                        'LACOSMIC': [True, 'TRUE if Laplacian cosmic ray rejection has been applied to the image']})
                else:
                    ntt.util.updateheader(nameout0, 0, {
                        'LACOSMIC': [False, 'TRUE if Laplacian cosmic ray rejection has been applied to the image']})
                if arcfile:
                    arcref = ntt.util.searcharc(nameout0, '')[0]
                    if arcref:
                        os.system('cp ' + arcref + ' .')
                        arcref = string.split(arcref, '/')[-1]
                        if not os.path.isdir('database/'):
                            os.mkdir('database/')
                        if os.path.isfile(ntt.util.searcharc(nameout0, '')[1] +
                                          '/database/id' + re.sub('.fits', '', arcref)):
                            os.system('cp ' + ntt.util.searcharc(nameout0, '')[1] +
                                      '/database/id' + re.sub('.fits', '', arcref) + ' database/')
                        identific = iraf.longslit.reidentify(referenc=arcref, images=arcfile, interac=_inter,
                                                             section='column 10',
                                                             coordli='direc$standard/ident/Lines_HgCdHeNeAr600.dat',
                                                             overrid='yes', step=0, newaps='no', nsum=5, nlost=2,
                                                             cradius=10, mode='h', verbose='yes', Stdout=1)
                        if _interactive:
                            answ = raw_input(
                                '### do you like the identification [[y]/n]')
                            if not answ:
                                answ = 'y'
                        else:
                            answ = 'y'
                        if answ in ['n', 'N', 'no', 'NO', 'No']:
                            yy1 = pyfits.open(arcref)[0].data[
                                :, 400:410].mean(1)
                            xx1 = np.arange(len(yy1))
                            yy2 = pyfits.open(arcfile)[0].data[
                                :, 400:410].mean(1)
                            xx2 = np.arange(len(yy2))
                            _shift = ntt.efoscspec2Ddef.checkwavelength_arc(xx1, yy1, xx2, yy2, 6000, 7500,
                                                                            inter=_interactive) * (-1)
                            identific = iraf.longslit.reidentify(referenc=arcref, images=arcfile, interac='YES',
                                                                 section='column 10', shift=_shift, overrid='yes',
                                                                 step=0,
                                                                 coordli='direc$standard/ident/Lines_HgCdHeNeAr600.dat',
                                                                 newaps='no', nsum=5, nlost=2, cradius=10, mode='h',
                                                                 verbose='yes', Stdout=1)
                            answ = raw_input('### is it ok now [[y]/n]')
                            if not answ:
                                answ = 'y'
                            if answ in ['n', 'N', 'no', 'NO', 'No']:
                                sys.exit(
                                    'warning: line identification with some problems')
                    else:
                        identific = iraf.longslit.identify(images=arcfile, section='column 10',
                                                           coordli='direc$standard/ident/Lines_HgCdHeNeAr600.dat',
                                                           nsum=10, fwidth=7, order=5,
                                                           functio='legendre', cradius=10, mode='h', Stdout=1)

                    iraf.longslit.reidentify(referenc=arcfile, images=arcfile, interac='NO', section='column 10',
                                             newaps='yes', nsum=5, nlost=2,
                                             coordli='direc$standard/ident/Lines_HgCdHeNeAr600.dat', overrid='yes',
                                             step=10, cradius=10, mode='h', verbose='no', Stdout=1)
                    iraf.longslit.dispaxi = 2
                    qqq = iraf.longslit.fitcoords(images=re.sub('.fits', '', arcfile),
                                                  fitname=re.sub('.fits', '', arcfile), interac='no', combine='yes',
                                                  databas='database',
                                                  function='legendre', yorder=4, logfile='', plotfil='', mode='h')
                    ntt.util.delete('t' + nameout0)
                    iraf.specred.transform(input=nameout0, output='t' + nameout0, minput='',
                                           fitnames=re.sub('.fits', '', arcfile), databas='database',
                                           x1='INDEF', x2='INDEF', y1='INDEF', y2='INDEF', flux='yes',
                                           logfile='logfile')  # , mode='h')
                    #####################
                    ntt.util.delete('t' + arcfile)
                    iraf.specred.transform(input=arcfile, output='t' + arcfile, minput='',
                                           fitnames=re.sub('.fits', '', arcfile), databas='database',
                                           x1='INDEF', x2='INDEF', y1='INDEF', y2='INDEF', flux='yes',
                                           logfile='logfile')  # , mode='h')
                    specred = ntt.util.spectraresolution2(arcfile)
                    if specred:
                        ntt.util.updateheader('t' + nameout0, 0,
                                              {'SPEC_RES': [specred, 'Spectral resolving power']})
                    ntt.util.delete('t' + arcfile)
                    #####################
                    if 't' + nameout0 not in outputlist:
                        outputlist.append('t' + nameout0)
                    ntt.util.updateheader('t' + nameout0, 0,
                                          {'FILETYPE': [22106, 'wavelength calibrated 2D spectroscopic frame']})
                    ntt.util.updateheader(
                        't' + nameout0, 0, {'TRACE1': [nameout0, 'Originating file']})
                    ntt.util.updateheader(
                        't' + nameout0, 0, {'ARC': [arcfile, '']})
                    # specred=ntt.util.spectraresolution('t'+nameout0)
                    # print specred
                    #######################  check wavelength calibration #####
                    _skyfile = ntt.__path__[
                        0] + '/standard/ident/sky_' + setup[0] + '_' + setup[1] + '.fits'
                    if glob.glob(_skyfile) and \
                            float(ntt.util.readkey3(ntt.util.readhdr('t' + nameout0), 'exptime')) > 300.:
                        shift = ntt.efoscspec2Ddef.skyfrom2d(
                            't' + nameout0, _skyfile, _interactive)
                        zro = pyfits.open(
                            't' + nameout0)[0].header.get('CRVAL2')
                        print('\n### check wavelengh calibration, found a shift of ' + str(shift) + ' Angstrom \n')
                        if _interactive:
                            answ = raw_input(
                                '### do you want to correct the wavelengh calibration with this shift: ' + str(
                                    shift) + ' [[y]/n] ? ')
                            if not answ:
                                answ = 'y'
                        else:
                            answ = 'y'
                        if answ.lower() in ['y', 'yes']:
                            ntt.util.updateheader(
                                't' + nameout0, 0, {'CRVAL2': [zro + int(shift), '']})
                            ntt.util.updateheader(
                                't' + nameout0, 0, {'shift': [float(shift), '']})
                    else:
                        print('\n### exposure too short, the sky lines could be not visible \n')

                    if identific:
                        _rms = float(identific[-1].split()[-1])
                        _num = float(identific[-1].split()[2].split('/')[0])
                        hdr = ntt.util.readhdr('t' + nameout0)
                        hedvec = {'LAMRMS': [_rms * .1, 'residual RMS [nm]'],
                                  'LAMNLIN': [_num, 'Nb of arc lines used in the fit of the wavel. solution'],
                                  'SPEC_ERR': [(_rms * .1) / np.sqrt(float(_num)), 'statistical uncertainty'],
                                  'SPEC_SYE': [0.1, 'systematic error']}
                        try:
                            wavelmin = float(ntt.util.readkey3(hdr, 'CRVAL2')) + \
                                (0.5 - float(ntt.util.readkey3(hdr, 'CRPIX2'))) * float(
                                ntt.util.readkey3(hdr, 'CDELT2'))
                            wavelmax = float(ntt.util.readkey3(hdr, 'CRVAL2')) + (
                                (float(ntt.util.readkey3(hdr, 'NAXIS2')) + 0.5 -
                                 float(ntt.util.readkey3(hdr, 'CRPIX2'))) *
                                float(ntt.util.readkey3(hdr, 'CDELT2')))
                            hedvec['WAVELMIN'] = [
                                wavelmin * .1, '[nm] minimum wavelength']
                            hedvec['WAVELMAX'] = [
                                wavelmax * .1, '[nm]  maximum wavelength']
                            hedvec['XMIN'] = [
                                wavelmin, '[A] minimum wavelength']
                            hedvec['XMAX'] = [
                                wavelmax, '[A]  maximum wavelength']
                            hedvec['SPEC_BW'] = [
                                (wavelmax * .1) - (wavelmin * .1), '[nm] Bandpass Width Wmax - Wmin']
                            hedvec['SPEC_VAL'] = [
                                ((wavelmax * .1) + (wavelmin * .1)) / 2., '[nm] Mean Wavelength']
                            hedvec['SPEC_BIN'] = [((wavelmax * .1) - (wavelmin * .1)) /
                                                  (float(ntt.util.readkey3(
                                                      hdr, 'NAXIS2')) - 1),
                                                  'average spectral coordinate bin size [nm/pix]']
                            #                        print wavelmin,wavelmax
                            #                        print float(readkey3(hdr,'NAXIS2'))
                            #                        print ((wavelmax*.1)-(wavelmin*.1))/float(readkey3(hdr,'NAXIS2'))
                            #                        raw_input('sss')
                            hedvec['VOCLASS'] = [
                                'SPECTRUM V1.0', 'VO Data Model']
                            hedvec['VOPUB'] = ['ESO/SAF',
                                               'VO Publishing Authority']
                            hedvec['APERTURE'] = [2.778e-4 * float(re.sub('slit', '', ntt.util.readkey3(hdr, 'slit'))),
                                                  '[deg] Aperture diameter']
                        except:
                            pass
                        ntt.util.updateheader('t' + nameout0, 0, hedvec)

    print('\n### adding keywords for phase 3 ....... ')
    reduceddata = ntt.util.rangedata(outputlist)
    f = open('logfile_spec2D_' + str(reduceddata) +
             '_' + str(datenow) + '.raw.list', 'w')
    for img in outputlist:
        if img[-4:] == 'fits':
            ################################################
            ntt.util.updateheader(img, 0, {'DETRON ': [11.6, 'Readout noise per output (e-)']})

            try:
                pyv = int(re.sub('\.', '', str(pyfits.__version__))[:2])
            except:
                pyv = 40 # astropy.pyfits do not have a version, set high number to avoud next if

            if pyv <= 30:
                ntt.util.updateheader(img, 0,
                                      {'HIERARCH ESO DET OUT1 GAIN': [1.18, 'Conversion from electrons to ADU']})
                ntt.util.updateheader(img, 0, {'HIERARCH ESO DET OUT1 RON': [
                                      11.6, 'Readout noise per output (e-)']})
            else:
                imm = pyfits.open(img, mode='update')
                header = imm[0].header
                try:
                    header.pop('ESO DET OUT1 GAIN')
                except:
                    pass
                try:
                    header.pop('ESO DET OUT1 RON')
                except:
                    pass
                header['ESO DET OUT1 GAIN'] = (
                    1.18, 'Conversion from electrons to ADU')
                header['ESO DET OUT1 RON'] = (
                    11.6, 'Readout noise per output (e-)')
                imm.flush()
                imm.close()
            #################################################
            hdr = ntt.util.readhdr(img)
            ntt.util.phase3header(img)  # phase 3 definitions
            ntt.util.airmass(img)  # phase 3 definitions

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
            ntt.util.updateheader(img, 0, {'EFFRON': [ntt.util.readkey3(hdr, 'ron') *
                                                      np.sqrt(
                                                          1. + 1. / nflat + 1. / nbias),
                                                      'Effective readout noise per output (e-)']})

            hedvec = {'quality': ['Final', 'Final or fast reduction'],
                      'BUNIT': ['ADU', 'Physical unit of array values'],
                      'TEXPTIME': [ntt.util.readkey3(hdr, 'EXPTIME'), 'Total integ. time of all exposure']}
            if ntt.util.readkey3(hdr, 'MJD-OBS'):
                mjdend = float(ntt.util.readkey3(hdr, 'MJD-OBS')) + \
                    (float(ntt.util.readkey3(hdr, 'exptime')) * 0.00001 / (0.864))
                # mjdend=float(readkey3(hdr,'MJD-OBS'))+float(readkey3(hdr,'exptime'))+1.8)/(60.*60.*24.)
                # ??

                hedvec['MJD-END'] = [mjdend, 'End of observations (days)']
                # hedvec['TELAPSE']=[86400*(mjdend-float(readkey3(hdr,'MJD-OBS'))),'Total
                # elapsed time [days]']  # second or days ?
                hedvec['TELAPSE'] = [(mjdend - float(ntt.util.readkey3(hdr, 'MJD-OBS'))) * 60. * 60 * 24.,
                                     'Total elapsed time [days]']  # second or days ?
                hedvec['TMID'] = [
                    (mjdend + float(ntt.util.readkey3(hdr, 'MJD-OBS'))) / 2., '[d] MJD mid exposure']
                hedvec['TITLE'] = [ntt.util.readkey3(
                    hdr, 'object'), 'Dataset title']
                # hedvec['TITLE'] = [str(hedvec['TMID'][0])[0:9] + ' ' + str(ntt.util.readkey3(hdr, 'object')) +
                #                   ' ' + str(ntt.util.readkey3(hdr, 'grism')) + ' ' +
                #                   str(ntt.util.readkey3(hdr, 'filter')) + ' ' + str(ntt.util.readkey3(hdr, 'slit')),
                #                   'Dataset title']
            if ntt.util.readkey3(hdr, 'tech'):
                hedvec['PRODCATG'] = ['SCIENCE.IMAGE', 'Data product category']
            hedvec['EXT_OBJ'] = [False, 'TRUE if extened']
            hedvec['CONTNORM'] = [False, 'TRUE if normalised to the continuum']
            hedvec['TOT_FLUX'] = [
                False, 'TRUE if phot. cond. and all src flux is captured']
            hedvec['FLUXCAL'] = ['ABSOLUTE',
                                 'Certifies the validity of PHOTZP']
            hedvec['FLUXERR'] = [
                15.2, 'Fractional uncertainty of the flux [%]']
            hedvec['SPECSYS'] = ['TOPOCENT', 'Observed frame']
            hedvec['DISPELEM'] = [
                'Gr#' + re.sub('Gr', '', ntt.util.readkey3(hdr, 'grism')), 'Dispersive element name']
            ntt.util.updateheader(img, 0, hedvec)
            f.write(ntt.util.readkey3(hdr, 'arcfile') + '\n')
    f.close()

    return outputlist, 'logfile_spec2D_' + str(reduceddata) + '_' + str(datenow) + '.raw.list'

# ###############################################################################################
