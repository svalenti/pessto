def skysofifrom2d(fitsfile, skyfile):
    # print "LOGX:: Entering `skysofifrom2d` method/function in %(__file__)s"
    # % globals()
    import ntt
    from ntt.util import readhdr, readkey3, delete
    from numpy import mean, arange, compress
    try:        import pyfits
    except:     from astropy.io import fits as pyfits

    from numpy import interp as ninterp

    hdr = readhdr(fitsfile)
    _grism = readkey3(hdr, 'grism')
    if _grism == 'GR':
        _order1 = 10
    else:
        _order1 = 6

    yy1 = pyfits.open(fitsfile)[0].data[:, :].mean(1)
    crval2 = readkey3(hdr, 'CRVAL2')
    cd2 = readkey3(hdr, 'CD2_2')
    xx1 = arange(len(yy1))
    aa1 = crval2 + (xx1) * cd2
    yy1cut = compress((aa1 < 18400) | (aa1 > 18650), yy1)
    aa1cut = compress((aa1 < 18400) | (aa1 > 18650), aa1)
    yy1cut1 = compress((aa1cut < 11600) | (aa1cut > 11800), yy1cut)
    aa1cut1 = compress((aa1cut < 11600) | (aa1cut > 11800), aa1cut)
    yy1interp = ninterp(aa1, aa1cut1, yy1cut1)
    delete('_new3.fits')
    hdu = pyfits.PrimaryHDU(yy1interp)
    hdulist = pyfits.HDUList([hdu])
    hdulist[0].header.update('CRVAL1', crval2)
    hdulist[0].header.update('CD1_1', cd2)
    hdulist.writeto('_new3.fits')
    hdulist.close()
    fitsfile = ntt.efoscspec2Ddef.continumsub('_new3.fits', _order1, 1)
    yy1 = pyfits.open(fitsfile)[0].data
    crval2 = pyfits.open(fitsfile)[0].header.get('CRVAL1')
    cd2 = pyfits.open(fitsfile)[0].header.get('CD1_1')
    xx1 = arange(len(yy1))
    aa1 = crval2 + (xx1) * cd2

    skyff = pyfits.open(skyfile)[0].data
    crval1 = pyfits.open(skyfile)[0].header.get('CRVAL1')
    cd1 = pyfits.open(skyfile)[0].header.get('CD1_1')
    skyxx = arange(len(skyff))
    skyaa = crval1 + (skyxx) * cd1
    shift = ntt.efoscspec2Ddef.checkwavelength_arc(
        aa1, yy1, skyaa, skyff, '', '')
    delete('_new3.fits')
    return shift


def findsubimage(imglist):
    # print "LOGX:: Entering `findsubimage` method/function in %(__file__)s" %
    # globals()
    import ntt
    from ntt.util import readhdr, readkey3
    from numpy import abs, argmin
    import string
    import os
    import sys
    import re

    imgsub = []
    for img in imglist:
        imglist2 = imglist[:]
        imglist2.remove(img)
        hdr = readhdr(img)
        JD0 = readkey3(hdr, 'JD')
        xcum0 = readkey3(hdr, 'xcum')
        distance = []
        imglist3 = []
        for img2 in imglist2:
            if abs(readkey3(readhdr(img2), 'xcum') - xcum0) > 10:
                distance.append(abs(JD0 - readkey3(readhdr(img2), 'JD')))
                imglist3.append(img2)
        if len(distance) > 0:
            imgsub0 = imglist3[argmin(distance)]
        else:
            imgsub0 = ''
        imgsub.append(imgsub0)
    return imgsub


def sofispecreduction(files, _interactive, _doflat, listflat, _docross, _verbose=False):
    # print "LOGX:: Entering `sofispecreduction` method/function in
    # %(__file__)s" % globals()
    import ntt
    from ntt.util import delete, readhdr, readkey3, correctcard, rangedata
    import string, re, sys, os, glob

    try:        import pyfits
    except:     from astropy.io import fits as pyfits

    from pyraf import iraf
    from numpy import argmin, array, min, isnan, arange, mean, sum
    from numpy import sqrt, pi

    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.ccdred(_doprint=0)
    iraf.twodspec(_doprint=0)
    iraf.longslit(_doprint=0)
    iraf.specred(_doprint=0)
    toforget = ['ccdred.flatcombine', 'ccdproc', 'specred.apall', 'longslit.identify', 'longslit.reidentify',
                'longslit.fitcoords', 'specred.transform', 'specred.response', 'imutil.hedit']
    for t in toforget:
        iraf.unlearn(t)
    iraf.longslit.dispaxi = 2
    iraf.longslit.mode = 'h'
    iraf.specred.dispaxi = 2
    iraf.specred.mode = 'h'
    iraf.ccdproc.darkcor = 'no'
    iraf.ccdproc.fixpix = 'no'
    iraf.ccdproc.flatcor = 'no'
    iraf.ccdproc.zerocor = 'no'
    iraf.ccdproc.overscan = 'no'
    iraf.ccdproc.ccdtype = ''
    iraf.ccdred.instrument = "/dev/null"

    iraf.set(direc=ntt.__path__[0] + '/')

    if _interactive:
        _interact = 'yes'
    else:
        _interact = 'no'
    if _verbose:
        iraf.ccdred.verbose = 'yes'
        iraf.specred.verbose = 'yes'
    else:
        iraf.specred.verbose = 'no'
        iraf.ccdred.verbose = 'no'
    import datetime
    import time

    now = datetime.datetime.now()
    datenow = now.strftime('20%y%m%d%H%M')
    MJDtoday = 55927 + (datetime.date.today() - datetime.date(2012, 01, 01)).days
    # if they are not sorted the fieldlist dict could crash
    files = ntt.sofiphotredudef.sortbyJD(files)
    outputlist = []
    setup = []
    fieldlist = {}
    OBID = {}
    RA = {}
    DEC = {}
    objects = {}
    flats = {}
    lamps1 = {}
    _rdnoise = readkey3(readhdr(re.sub('\n', '', files[0])), 'ron')
    _gain = readkey3(readhdr(re.sub('\n', '', files[0])), 'gain')
    for img in files:
        img = re.sub('\n', '', img)
        hdr = readhdr(img)
        _object = readkey3(hdr, 'object')
        _filter = readkey3(hdr, 'filter')
        _date = readkey3(hdr, 'date-night')
        _exptime = readkey3(hdr, 'exptime')
        _grism = readkey3(hdr, 'grism')
        _obsmode = readkey3(hdr, 'obsmode')
        _type = ''
        if _grism.lower() not in ['gr', 'gb']:
            _type = 'image'
        if not _type:
            if _object.lower() == 'flat':
                _type = 'flat'
                if _date not in flats:
                    flats[_date] = {}
                if _grism not in flats[_date]:
                    flats[_date][_grism] = [img]
                else:
                    flats[_date][_grism].append(img)
            elif _object.lower() == 'lamp':
                _lampid = (readkey3(hdr, 'esoid'), readkey3(hdr, 'grism'))
                if _lampid not in lamps1:
                    lamps1[_lampid] = [None, None]
                if readkey3(hdr, 'lamp1') == 'Xenon':
                    lamps1[_lampid][0] = img
                else:
                    lamps1[_lampid][1] = img
                _type = 'lamp'
                # if readkey3(hdr,'lamp1')=='Xenon':
            #                     _type='lamp'
            #                     if _grism not in lamps:
            #                         lamps[_grism]=[img]
            #                     else:
            #                         lamps[_grism].append(img)
            #                 else:
            #                     _type='notgood'
        if not _type:
            _ra = readkey3(hdr, 'RA')
            _dec = readkey3(hdr, 'DEC')
            _object_name = readkey3(hdr, 'object')
            _OBID = (readkey3(hdr, 'esoid'), _grism)
            if string.count(_object_name, '/') or string.count(_object_name, '.') or string.count(_object_name, ' '):
                nameobj = string.split(_object_name, '/')[0]
                nameobj = string.split(nameobj, ' ')[0]
                nameobj = string.split(nameobj, '.')[0]
            else:
                nameobj = _object_name
            if _grism not in fieldlist:
                fieldlist[_grism] = {}
            if _OBID not in OBID:
                count = 1
                nameobj0 = nameobj + '_' + str(count)
                answ = 'yes'
                while answ == 'yes':
                    if nameobj0 in fieldlist[_grism]:
                        count = count + 1
                        nameobj0 = nameobj + '_' + str(count)
                    else:
                        answ = 'no'
                fieldlist[_grism][nameobj0] = []
                OBID[readkey3(hdr, 'esoid'), _grism] = nameobj0
            fieldlist[_grism][nameobj0].append(img)

        if _verbose:
            print img
            print _type, _object, _filter
            print 'lamps', lamps1

    lamps = {}
    for _lampid in lamps1:
        lamp = ''
        output = 'arc_' + str(_lampid[0]) + '_' + str(_lampid[1]) + '.fits'
        if lamps1[_lampid][0] and lamps1[_lampid][1]:
            print lamps1[_lampid][0], lamps1[_lampid][1]
            # try:
            ntt.util.delete(output)
            iraf.imarith(lamps1[_lampid][0], '-', lamps1[_lampid]
                         [1], result=output, verbose='yes')
            #            except:
            #                print 'warning, lamp file not ON/OFF'
            #                os.system('cp '+lamps1[_lampid][0]+' '+output)

            lamp = output
        elif lamps1[_lampid][0] and not lamps1[_lampid][1]:
            os.system('cp ' + lamps1[_lampid][0] + ' ' + output)
            lamp = output
        if lamp:
            if _lampid[1] not in lamps:
                lamps[_lampid[1]] = [lamp]
            else:
                lamps[_lampid[1]].append(lamp)

    if _verbose:
        print '\n### FIELDS\n', fieldlist
        print '\n### OBID\n', OBID
        print '\n### FLATS\n', flats
        print '\n### LAMPS\n', lamps

#    if not flats:
#        sys.exit('\n### error: spectroscopic flat not available, add flats in the directory and try again')
#    if not lamps:
# sys.exit('\n### error: spectroscopic lamp not available, add lamps in
# the directory and try again')

    if not listflat:
        print '\n### list of available spectroscopic flats (ON,OFF):'
        for _date in flats:
            for _grism in flats[_date]:
                for img in flats[_date][_grism]:
                    if pyfits.open(img)[0].data.mean() >= 2000:
                        print img, _grism, _date, 'ON ? '
                    else:
                        print img, _grism, _date, 'OFF ? '
        for _date in flats:
            for _grism in flats[_date]:
                flat = {'ON': [], 'OFF': []}
                for img in flats[_date][_grism]:
                    _type = ''
                    if readkey3(hdr, 'lamp3'):
                        print '\n### header lamp3 found: flat ON ', str(img)
                        _type = 'ON'
                    else:
                        if pyfits.open(img)[0].data.mean() >= 2000:
                            _type = 'ON'
                        else:
                            _type = 'OFF'
                    aa, bb, cc = ntt.util.display_image(img, 1, '', '', False)
                    print '\n### number of flat already selected (ON,OFF): \n ### please select same number ' \
                          'of ON and OFF flats \n' + \
                        str(len(flat['ON'])) + '  ' + str(len(flat['OFF']))
                    print '\n### image ' + str(img)
                    answ = raw_input(
                        'ON/OFF/REJECT/STOP [' + str(_type) + ']  ok (ON[n]/OFF[f]/r/s) [' + _type + '] ? ')
                    if not answ:
                        answ = _type
                    if answ in ['ON', 'on', 'n']:
                        _type = 'ON'
                    if answ in ['OFF', 'off', 'f']:
                        _type = 'OFF'
                    if answ in ['s', 'S', 'STOP', 'stop', 'Stop']:
                        _type = 'stop'
                    if answ in ['r', 'R', 'reject']:
                        _type = 'r'
                    if _type in ['ON', 'OFF']:
                        flat[_type].append(img)
                    elif _type == 'stop':
                        if len(flat['ON']) == len(flat['OFF']) and len(flat['OFF']) >= 2:
                            break
                        elif len(flat['ON']) == len(flat['OFF']) and len(flat['OFF']) == 0:
                            break
                        else:
                            print '\n### Warning: you can stop only if the numbers of ON and OFF are the same'
                print len(flat['ON']), len(flat['OFF'])
                if len(flat['ON']) == len(flat['OFF']) and len(flat['OFF']) >= 2:
                    ff = open('_flatlist', 'w')
                    for ii in range(0, len(flat['OFF'])):
                        delete('flat_' + str(_date) + '_' + str(_grism) +
                               '_' + str(MJDtoday) + '_' + str(ii) + '.fits')
                        iraf.imarith(flat['ON'][ii], '-', flat['OFF'][ii],
                                     result='flat_' + str(_date) + '_' + str(_grism) + '_' + str(MJDtoday) + '_' + str(
                                         ii) + '.fits', verbose='no')
                        ff.write(
                            'flat_' + str(_date) + '_' + str(_grism) + '_' + str(MJDtoday) + '_' + str(ii) + '.fits\n')
                    ff.close()
                    masterflat = 'flat_' + \
                        str(_date) + '_' + str(_grism) + \
                        '_' + str(MJDtoday) + '.fits'
                    delete(masterflat)
                    _order = '80'
                    iraf.ccdred.flatcombine(input='@_flatlist', output=masterflat, combine='median', rdnoise=_rdnoise,
                                            gain=_gain, ccdtype='')
                    hdr = readhdr(masterflat)
                    matching = [s for s in hdr.keys() if "IMCMB" in s]
                    for imcmb in matching:
                        aaa = iraf.hedit(masterflat, imcmb, delete='yes', update='yes',
                                         verify='no', Stdout=1)
                    delete('_flatlist')
                    print masterflat
                    correctcard(masterflat)
                    if masterflat not in outputlist:
                        outputlist.append(masterflat)
                    ntt.util.updateheader(masterflat, 0, {'FILETYPE': [41102, 'flat field'],
                                                          'SINGLEXP': [False, 'TRUE if resulting from single exposure'],
                                                          'M_EPOCH': [False, 'TRUE if resulting from multiple epochs']})

                    print '\n###  master flat ........... done '
                    delete('n' + masterflat)
                    iraf.specred.response(masterflat, normaliz=masterflat + '[100:900,*]',
                                          response='n' + masterflat, interac=_interact, thresho='INDEF', sample='*',
                                          naverage=2,
                                          function='spline3', low_rej=3, high_rej=3, order=_order, niterat=20, grow=0,
                                          graphic='stdgraph', mode='q')
                    listflat.append('n' + masterflat)
                    if 'n' + masterflat not in outputlist:
                        outputlist.append('n' + masterflat)
                    ntt.util.updateheader('n' + masterflat, 0, {'FILETYPE': [41203, 'normalized flat field'],
                                                                'TRACE1': [masterflat, 'Originating file']})
                    # ntt.util.updateheader('n'+masterflat,0,{'TRACE1':[masterflat,'']})

                    flattot = flat['ON'] + flat['OFF']
                    num = 0
                    for img in flattot:
                        num = num + 1
                        ntt.util.updateheader(masterflat, 0, {
                            'PROV' + str(num): [readkey3(readhdr(img), 'ARCFILE'), 'Originating file'],
                            'TRACE' + str(num): [readkey3(readhdr(img), 'ARCFILE'), 'Originating file']})
                        ntt.util.updateheader('n' + masterflat, 0, {
                            'PROV' + str(num): [readkey3(readhdr(img), 'ARCFILE'), 'Originating file']})

                    if listflat:
                        print '\n### flat available:\n### ' + str(listflat), '\n'
                elif len(flat['ON']) == len(flat['OFF']) and len(flat['OFF']) == 0:
                    print '\n### no good flats in this set ......'
                else:
                    sys.exit('\n### Error: number of ON and OFF not the same')

    for _grism in fieldlist:
        obj0 = fieldlist[_grism][fieldlist[_grism].keys()[0]][0]
        # #############              arc              #########################
        if _grism not in lamps:
            print '\n### take arc from archive '
            arcfile = ntt.util.searcharc(obj0, '')[0]
            if arcfile[0] == '/':
                os.system('cp ' + arcfile + ' ' +
                          string.split(arcfile, '/')[-1])
                arcfile = string.split(arcfile, '/')[-1]
            lamps[_grism] = [arcfile]

        if _grism in lamps:
            arclist = lamps[_grism]
            if arclist:
                arcfile = ntt.util.searcharc(obj0, arclist)[0]
            else:
                arcfile = ntt.util.searcharc(obj0, '')[0]

            print arcfile
            if arcfile:
                print arcfile
                datea = readkey3(readhdr(arcfile), 'date-night')
                if arcfile[0] == '/':
                    os.system('cp ' + arcfile + ' ' +
                              string.split(arcfile, '/')[-1])
                    arcfile = string.split(arcfile, '/')[-1]

                if _doflat:
                    if listflat:
                        flat0 = ntt.util.searchflat(arcfile, listflat)[0]
                    else:
                        flat0 = ''
                else:
                    flat0 = ''

                if flat0:
                    _flatcor = 'yes'
                else:
                    _flatcor = 'no'
                    _doflat = False

                ntt.util.delete('arc_' + datea + '_' + _grism +
                                '_' + str(MJDtoday) + '.fits')

                print arcfile, flat0, _flatcor, _doflat

                if _doflat:
                    iraf.noao.imred.ccdred.ccdproc(arcfile,
                                                   output='arc_' + datea + '_' + _grism +
                                                   '_' +
                                                   str(MJDtoday) + '.fits',
                                                   overscan='no', trim='no', zerocor='no', flatcor=_flatcor, flat=flat0)
                else:
                    os.system('cp ' + arcfile + ' ' + 'arc_' + datea +
                              '_' + _grism + '_' + str(MJDtoday) + '.fits')

                iraf.noao.imred.ccdred.ccdproc('arc_' + datea + '_' + _grism + '_' + str(MJDtoday) + '.fits', output='',
                                               overscan='no', trim='yes', zerocor='no', flatcor='no', flat='',
                                               trimsec='[30:1000,1:1024]')

                arcfile = 'arc_' + datea + '_' + \
                    _grism + '_' + str(MJDtoday) + '.fits'

                ntt.util.correctcard(arcfile)
                print arcfile

                if arcfile not in outputlist:
                    outputlist.append(arcfile)

                ntt.util.updateheader(arcfile, 0, {'FILETYPE': [41104, 'pre-reduced 2D arc'],
                                                   'SINGLEXP': [True, 'TRUE if resulting from single exposure'],
                                                   'M_EPOCH': [False, 'TRUE if resulting from multiple epochs'],
                                                   'PROV1': [readkey3(readhdr(arcfile), 'ARCFILE'), 'Originating file'],
                                                   'TRACE1': [readkey3(readhdr(arcfile), 'ARCFILE'),
                                                              'Originating file']})

                arcref = ntt.util.searcharc(obj0, '')[0]
                if not arcref:
                    identific = iraf.longslit.identify(images=arcfile, section='column 10',
                                                       coordli='direc$standard/ident/Lines_XeAr_SOFI.dat', nsum=10,
                                                       fwidth=7, order=3, mode='h', Stdout=1, verbose='yes')
                else:
                    print arcref
                    os.system('cp ' + arcref + ' .')
                    arcref = string.split(arcref, '/')[-1]
                    if not os.path.isdir('database/'):
                        os.mkdir('database/')
                    if os.path.isfile(ntt.util.searcharc(obj0, '')[1] + '/database/id' + re.sub('.fits', '', arcref)):
                        os.system('cp ' + ntt.util.searcharc(obj0, '')[1] + '/database/id' + re.sub('.fits', '',
                                                                                                    arcref) + ' database/')

                    print arcref, arcfile
                    #                        time.sleep(5)
                    #                        os.system('rm -rf database/idarc_20130417_GR_56975')
                    #                        raw_input('ddd')
                    identific = iraf.longslit.reidentify(referenc=arcref, images=arcfile, interac='NO',  # _interact,
                                                         section='column 10', shift=0.0,
                                                         coordli='direc$standard/ident/Lines_XeAr_SOFI.dat',
                                                         overrid='yes', step=0, newaps='no', nsum=5, nlost=2,
                                                         mode='h', verbose='yes', Stdout=1)
                    #                        print identific
                    #                        raw_input('ddd')
                    identific = iraf.longslit.reidentify(referenc=arcref, images=arcfile, interac=_interact,
                                                         section='column 10', shift=1.0,
                                                         coordli='direc$standard/ident/Lines_XeAr_SOFI.dat',
                                                         overrid='yes', step=0, newaps='no', nsum=5, nlost=2,
                                                         mode='h', verbose='yes', Stdout=1)
                    #                        fitsfile = ntt.efoscspec2Ddef.continumsub('new3.fits', 6, 1)
                    # I need to run twice I don't know why
                    #                        print identific
                    #                        raw_input('ddd')
                    if _interactive:
                        answ = raw_input(
                            '\n### do you like the identification [[y]/n]')
                        if not answ:
                            answ = 'y'
                    else:
                        answ = 'y'
                    if answ in ['n', 'N', 'no', 'NO', 'No']:
                        yy1 = pyfits.open(arcref)[0].data[:, 10:20].mean(1)
                        xx1 = arange(len(yy1))
                        yy2 = pyfits.open(arcfile)[0].data[:, 10:20].mean(1)
                        xx2 = arange(len(yy2))

                        ntt.util.delete('_new3.fits')
                        hdu = pyfits.PrimaryHDU(yy1)
                        hdulist = pyfits.HDUList([hdu])
                        hdulist.writeto('_new3.fits')

                        fitsfile = ntt.efoscspec2Ddef.continumsub('_new3.fits', 4, 1)
                        yy1 = pyfits.open(fitsfile)[0].data

                        ntt.util.delete('_new3.fits')
                        hdu = pyfits.PrimaryHDU(yy2)
                        hdulist = pyfits.HDUList([hdu])
                        hdulist.writeto('_new3.fits')

                        fitsfile = ntt.efoscspec2Ddef.continumsub('_new3.fits', 4, 1)
                        yy2 = pyfits.open(fitsfile)[0].data

                        _shift = ntt.efoscspec2Ddef.checkwavelength_arc(
                            xx1, yy1, xx2, yy2, '', '') * (-1)

                        print arcref, arcfile, _shift
                        identific = iraf.longslit.reidentify(referenc=arcref, images=arcfile, interac='YES',
                                                             section='column 10', shift=_shift,
                                                             coordli='direc$standard/ident/Lines_XeAr_SOFI.dat',
                                                             overrid='yes', step=0, newaps='no', nsum=5, nlost=2,
                                                             mode='h', verbose='yes', Stdout=1)

                        answ = raw_input('\n### is it ok now ? [[y]/n] ')
                        if not answ:
                            answ = 'y'
                        if answ in ['n', 'N', 'no', 'NO', 'No']:
                            sys.exit(
                                '\n### Warning: line identification with some problems')
                iraf.longslit.reidentify(referenc=arcfile, images=arcfile, interac='NO', section='column 10',
                                         coordli='direc$standard/ident/Lines_XeAr_SOFI.dat', overrid='yes', step=10,
                                         newaps='yes', nsum=5, nlost=2, mode='h', verbose='no')
                iraf.longslit.fitcoords(images=re.sub('.fits', '', arcfile), fitname=re.sub('.fits', '', arcfile),
                                        interac='no', combine='yes', databas='database',
                                        function='legendre', yorder=4, logfile='', plotfil='', mode='h')
                if identific:
                    _rms = float(identific[-1].split()[-1])
                    _num = float(identific[-1].split()[2].split('/')[0])
                    hdr = ntt.util.readhdr(arcfile)
                    hedvec = {'LAMRMS': [_rms * .1, 'residual RMS [nm]'],
                              'LAMNLIN': [_num, 'Nb of arc lines used in the fit of the wavel. solution'],
                              'SPEC_ERR': [(_rms * .1) / sqrt(float(_num)), 'statistical uncertainty'],
                              'SPEC_SYE': [0.1, 'systematic error']}
                    ntt.util.updateheader(arcfile, 0, hedvec)
            else:
                sys.exit('Warning: arcfile not found')
        else:
            print 'here'
        # ########################################################################################################
        for field in fieldlist[_grism]:
            listaobj = fieldlist[_grism][field]
            listaobj = ntt.sofiphotredudef.sortbyJD(listaobj)
            listatemp = listaobj[:]
            # ##############             flat            ######################
            if listflat and _doflat:
                flat0 = ntt.util.searchflat(listaobj[0], listflat)[0]
            else:
                flat0 = ''
            if flat0:
                _flatcor = 'yes'
            else:
                _flatcor = 'no'

            ##########   crosstalk        ###########################

            listatemp2 = []
            _date = readkey3(readhdr(listatemp[0]), 'date-night')
            for img in listatemp:
                #                    num2=listatemp.index(listasub[j])
                imgout = field + '_' + str(_date) + '_' + str(_grism) + '_' + str(MJDtoday) + '_' + str(
                    listatemp.index(img)) + '.fits'
                print '\n### input image: ' + str(img)
                delete(imgout)
                listatemp2.append(imgout)
                if _docross:
                    print '### correct for cross talk   .....   done'
                    ntt.sofiphotredudef.crosstalk(img, imgout)
                    correctcard(imgout)
                    ntt.util.updateheader(
                        imgout, 0, {'CROSSTAL': ['True', '']})
                else:
                    os.system('cp ' + img + ' ' + imgout)
                    correctcard(imgout)
                if _flatcor == 'yes':
                    print '### correct for flat field   .....   done'
                    try:
                        iraf.noao.imred.ccdred.ccdproc(imgout, output='', overscan='no', trim='no', zerocor='no',
                                                       flatcor=_flatcor, flat=flat0)
                    except:
                        iraf.imutil.imreplace(
                            images=flat0, value=0.01, lower='INDEF', upper=0.01, radius=0)
                        iraf.noao.imred.ccdred.ccdproc(imgout, output='', overscan='no', trim='no', zerocor='no',
                                                       flatcor=_flatcor, flat=flat0)
                iraf.noao.imred.ccdred.ccdproc(imgout, output='', overscan='no', trim='yes', zerocor='no',
                                               flatcor='no', flat='', trimsec='[30:1000,1:1024]')
                ntt.util.updateheader(
                    imgout, 0, {'FLATCOR': [flat0, 'flat correction']})

                if imgout not in outputlist:
                    outputlist.append(imgout)
                ntt.util.updateheader(imgout, 0, {'FILETYPE': [42104, 'pre-reduced frame'],
                                                  'SINGLEXP': [True, 'TRUE if resulting from single exposure'],
                                                  'M_EPOCH': [False, 'TRUE if resulting from multiple epochs'],
                                                  'PROV1': [readkey3(readhdr(imgout), 'ARCFILE'), 'Originating file'],
                                                  'TRACE1': [readkey3(readhdr(imgout), 'ARCFILE'), 'Originating file']})
                print '### output image: ' + str(imgout)

            listatemp = listatemp2[:]
            #########    differences object images  #####################
            listasub = ntt.sofispec2Ddef.findsubimage(listatemp)
            reduced = []
            print '\n### Select Frames to be subtracted (eg A-B, B-A, C-D, D-C, ....) '
            print '###    frame1 \t  frame2  \t   offset1  \t   offset2  \t  JD1  \t    JD2\n'
            if len(listatemp) >= 2 and len(listasub) >= 2:
                for j in range(0, len(listatemp)):
                    print '### ', listatemp[j], listasub[j], str(readkey3(readhdr(listatemp[j]), 'xcum')), str(
                        readkey3(readhdr(listasub[j]), 'xcum')), \
                        str(readkey3(readhdr(listatemp[j]), 'JD')), str(
                            readkey3(readhdr(listatemp[j]), 'JD'))
                    if _interactive:
                        answ = raw_input('\n### ok [[y]/n] ? ')
                        if not answ:
                            answ = 'y'
                    else:
                        answ = 'y'
                    num1 = j
                    image1 = listatemp[j]
                    _date = readkey3(readhdr(image1), 'date-night')
                    if answ == 'y':
                        num2 = listatemp.index(listasub[j])
                        image2 = listasub[j]
                    else:
                        image2 = raw_input(
                            'which image do you want to subtract')
                        num2 = listatemp.index(image2)
                    imgoutsub = field + '_' + str(_date) + '_' + str(_grism) + '_' + str(MJDtoday) + '_' + str(
                        num1) + '_' + str(num2) + '.fits'
                    delete(imgoutsub)
                    iraf.images.imutil.imarith(
                        operand1=image1, op='-', operand2=image2, result=imgoutsub, verbose='no')
                    ntt.util.updateheader(imgoutsub, 0, {'skysub': [image2, 'sky image subtracted'],
                                                         'FILETYPE': [42115, 'pre-reduced frame sky subtracted'],
                                                         'TRACE1': [image1, 'Originating file'],
                                                         'PROV2': [readkey3(readhdr(image2), 'ARCFILE'),
                                                                   'Originating file'],
                                                         'TRACE2': [image2, 'Originating file']})

                    reduced.append(imgoutsub)
                    if imgoutsub not in outputlist:
                        outputlist.append(imgoutsub)
            ########################     2D wavelengh calibration      ########
            for img in reduced:
                if arcfile:
                    hdra = ntt.util.readhdr(arcfile)
                    delete('t' + img)
                    iraf.specred.transform(input=img, output='t' + img, minput='',
                                           fitnames=re.sub('.fits', '', arcfile), databas='database',
                                           x1='INDEF', x2='INDEF', y1='INDEF', y2='INDEF', flux='yes', mode='h',
                                           logfile='logfile')
                    ntt.util.updateheader('t' + img, 0,
                                          {'ARC': [arcfile, ''], 'FILETYPE': [42106, 'wavelength calibrate 2D frames'],
                                           'TRACE1': [img, 'Originating file']})
                    ntt.util.updateheader(
                        't' + img, 0, {'TRACE1': [img, 'Originating file']})
                    ntt.util.updateheader('t' + img, 0,
                                          {'LAMRMS': [ntt.util.readkey3(hdra, 'LAMRMS'), 'residual RMS [nm]'],
                                           'LAMNLIN': [ntt.util.readkey3(hdra, 'LAMNLIN'), 'number of arc lines'],
                                           'SPEC_ERR': [ntt.util.readkey3(hdra, 'SPEC_ERR'), 'statistical uncertainty'],
                                           'SPEC_SYE': [ntt.util.readkey3(hdra, 'SPEC_SYE'), 'systematic error']})
                    ###########################
                    delete('t' + arcfile)
                    iraf.specred.transform(input=arcfile, output='t' + arcfile, minput='',
                                           fitnames=re.sub('.fits', '', arcfile), databas='database',
                                           x1='INDEF', x2='INDEF', y1='INDEF', y2='INDEF', flux='yes', mode='h',
                                           logfile='logfile')
                    specred = ntt.util.spectraresolution2(arcfile, 50)
                    if specred:
                        ntt.util.updateheader(
                            't' + img, 0, {'SPEC_RES': [specred, 'Spectral resolving power']})
                    delete('t' + arcfile)
                    ###########################
                    iraf.hedit('t' + img, 'TRACE2', delete='yes',
                               update='yes', verify='no', Stdout=1)

                    if 't' + img not in outputlist:
                        outputlist.append('t' + img)
                    print '\n### 2D frame t' + str(img) + ' wavelengh calibrated  ............ done'

                    _skyfile = ntt.__path__[
                        0] + '/standard/ident/sky_' + _grism + '.fits'  # check in wavelengh   #########
                    hdr = ntt.util.readhdr(img)
                    if glob.glob(_skyfile) and readkey3(hdr, 'exptime') > 20.:
                        _original = readkey3(hdr, 'ORIGFILE')
                        _archive = readkey3(hdr, 'ARCFILE')
                        if os.path.isfile(_archive):
                            imgstart = _archive
                        elif os.path.isfile(_original):
                            imgstart = _original
                        else:
                            imgstart = ''
                        if imgstart:
                            delete('_tmp.fits')
                            print imgstart, arcfile
                            iraf.specred.transform(input=imgstart, output='_tmp.fits', minput='',
                                                   fitnames=re.sub('.fits', '', arcfile), databas='database',
                                                   x1='INDEF', x2='INDEF', y1='INDEF', y2='INDEF', flux='yes', mode='h',
                                                   logfile='logfile')

                            shift = ntt.sofispec2Ddef.skysofifrom2d('_tmp.fits', _skyfile)
                            zro = pyfits.open('_tmp.fits')[0].header.get('CRVAL2')

                            delete('_tmp.fits')
                            if _interactive:
                                answ = raw_input(
                                    'do you want to correct the wavelengh calibration with this shift: ' + str(
                                        shift) + ' [[y]/n] ? ')
                                if not answ:
                                    answ = 'y'
                            else:
                                answ = 'y'
                            if answ.lower() in ['y', 'yes']:
                                ntt.util.updateheader('t' + img, 0,
                                                      {'CRVAL2': [zro + int(shift), ''], 'shift': [float(shift), '']})
                            #                                    ntt.util.updateheader('t'+img,0,{'shift':[float(shift),'']})
                            print '\n### check wavelengh calibration with sky lines ..... done'
                    try:
                        hdrt = ntt.util.readhdr('t' + img)
                        wavelmin = float(readkey3(hdrt, 'CRVAL2')) + (0.5 - float(readkey3(hdrt, 'CRPIX2'))) * float(
                            readkey3(hdrt, 'CDELT2'))
                        wavelmax = float(readkey3(hdrt, 'CRVAL2')) + (
                            (float(readkey3(hdrt, 'NAXIS2')) + 0.5 - float(readkey3(hdrt, 'CRPIX2'))) * float(
                                readkey3(hdrt, 'CDELT2')))
                        hedvec = {}
                        hedvec['WAVELMIN'] = [
                            wavelmin * .1, '[nm] minimum wavelength']
                        hedvec['WAVELMAX'] = [
                            wavelmax * .1, ' [nm] maximum wavelength']
                        hedvec['XMIN'] = [wavelmin, '[A] minimum wavelength']
                        hedvec['XMAX'] = [wavelmax, '[A]  maximum wavelength']
                        hedvec['SPEC_BW'] = [
                            (wavelmax * .1) - (wavelmin * .1), '[nm] Bandpass Width Wmax - Wmin']
                        hedvec['SPEC_VAL'] = [
                            ((wavelmax * .1) + (wavelmin * .1)) / 2., '[nm] Mean Wavelength']
                        hedvec['SPEC_BIN'] = [
                            ((wavelmax * .1) - (wavelmin * .1)) /
                            (float(readkey3(hdr, 'NAXIS2')) - 1),
                            'Wavelength bin size [nm/pix]']
                        hedvec['VOCLASS'] = ['SPECTRUM V1.0', 'VO Data Model']
                        hedvec['VOPUB'] = ['ESO/SAF',
                                           'VO Publishing Authority']
                        #                            hedvec['APERTURE']=[float(re.sub('slit','',readkey3(hdrt,'slit'))),'aperture width']
                        ntt.util.updateheader('t' + img, 0, hedvec)
                    except:
                        pass
                else:
                    print '\n### Warning: arc not found for the image ' + str(img) + ' with setup ' + str(_grism)

    reduceddata = rangedata(outputlist)
    print '\n### adding keywords for phase 3 ....... '
    f = open('logfile_spec2d_' + str(reduceddata) +
             '_' + str(datenow) + '.raw.list', 'w')
    for img in outputlist:
        if img[-4:] == 'fits':
            hdr = readhdr(img)
            # ###############################################
            # cancel pc matrix
            if 'PC1_1' in hdr.keys():
                aaa = iraf.hedit(img, 'PC1_1', delete='yes',
                                 update='yes', verify='no', Stdout=1)
            if 'PC2_2' in hdr.keys():
                aaa = iraf.hedit(img, 'PC2_2', delete='yes',
                                 update='yes', verify='no', Stdout=1)
            if 'PC1_2' in hdr.keys():
                aaa = iraf.hedit(img, 'PC1_2', delete='yes',
                                 update='yes', verify='no', Stdout=1)
            if 'PC2_1' in hdr.keys():
                aaa = iraf.hedit(img, 'PC2_1', delete='yes',
                                 update='yes', verify='no', Stdout=1)
            #################
            # added for DR2
            print img

            if 'NCOMBINE' in hdr:
                _ncomb = readkey3(hdr, 'NCOMBINE')
            else:
                _ncomb = 1.0

            ntt.util.updateheader(
                img, 0, {'DETRON ': [12, 'Readout noise per output (e-)']})
            ntt.util.updateheader(img, 0, {'EFFRON': [12. * (1 / sqrt(readkey3(hdr, 'ndit') * _ncomb)) * sqrt(pi / 2),
                                                      'Effective readout noise per output (e-)']})
            ntt.util.phase3header(img)  # phase 3 definitions
            ############################
            #  change for DR2
            ############################
            texp = float(readkey3(hdr, 'dit')) * float(readkey3(hdr, 'ndit'))
            mjdend = float(readkey3(hdr, 'MJD-OBS')) + (float(readkey3(hdr, 'ndit')) * (
                float(readkey3(hdr, 'dit')) + 1.8)) / (60. * 60. * 24.)
            strtexp = time.strftime('%H:%M:%S', time.gmtime(texp))
            _telapse = (mjdend - float(readkey3(hdr, 'MJD-OBS'))) * \
                60. * 60 * 24.
            # tmid=_telapse/2.
            tmid = (mjdend + float(readkey3(hdr, 'MJD-OBS'))) / 2
            ntt.util.updateheader(img, 0, {'quality': ['Final', 'fast or rapid reduction'],
                                           'BUNIT': ['ADU', 'Physical unit of array values'],
                                           'DIT': [readkey3(hdr, 'dit'), 'Detector Integration Time'],
                                           'NDIT': [readkey3(hdr, 'ndit'), 'Number of sub-integrations'],
                                           'TEXPTIME': [texp, 'Total integration time of all exposures (s)'],
                                           'EXPTIME': [texp, 'Total integration time. ' + strtexp],
                                           'MJD-END': [mjdend, 'End of observations (days)'],
                                           'TELAPSE': [_telapse, 'Total elapsed time [days]'],
                                           'TMID': [tmid, '[d] MJD mid exposure'],
                                           'TITLE': [readkey3(hdr, 'object'), 'Dataset title'],
                                           #'TITLE':[str(tmid)[0:9]+' '+str(readkey3(hdr,'object'))+' '+str(readkey3(hdr,'grism'))+' '+\
                                           # str(readkey3(hdr,'filter'))+'
                                           # '+str(readkey3(hdr,'slit')),'Dataset
                                           # title'],\
                                           'EXT_OBJ': [False, 'TRUE if extended'],
                                           'CONTNORM': [False, 'spectrum normalized to the continuum'],
                                           'TOT_FLUX': [False, 'TRUE if phot cond and all src flux is captured'],
                                           'SPECSYS': ['TOPOCENT', 'Reference frame for spectral coordinate'],
                                           'FLUXCAL': ['ABSOLUTE', 'type of flux calibration'],
                                           'FLUXERR': [34.7, 'Fractional uncertainty of the flux [%]'],
                                           'DISPELEM': ['Gr#' + re.sub('Gr', '', readkey3(hdr, 'grism')),
                                                        'Dispersive element name']})
            if readkey3(hdr, 'tech'):
                ntt.util.updateheader(
                    img, 0, {'PRODCATG': ['SCIENCE.IMAGE', 'Data product category']})
            aaa = str(readkey3(hdr, 'arcfiles')) + '\n'
            f.write(aaa)
            try:
                ntt.util.airmass(img)  # phase 3 definitions
            except:
                print '\n### airmass not computed for image: ', img
        else:
            print img + ' is not a fits image'
    f.close()
    return outputlist, 'logfile_spec2d_' + str(reduceddata) + '_' + str(datenow) + '.raw.list'
