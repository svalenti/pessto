try:       from astropy.io import fits as pyfits
except:    import pyfits

def efoscreduction(imglist, _interactive, _doflat, _dobias, listflat, listbias, _dobadpixel, badpixelmask,
                   fringingmask, _archive, typefile, filenameobjects, _system, _cosmic, _verbose=False, method='iraf'):
    import ntt
    from ntt.efoscphotredudef import searchbias
    from ntt.util import delete, readhdr, readkey3, display_image, searchflat, rangedata, correctcard
    from numpy import argmin, min, abs, sqrt
    import string, os, re, math, sys
    from pyraf import iraf
    # ##   Call and set parameters for useful iraf tasks
    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.ccdred(_doprint=0)
    iraf.proto(_doprint=0)

    toforget = ['ccdproc', 'zerocombine', 'flatcombine', 'imreplace', 'proto.fixpix']
    for t in toforget:
        iraf.unlearn(t)

    iraf.ccdproc.darkcor = 'no'
    iraf.ccdproc.fixpix = 'no'
    iraf.ccdproc.flatcor = 'no'
    iraf.ccdproc.zerocor = 'no'
    iraf.ccdproc.overscan = 'no'
    iraf.ccdproc.ccdtype = ''
    iraf.ccdproc.biassec = ''
    iraf.ccdred.instrument = "/dev/null"

    if _verbose:
        iraf.ccdred.verbose = 'yes'
    else:
        iraf.ccdred.verbose = 'no'
    import datetime, time
    #      starttime=time.time()
    now = datetime.datetime.now()
    datenow = now.strftime('20%y%m%d%H%M')
    MJDtoday = 55927 + (datetime.date.today() - datetime.date(2012, 01, 01)).days
    outputfile = []
    reduceddata = rangedata(imglist)
    img = re.sub('\n', '', imglist[0])
    hdr = readhdr(img)
    _gain = readkey3(hdr, 'gain')
    _rdnoise = readkey3(hdr, 'ron')
    _instrume = readkey3(hdr, 'instrume')
    _trimsec = '[3:1010,1:1015]'
    biaslist = {}
    flatlist1 = {}
    flatlist2 = {}
    objectlist = {}
    filterlist1 = []
    filterlist2 = []
    for img in imglist:
        _type = ''
        img = re.sub('\n', '', img)
        hdr = readhdr(img)
        _naxis1 = readkey3(hdr, 'NAXIS1')
        _naxis2 = readkey3(hdr, 'NAXIS2')
        if _naxis1 != 1030 or _naxis2 != 1030:
            ntt.util.writeinthelog('image ' + str(img) + ' different dimension =\n', './logNTT.txt')
            _type = 'not good'
        if not _type and readkey3(hdr, 'speed') != 'fastL':
            _type = 'not good'
        if not _type and readkey3(hdr, 'instrume') != 'efosc':
            _type = 'not good'
        _imagetype = readkey3(hdr, 'tech')
        if not _type and _imagetype == 'SPECTRUM':
            _type = 'spectroscopic data'
        if not _type:
            _exptime = readkey3(hdr, 'exptime')
            _date = readkey3(hdr, 'date-night')
            _filter = readkey3(hdr, 'filter')
            if float(_exptime) == 0.0:
                if _date not in biaslist:        biaslist[_date] = []
                biaslist[_date].append(img)
                _type = 'bias'
            if not _type:
                _object = readkey3(hdr, 'object')
                if _filter.lower() in ['g782', 'r784', 'z623', 'u640', 'b639', 'v641', 'r642',
                                       'i705'] and _imagetype == 'IMAGE':
                    if 'sky,flat' in _object.lower():
                        _type = 'flat'
                    elif 'dome' in _object.lower() or 'flat' in _object.lower():
                        _type = 'flat dome'
                    if _type == 'flat':
                        if _filter not in filterlist1:
                            filterlist1.append(_filter)
                            flatlist1[_filter] = []
                        flatlist1[_filter].append(img)
                    if _type == 'flat dome':
                        if _filter not in filterlist2:
                            filterlist2.append(_filter)
                            flatlist2[_filter] = []
                        flatlist2[_filter].append(img)
            if not _type:
                _catg = readkey3(hdr, 'catg')
                if 'science' in _catg.lower() or 'acquisition' in _catg.lower():
                    _type = 'object'
                    if _filter not in objectlist: objectlist[_filter] = []
                    objectlist[_filter].append(img)
                    if 'acquisition' in _catg.lower():
                        try:
                            correctcard(img)
                            _ra1, _dec1, _name = ntt.util.correctobject(img, 'standard_efosc_mab.txt')
                            _ra1, _dec1, _name = ntt.util.correctobject(img, filenameobjects)
                        except:
                            pass

                elif 'focus' in _object.lower():
                    _type = 'not good'
            if not _type:
                print '\n### warning: object not recognized '
                _object = readkey3(hdr, 'object')
                print img, _object, _imagetype
                answ = raw_input('what is it: bias [1], flat [3], object[4], test [5] ?  [5] ')
                if not answ:
                    answ = '5'
                if answ == '1':
                    if _date not in biaslist:
                        biaslist[_date] = ()
                    biaslist[_date].append(img)
                elif answ == '4':
                    if _filter not in objectlist:
                        objectlist[_filter] = []
                    objectlist[_filter].append(img)
                elif answ == '3':
                    tt = raw_input('dome or sky [d/[s]] ? ')
                    if tt == 's':
                        _type = 'flat'
                        _filter = readkey3(hdr, 'filter')
                        if _filter not in filterlist1:
                            filterlist1.append(_filter)
                            flatlist1[_filter] = []
                        flatlist1[_filter].append(img)
                    elif tt == 'd':
                        _type = 'flat dome'
                        _filter = readkey3(hdr, 'filter')
                        if _filter not in filterlist2:
                            filterlist2.append(_filter)
                            flatlist2[_filter] = []
                        flatlist2[_filter].append(img)
                elif answ == '5':
                    _type = 'not good'

    filterlist = list(set(filterlist1 + filterlist2))
    if _verbose:
        print filterlist1
        print filterlist2
        print flatlist1
        print flatlist2
    flatlist = {}
    for _filt in filterlist:
        if _filt not in flatlist1.keys():
            if _filt in flatlist2.keys():
                if len(flatlist2[_filt]) >= 3:
                    flatlist[_filt] = flatlist2[_filt]
        elif len(flatlist1[_filt]) < 3:
            if _filt in flatlist2.keys():
                if len(flatlist2[_filt]) >= 3:
                    flatlist[_filt] = flatlist2[_filt]
        elif _filt in flatlist1.keys():
            if len(flatlist1[_filt]) >= 3:
                flatlist[_filt] = flatlist1[_filt]

    listaout = []
    if _verbose:
        print '\n### flat ', str(flatlist), '\n'
        print '\n### bias ', str(biaslist), '\n'
        print '\n### object ', str(objectlist), '\n'
        ###### masterbias  #################
    if _dobias:
        if not _archive:
            if listbias:
                masterbiaslist = listbias
            else:
                masterbiaslist = []
                if biaslist:
                    for _date in biaslist:
                        print '\n do bias ' + str(_date) + '\n'
                        biaslist[_date] = rejectbias(biaslist[_date], False, 10)
                        if len(biaslist[_date]) >= 3:
                            masterbiasfile = 'bias_' + str(_date) + '_' + str(MJDtoday) + '.fits'
                            delete(masterbiasfile)
                            f = open('biaslist', 'w')
                            h = open('obiaslist', 'w')
                            for img in biaslist[_date]:
                                f.write(img + '\n')
                                h.write('o' + img + '\n')
                                delete('o' + img)
                            f.close()
                            h.close()
                            try:
                                print 'processing bias .....'
                                iraf.ccdproc('@biaslist', output='@obiaslist', overscan="no", trim="yes", zerocor='no',
                                             fixpix='no', ccdtype='', flatcor='no', darkcor='no', biassec='',
                                             trimsec=str(_trimsec), readaxi='column', Stdout=1)
                                iraf.zerocombine('@obiaslist', output=masterbiasfile, combine='median',
                                                 reject='ccdclip', ccdtype='', process='no',
                                                 rdnoise=_rdnoise, gain=_gain, Stdout=1)
                                correctcard(masterbiasfile)
                                num = 0
                                for img in biaslist[_date]:
                                    num = num + 1
                                    ntt.util.updateheader(masterbiasfile, 0, {
                                    'PROV' + str(num): [readkey3(readhdr(img), 'ARCFILE'), 'Originating file']})
                                    ntt.util.updateheader(masterbiasfile, 0, {
                                    'TRACE' + str(num): [readkey3(readhdr(img), 'ARCFILE'), 'Originating file']})
                                    delete('o' + img)
                                ntt.util.updateheader(masterbiasfile, 0,
                                                      {'M_EPOCH': [False, 'TRUE if resulting from multiple epochs']})
                                ntt.util.updateheader(masterbiasfile, 0,
                                                      {'SINGLEXP': [False, 'TRUE if resulting from single exposure']})
                                ntt.util.updateheader(masterbiasfile, 0, {'FILETYPE': [11201, 'bias']})
                                masterbiaslist.append(masterbiasfile)

                                if masterbiasfile not in outputfile:
                                    outputfile.append(masterbiasfile)
                            except:
                                ntt.util.writeinthelog(
                                    'Warning ' + str(biaslist[_date]) + ' problem with this list of bias \n',
                                    './logNTT.txt')
                            if masterbiasfile and _interactive:
                                aa, bb, cc = display_image(masterbiasfile, 1, '', '', False)
                                answ = raw_input('is the masterbias ok [[y]/n] ?')
                                if not answ:
                                    answ = 'y'
                                if answ in ['n', 'no']:
                                    sys.exit('remove bad bias from input list and restart')
        else:
            masterbiaslist = []

    ########## masterflat   #########################

    if _doflat:
        if not _archive:
            if listflat:
                masterflatlist = listflat
            else:
                masterflatlist = []
                if flatlist:
                    for _filter in flatlist:
                        print '\n do flat ' + str(_filter) + '\n'
                        flatlist[_filter] = rejectflat(flatlist[_filter], False)
                        if len(flatlist[_filter]) >= 3:
                            _date = readkey3(readhdr(flatlist[_filter][0]), 'date-night')
                            masterflat = 'flat_' + str(_date) + '_' + str(_filter) + '_' + str(MJDtoday) + '.fits'
                            listaflat = 'flatlist_' + str(_date) + '_' + str(_filter)
                            _bias = ''
                            if masterbiaslist:
                                _bias = searchbias(flatlist[_filter][0], masterbiaslist)[0]
                            if not _bias:
                                _bias = searchbias(flatlist[_filter][0], '')[0]
                            if _bias:
                                if _bias[0] == '/':
                                    os.system('cp ' + _bias + ' .')
                                    _bias = string.split(_bias, '/')[-1]
                                    _zerocor = 'yes'
                                else:
                                    _zerocor = 'yes'
                            else:
                                _zerocor = 'no'
                            answ0 = 'n'
                            while answ0 != 'y':
                                f = open(listaflat, 'w')
                                h = open('o' + listaflat, 'w')
                                for img in flatlist[_filter]:
                                    f.write(img + '\n')
                                    h.write('o' + img + '\n')
                                    delete('o' + img)
                                f.close()
                                h.close()
                                try:
                                    print 'processing flat .....'
                                    iraf.ccdproc('@' + listaflat, output='@o' + listaflat, overscan='no', trim='yes',
                                                 darkcor='no', fixpix='no',
                                                 zerocor=_zerocor, flatcor='no', trimsec=str(_trimsec), biassec='',
                                                 zero=_bias, readaxi='column', ccdtype='', Stdout=1)
                                    delete(masterflat)
                                    iraf.flatcombine('@o' + listaflat, output=masterflat, combine='average',
                                                     reject='avsigclip', ccdtype='', process='no',
                                                     rdnoise=_rdnoise, gain=_gain, statsec='[100:800,100:800]',
                                                     lsigma=3, hsigma=2, Stdout=1)
                                    masterflatlist.append(masterflat)
                                    correctcard(masterflat)
                                    num = 0
                                    for img in flatlist[_filter]:
                                        num = num + 1
                                        ntt.util.updateheader(masterflat, 0, {
                                        'PROV' + str(num): [readkey3(readhdr(img), 'ARCFILE'), 'Originating file']})
                                        ntt.util.updateheader(masterflat, 0, {
                                        'TRACE' + str(num): [readkey3(readhdr(img), 'ARCFILE'), 'Originating file']})
                                        delete('o' + img)
                                    ntt.util.updateheader(masterflat, 0, {'ZEROCOR': [_bias, '']})
                                    ntt.util.updateheader(masterflat, 0, {
                                    'M_EPOCH': [False, 'TRUE if resulting from multiple epochs']})
                                    ntt.util.updateheader(masterflat, 0, {
                                    'SINGLEXP': [False, 'TRUE if resulting from single exposure']})
                                    ntt.util.updateheader(masterflat, 0, {'FILETYPE': [11202, 'flat field']})

                                    if masterflat not in outputfile:
                                        outputfile.append(masterflat)
                                except:
                                    ntt.util.writeinthelog(
                                        'Warning ' + str(flatlist[_filter]) + ' problem with this list of flat \n',
                                        './logNTT.txt')
                                aa, bb, cc = display_image(masterflat, 1, '', '', False)
                                if masterflat and _interactive:
                                    answ = raw_input('is the masterflat ok [[y]/n] ?')
                                    if not answ: answ = 'y'
                                    if answ.lower() in ['n', 'no']:
                                        answ1 = raw_input('try again [[y]/n] ?')
                                        if not answ1: answ1 = 'y'
                                        if answ1.lower() in ['y', 'yes']:
                                            flatlist[_filter] = ntt.efoscphotredudef.rejectflat(flatlist[_filter], True)
                                        else:
                                            sys.exit('error: problem with flat .... exit')
                                    else:
                                        answ0 = 'y'
                                else:
                                    answ0 = 'y'
        else:
            masterflatlist = []
    ##########################################################################
    if len(masterbiaslist) == 0:
        masterbiaslist = ''
    if len(masterflatlist) == 0:
        masterflatlist = ''
    ######################################
    if _verbose:
        print ''
        print '#############################'
        print masterflatlist
        print masterbiaslist
        print '#############################'
        print ''
    if masterflatlist:
        listaout = listaout + masterflatlist
    if masterbiaslist:
        listaout = listaout + masterbiaslist
    if typefile == 'calib':
        objectlist = {}
    for _filter in objectlist:
        for img in objectlist[_filter]:
            hdr = readhdr(img)
            print '\n#####################################################################\n'
            _object = readkey3(hdr, 'object')
            _object = re.sub(' ', '', _object)
            _object = re.sub('/', '_', _object)
            _object = re.sub('\n', '', _object)
            _exptime = readkey3(hdr, 'exptime')
            _date = readkey3(hdr, 'date-night')
            nameout = ntt.util.name_duplicate(img, str(_object) + '_' + str(_date) + '_' + str(_filter) + '_' + str(
                MJDtoday), '')
            _bias = ''
            if _dobias:
                if masterbiaslist:
                    _bias = searchbias(img, masterbiaslist)[0]
                if not _bias:
                    _bias = searchbias(img, '')[0]
            _flat = ''
            if _doflat:
                if masterflatlist:
                    _flat = searchflat(img, masterflatlist)[0]
                if not _flat:
                    _flat = searchflat(img, '')[0]
                if _bias:  ### bias  ###
                    if _bias[0] == '/':
                        os.system('cp ' + _bias + ' .')
                        _bias = string.split(_bias, '/')[-1]
                    _zerocor = 'yes'
                else:
                    _zerocor = 'no'
            else:
                _zerocor = 'no'

            if _flat:  ###  flat  ###
                if _flat[0] == '/':
                    os.system('cp ' + _flat + ' .')
                    _flat = string.split(_flat, '/')[-1]
                _flatcor = 'yes'
            else:
                _flatcor = 'no'
            sss = str(_object) + '_' + str(_date) + '_' + str(_filter)
            print '### input', img, sss
            print '### bias ', _zerocor, _bias
            print '### flat ', _flatcor, _flat
            print '### name ', nameout
            delete(nameout)
            try:
                iraf.ccdproc(img, output=nameout, overscan="no", trim="yes", zerocor=_zerocor, flatcor='no',
                             darkcor='no', trimsec=str(_trimsec), zero=_bias, biassec='', readaxi='column', Stdout=1)
                try:
                    iraf.ccdproc(nameout, output='', overscan="no", trim="no", zerocor='no', flatcor=_flatcor,
                                 darkcor='no', flat=_flat, readaxi='column', ccdtype='', Stdout=1)
                except:
                    iraf.imrepla(images=_flat, value=0.01, lower='INDEF', upper=0.01, radius=0)
                    iraf.ccdproc(nameout, output='', overscan="no", trim="no", zerocor='no', flatcor=_flatcor,
                                 darkcor='no', flat=_flat, readaxi='column', ccdtype='', Stdout=1)
                correctcard(nameout)
                ntt.util.updateheader(nameout, 0, {'FILTER': [readkey3(readhdr(nameout), 'filter'), 'Filter name'],
                                                   'SINGLEXP': [True, 'TRUE if resulting from single exposure'],
                                                   'M_EPOCH': [False, 'TRUE if resulting from multiple epochs'],
                                                   'FLATCOR': [_flat, ''],
                                                   'ZEROCOR': [_bias, ''], 'FILETYPE': [12204, 'pre-reduced image'],
                                                   'PROV1': [readkey3(readhdr(nameout), 'ARCFILE'), 'Originating file'],
                                                   'NCOMBINE': [1, 'Number of raw science data'],
                                                   'TRACE1': [readkey3(readhdr(nameout), 'ARCFILE'),
                                                              'Originating file']})
                ntt.util.airmass(nameout)  #  phase 3 definitions

                ntt.util.writeinthelog('\n', './logNTT.txt')
                ntt.util.writeinthelog('image= ' + str(img) + ' output= ' + str(nameout) + '\n', './logNTT.txt')
                ntt.util.writeinthelog('bias= ' + str(_bias) + ', flat= ' + str(_flat) + '\n', './logNTT.txt')
                ntt.util.writeinthelog('\n', './logNTT.txt')
                if nameout not in outputfile:
                    outputfile.append(nameout)
            except:
                ntt.util.writeinthelog('image ' + str(img) + ' probably corrupted\n', './logNTT.txt')
            if _dobadpixel:
                if not badpixelmask:
                    badpixelmask = 'bad_pixel_mask.fits'
                    delete(badpixelmask)
                    os.system('cp ' + ntt.__path__[0] + '/archive/' + str(
                        _instrume) + '/badpixels/badpixel.fits ' + badpixelmask)
                iraf.proto.fixpix(images=nameout, masks=badpixelmask, linterp='INDEF', cinterp='INDEF', verbose='no')
                ntt.util.updateheader(nameout, 0, {'FIXPIX': [badpixelmask, '']})
                ntt.util.writeinthelog('image ' + str(nameout) + ' bad pixel corrected with ' + badpixelmask + '\n',
                                       './logNTT.txt')
                print '\n### bad pixel mask correction ..... done'
            else:
                ntt.util.writeinthelog('image ' + str(nameout) + ' bad pixel NOT corrected\n', './logNTT.txt')
            if _cosmic:
                try:
                    print '\n### cosmic  ..... '
                    ntt.cosmics.lacos_im(nameout, _output='', gain=_gain, readn=_rdnoise, xorder=9, yorder=9,
                                         sigclip=4.5, sigfrac=0.5, objlim=1, skyval=0, niter=0, verbose=True,
                                         interactive=False)
                    ntt.util.updateheader(nameout, 0, {
                    'LACOSMIC': [True, 'TRUE if Laplacian cosmic ray rejection has been applied to the image']})
                    print '\n### cosmic  .....  removed '
                except Exception, e:
                    print e
            else:
                ntt.util.updateheader(nameout, 0, {
                'LACOSMIC': [False, 'TRUE if Laplacian cosmic ray rejection has been applied to the image']})
            try:
                ##########################
                sexvec = ntt.efoscastrodef.sextractor(nameout)
                for cat in ['2mass', 'usnoa2', 'usnob1']:
                    rmsx3, rmsy3, num3, fwhmgess, ellgess, ccc, rasys3, decsys3, magsat3 = ntt.efoscastrodef.efoscastroloop(
                        [nameout], cat, False, 40, 40, 100, 'rxyscale', 100, 30, sexvec, True, 10, method)
                    if rmsx3 <= 2 and rmsy3 <= 2: break
                if rmsx3 > 2 and rmsy3 > 2:
                    for cat in ['2mass', 'usnoa2', 'usnob1']:
                        rmsx3, rmsy3, num3, fwhmgess, ellgess, ccc, rasys3, decsys3, magsat3 = ntt.efoscastrodef.efoscastroloop(
                            [nameout], cat, False, 20, int(20), int(50), 'rxyscale', 100, 30, sexvec, True, 5, method)
                        if rmsx3 <= 2 and rmsy3 <= 2: break
                    if rmsx3 > 2 and rmsy3 > 2:
                        for cat in ['2mass', 'usnoa2', 'usnob1']:
                            rmsx3, rmsy3, num3, fwhmgess, ellgess, ccc, rasys3, decsys3, magsat3 = ntt.efoscastrodef.efoscastroloop(
                                [nameout], cat, False, int(10), int(10),
                                int(25), 'rxyscale', 100, 30, sexvec, True, int(3), method)
                        ##########################
                astrostring = str(rmsx3) + ' ' + str(rmsy3) + ' ' + str(num3)
                ntt.util.updateheader(nameout, 0, {'ASTROMET': [astrostring, 'rmsx rmsy nstars']})
                print '\n### check astrometry: fine \n### rmsx rmsy nstars: ' + astrostring
            except Exception, e:
                print e
                rmsx3, rmsy3, num3, fwhmgess, ellgess, ccc, rasys3, decsys3, magsat3 = '', '', '', '', '', '', '', '', ''
                print '\n### problem with astrometry, do you have network ? '
            if fwhmgess and fwhmgess < 99:
                ntt.util.updateheader(nameout, 0, {'PSF_FWHM': [fwhmgess, 'Spatial resolution (arcsec)'],
                                                   'ELLIPTIC': [ellgess, 'Average ellipticity of point sources'],
                                                   'CRDER1': [(1 / sqrt(2.)) * float(rmsx3) * (1. / 3600.),
                                                              'Random error (degree)'],
                                                   'CRDER2': [(1 / sqrt(2.)) * float(rmsy3) * (1. / 3600.),
                                                              'Random error (degree)'],
                                                   'CUNIT1': ['deg', 'unit of the coord. trans.'],
                                                   'CUNIT2': ['deg', 'unit of the coord. trans.'],
                                                   'CSYER1': [rasys3, 'Systematic error (RA_m - Ra_ref)'],
                                                   'CSYER2': [decsys3, 'Systematic error (DEC_m - DEC_ref)']})
            else:
                ntt.util.updateheader(nameout, 0, {'PSF_FWHM': [9999., 'FHWM (arcsec) - computed with sectractor'],
                                                   'ELLIPTIC': [9999., 'ellipticity of point sources (1-b/a)'],
                                                   'CRDER1': [9999., 'Random error in axis 1'],
                                                   'CRDER2': [9999., 'Random error in axis 2'],
                                                   'CUNIT1': ['deg', 'unit of the coord. trans.'],
                                                   'CUNIT2': ['deg', 'unit of the coord. trans.'],
                                                   'CSYER1': [9999., 'Systematic error (RA_m - Ra_ref)'],
                                                   'CSYER2': [9999., 'Systematic error (DEC_m - DEC_ref)']})

            try:
                result = ntt.efoscastrodef.zeropoint(nameout, _system, False, False)
            except:
                result = ''
            if result:
                if os.path.isfile(re.sub('.fits', '.ph', nameout)):
                    if re.sub('.fits', '.ph', nameout) not in outputfile: outputfile.append(
                        re.sub('.fits', '.ph', nameout))
                print '\n### zeropoint ..... done'
                for ll in result:
                    valore = '%3.3s %6.6s %6.6s' % (str(ll), str(result[ll][1]), str(result[ll][0]))
                    print '### ', valore
                    ntt.util.updateheader(nameout, 0, {'zp' + ll: [str(valore), '']})
            if magsat3:
                if readkey3(readhdr(nameout), 'FLUXCAL') == 'ABSOLUTE':
                    try:
                        ntt.util.updateheader(nameout, 0, {
                        'ABMAGSAT': [float(magsat3) + float(readkey3(readhdr(nameout)), 'PHOTZP'),
                                     'Saturation limit for point sources (AB mags)']})
                    except:
                        ntt.util.updateheader(nameout, 0, {
                        'ABMAGSAT': [float(magsat3), 'Saturation limit for point sources (AB mags)']})
                else:
                    ntt.util.updateheader(nameout, 0, {
                    'ABMAGSAT': [float(magsat3), 'Saturation limit for point sources (AB mags)']})
            else:
                ntt.util.updateheader(nameout, 0, {'ABMAGSAT': [9999., 'Saturation limit for point sources (AB mags)']})

            maglim = ntt.util.limmag(nameout)
            if maglim:
                ntt.util.updateheader(nameout, 0,
                                      {'ABMAGLIM': [maglim, '5-sigma limiting AB magnitude for point sources']})
            else:
                ntt.util.updateheader(nameout, 0,
                                      {'ABMAGLIM': [9999., '5-sigma limiting AB magnitude for point sources']})

            if readkey3(readhdr(nameout), 'filter') in ['i705']:
                try:
                    nameout, maskname = ntt.efoscphotredudef.fringing2(nameout, fringingmask, _interactive, False)
                    if nameout not in outputfile:
                        outputfile.append(nameout)
                        if maskname not in outputfile:    outputfile.append(maskname)

                except:
                    ntt.util.writeinthelog('image ' + str(nameout) + ' probably corrupted\n', './logNTT.txt')
                    print '\n### problem with fringing correction'

    print outputfile
    print '\n### adding keywords for phase 3 ........'
    for img in outputfile:
        if str(img)[-5:] == '.fits':
            ################################################
            ntt.util.updateheader(img, 0, {'DETRON ': [11.6, 'Readout noise per output (e-)']})

            try:
                pyv = int(re.sub('\.', '', str(pyfits.__version__))[:2])
            except:
                pyv = 40

            if pyv <= 30:
                print 'here'
                ntt.util.updateheader(img, 0,
                                      {'HIERARCH ESO DET OUT1 GAIN': [1.18, 'Conversion from electrons to ADU']})
                ntt.util.updateheader(img, 0, {'HIERARCH ESO DET OUT1 RON': [11.6, 'Readout noise per output (e-)']})
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
                header['ESO DET OUT1 GAIN'] = (1.18, 'Conversion from electrons to ADU')
                header['ESO DET OUT1 RON'] = (11.6, 'Readout noise per output (e-)')
                imm.flush()
                imm.close()
            #################################################
            hdr = readhdr(img)
            if not readkey3(hdr, 'pixscale'):   ntt.util.updateheader(img, 0, {
            'pixscale': [readkey3(hdr, 'binx') * .12, 'pixel/scale (arcsec)']})

            if 'FLATCOR' in hdr and os.path.isfile(hdr['FLATCOR']):
                hdrn = readhdr(hdr['FLATCOR'])
                if 'NCOMBINE' in hdrn:
                    nflat = hdrn['NCOMBINE']
                else:
                    nflat = 1
            else:
                nflat = 1
            if 'ZEROCOR' in hdr and os.path.isfile(hdr['ZEROCOR']):
                hdrb = readhdr(hdr['ZEROCOR'])
                if 'NCOMBINE' in hdrb:
                    nbias = hdrb['NCOMBINE']
                else:
                    nbias = 1
            else:
                nbias = 1
            ntt.util.updateheader(img, 0, {'EFFRON': [readkey3(hdr, 'ron') * sqrt(1. + 1. / nflat + 1. / nbias),
                                                      'Effective readout noise per output (e-)']})

            try:
                ntt.util.phase3header(img)  #  phase 3 definitions
                hedvec = {'quality': ['Final', 'Final or Fast reduction'],
                          'BUNIT': ['ADU', 'pixel units(ADU,electrons)'],
                          'TEXPTIME': [readkey3(hdr, 'EXPTIME'), 'Total integ. time of all exposure']}
                if readkey3(hdr, 'MJD-OBS'):
                    hedvec['MJD-END'] = [
                        float(readkey3(hdr, 'MJD-OBS')) + (float(readkey3(hdr, 'exptime')) * 0.00001 / (0.864)),
                        'MJD-END = MJD-OBS + JD(exptime)']
                if readkey3(hdr, 'tech'):
                    hedvec['PRODCATG'] = ['SCIENCE.' + readkey3(hdr, 'tech').upper(), 'Data product category']
                ntt.util.updateheader(img, 0, hedvec)

                try:
                    pyv = int(re.sub('\.', '', str(pyfits.__version__))[:2])
                except:
                    pyv = 40

                if pyv <= 30:
                    aa = 'HIERARCH '
                else:
                    aa = ''
                imm = pyfits.open(img, mode='update')
                hdr = imm[0].header
                if aa + 'ESO DPR CATG' in hdr:
                    hdr.pop(aa + 'ESO DPR CATG')
                if aa + 'ESO DPR TECH' in hdr:
                    hdr.pop(aa + 'ESO DPR TECH')
                if aa + 'ESO DPR TYPE' in hdr:
                    hdr.pop(aa + 'ESO DPR TYPE')
                imm.flush()
                imm.close()

            except:
                print 'Warning: ' + img + ' is not a fits file'

    f = open('logfile_phot_' + str(reduceddata) + '_' + str(datenow) + '.raw.list', 'w')
    for img in outputfile:
        try:
            f.write(readkey3(readhdr(img), 'arcfile') + '\n')
        except:
            pass
    f.close()
    return outputfile, 'logfile_phot_' + str(reduceddata) + '_' + str(datenow) + '.raw.list'


# ##############################################################################################
def searchbias(img, listbias):
    from ntt.util import readhdr, readkey3
    import ntt
    import glob, re
    from numpy import argmin
    from numpy import abs

    hdr = readhdr(img)
    JD = readkey3(hdr, 'JD')
    _instrume = readkey3(hdr, 'instrume')
    if not listbias:
        directory = ntt.__path__[0] + '/archive/' + str(_instrume) + '/bias'
        listbias = glob.glob(directory + '/*fits')
    else:
        directory = ''
    if listbias:
        biasfile = ''
        distance = []
        for bias in listbias:
            JDbias = readkey3(readhdr(bias), 'JD')
            distance.append(abs(JD - JDbias))
        if len(distance) >= 1:
            biasfile = listbias[argmin(distance)]
        else:
            biasfile = ''
    else:
        biasfile = ''
    return biasfile, directory


def searchfringe(img, listfringe):
    from ntt.util import readhdr, readkey3
    import ntt
    import glob, re
    from numpy import argmin
    from numpy import abs

    hdr = readhdr(img)
    JD = readkey3(hdr, 'JD')
    _filter = readkey3(hdr, 'filter')
    _instrume = readkey3(hdr, 'instrume')
    if not listfringe:
        directory = ntt.__path__[0] + '/archive/' + str(_instrume) + '/fringing/' + _filter
        listfringe = glob.glob(directory + '/*fits')
    else:
        directory = ''
    if listfringe:
        fringefile = ''
        distance = []
        for fringe in listfringe:
            hdrf = readhdr(fringe)
            JDfringe = readkey3(hdrf, 'JD')
            filterfringe = readkey3(hdrf, 'filter')
            if filterfringe == _filter:
                distance.append(abs(JD - JDfringe))
        if len(distance) >= 1:
            fringefile = listfringe[argmin(distance)]
        else:
            fringefile = ''
    else:
        fringefile = ''
    return fringefile, directory


###################################################################

def rejectflat(lista, _interactive):
    from numpy import where, size
    from ntt.util import display_image
    from pyraf import iraf

    iraf.images(_doprint=0)
    import os, string

    listone = []
    for img in lista:
        dataimg = pyfits.open(img)[0].data
        numberpixels = size(dataimg)
        indices = where(dataimg > 60000)
        numbersaturated = size(dataimg[indices])
        indices2 = where(dataimg < 1000)
        numberlow = size(dataimg[indices2])
        if 100 * float(numbersaturated) / float(numberpixels) <= 5 and \
                                        100 * float(numberlow) / float(numberpixels) <= 10:
            listone.append(img)
    listgood = []
    if _interactive:
        for img in listone:
            aa, bb, cc = display_image(img, 1, '', '', False)
            titolo, result = iraf.imstat(img, Stdout=1)
            print titolo[1:]
            print result
            answ = 'nn'
            while answ not in ['g', 'G', 'b', 's']:
                answ = raw_input('good/bad  [[g]/b/G(all good)/s(stop)]? ')
                if not answ: answ = 'g'
                if answ not in ['g', 'G', 'b', 's']: print 'warning: value not valid, try again'
            if answ == 'g':
                listgood.append(img)
            elif answ == 'G':
                listgood = listone
                break
            elif answ == 's':
                break
    else:
        listgood = listone
    return listgood


#####################################
def rejectbias(lista, _interactive, nn=10):
    import numpy as np
    import ntt
    from pyraf import iraf

    iraf.images(_doprint=0)

    listgood = []
    f = open('_listgoodbias', 'w')
    for img in lista:
        f.write(img + '\n')
    f.close()
    biasstd = np.array(iraf.imstat('@_listgoodbias', fields='stddev', format='no', Stdout=1), float)
    #     biasmedian=np.array(iraf.imstat('@_listgoodbias', fields='mean',format='no',Stdout=1),float)

    #     print lista
    median = np.median(biasstd)
    sigma = (np.percentile(biasstd, 75) - np.percentile(biasstd, 25)) * 1.349
    lista1 = np.compress((biasstd < (median + nn * sigma)) & (biasstd > (median - nn * sigma)), np.array(lista))
    #     lista1=np.compress((np.mean(biasstd)-2*np.std(biasstd)<=biasstd)&(np.mean(biasstd)+2*np.std(biasstd)>=biasstd)&
    #                        (np.mean(biasmedian)-2*np.std(biasmedian)<=biasmedian)&
    #                        (np.mean(biasmedian)+2*np.std(biasmedian)>=biasmedian), np.array(lista))
    ntt.util.delete('_listgoodbias')
    for i in range(0, len(lista)):
        if lista[i] not in lista1:
            print lista[i] + ' ' + str(biasstd[i]) + ' rejected'
        else:
            print lista[i], biasstd[i]
    for img in lista1:
        if _interactive:
            aa, bb, cc = ntt.util.display_image(img, 1, '', '', False)
            titolo, result = iraf.imstat(img, Stdout=1)
            print titolo[1:]
            print result
            answ = 'nn'
            while answ not in ['g', 'G', 'b', 's']:
                answ = raw_input('good/bad  [[g]/b/G(all good)/s(stop) ] ? ')
                if not answ:
                    answ = 'g'
                if answ not in ['g', 'G', 'b', 's']:
                    print 'warning: value not valid, try again'
            if answ == 'g':
                listgood.append(img)
            elif answ == 'G':
                listgood = lista1
                break
            elif answ == 's':
                break
        else:
            listgood.append(img)
    return listgood


################################################################
def fringing2(img, fmask, _interactive, _verbose=False):
    from ntt.util import delete, display_image, updateheader, readhdr, readkey3
    from ntt.efoscphotredudef import searchfringe
    import datetime

    MJDtoday = 55927 + (datetime.date.today() - datetime.date(2012, 01, 01)).days
    import ntt
    import os, string, re
    from numpy import median, where
    from pyraf import iraf

    iraf.nproto(_doprint=0)
    iraf.unlearn('nproto.objmasks')
    if _verbose:
        ver = 'yes'
    else:
        ver = 'no'
    hdr = readhdr(img)
    _filter = readkey3(hdr, 'filter')
    _exptime = readkey3(hdr, 'exptime')
    _date = readkey3(hdr, 'date-night')
    imgout = img
    maskname = ''
    if not readkey3(hdr, 'FRICOR'):
        if fmask:
            fmask = searchfringe(img, fmask)[0]
        if fmask:
            print '###### use fringing mask from user ' + fmask
        else:
            fmask = searchfringe(img, '')[0]
            if fmask:
                print '###### use fringing mask from archive ' + fmask
        if fmask:
            imgout = re.sub('.fits', '_fr.fits', img)
            _trim = readkey3(hdr, 'TRIM')
            _trimmask = readkey3(readhdr(fmask), 'TRIM')
            maskname = 'fmask_' + str(_date) + '_' + _filter + '_' + str(MJDtoday) + '.fits'
            delete(maskname)
            if _trim and not _trimmask:
                _trim = '[' + string.split(_trim, '[')[1]
                iraf.ccdred.ccdproc(fmask, output=maskname, overscan="no", trim="yes", zerocor="no",
                                    trimsec=_trim, flatcor="no", zero="", Stdout=1)
            elif _trim and _trimmask:
                if _trim == _trimmask:
                    os.system('cp ' + str(fmask) + ' ' + str(maskname))
                else:
                    sys.exit('ERROR: fringing correction can be applied only to UNTRIMMED images or '
                             'images with this trim ' + str(_trimmask))
            elif not _trim and _trimmask:
                sys.exit('ERROR: image is not trimmed while selected fringing mask is trimmed: ' + str(_trimmask))
            else:
                os.system('cp ' + str(fmask) + ' ' + str(maskname))

            iraf.nproto.objmasks1.fitxord = 1
            iraf.nproto.objmasks1.fityord = 1
            delete('mask_' + img)
            imgcut = re.sub('.fits', '', img)
            xxx = iraf.nproto.objmasks(images=imgcut, objmasks='mask_' + img, omtype='boolean',
                                       blksize=-16, convolv='block 3 3', hsigma=5, lsigma=3, minpix=10, ngrow=2,
                                       agrow=4., Stdout=1)
            matimg = pyfits.open(img)[0].data
            matmask = pyfits.open('mask_' + img)[1].data
            matfrin = pyfits.open(maskname)[0].data
            indices = where(matmask < 1)
            scalevalue = median((matimg[indices] - median(matimg[indices])) / (matfrin[indices] - median(matfrin)))
            delete(imgout)
            iraf.imutil.imexpr(expr='a - (' + str(scalevalue) + '* (b - ' + str(median(matfrin)) + ') )', a=img,
                               b=maskname, output=imgout, verbose='no')
            print '\n### fringing correction  ..... done'
            iraf.hedit(images=imgout, fields='OBJMASK', value='', delete='yes', update='yes', verify='no')
            ntt.util.updateheader(imgout, 0, {
            'FRICOR': [str(scalevalue) + ' * (' + str(maskname) + ' - ' + str(median(matfrin)) + ')', '']})
            ntt.util.updateheader(imgout, 0, {'FILETYPE': [12205, 'pre-reduced image fringing corrected']})
            ntt.util.updateheader(imgout, 0, {'PROV1': [readkey3(readhdr(imgout), 'ARCFILE'), 'Originating file']})
            ntt.util.updateheader(imgout, 0, {'TRACE1': [img, 'Originating file']})
            stringa = '%7.7s  * (fmask_%8s.fits - %5.5s )' % (str(scalevalue), _date, str(median(matfrin)))
        else:
            print '\n### fringing mask not available for this filter'
    else:
        print '\n### fringing correction already applyed to this image'
    return imgout, maskname

#########################################################
