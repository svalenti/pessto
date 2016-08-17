
def findaperture(img, _interactive=False):
    # print "LOGX:: Entering `findaperture` method/function in %(__file__)s" %
    # globals()
    import re
    import string
    import os
    from pyraf import iraf
    import ntt
    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.specred(_doprint=0)
    toforget = ['specred.apfind']
    for t in toforget:
        iraf.unlearn(t)

    iraf.specred.databas = 'database'
    iraf.specred.dispaxi = 2
    iraf.specred.apedit.thresho = 0

    dv = ntt.dvex()
    grism = ntt.util.readkey3(ntt.util.readhdr(img), 'grism')
    if _interactive:
        _interac = 'yes'
        _edit = 'yes'
    else:
        _interac = 'no'
        _edit = 'no'
    if os.path.isfile('database/ap' + re.sub('.fits', '', img)):
        ntt.util.delete('database/ap' + re.sub('.fits', '', img))
    xx = iraf.specred.apfind(img, interac=_interac, find='yes', recenter='yes', edit=_edit, resize='no',
                             aperture=1, Stdout=1, nfind=1, line=dv['line'][grism], nsum=50, mode='h')
    try:
        for line in open('database/ap' + re.sub('.fits', '', img)):
            if "center" in line:
                center = float(string.split(line)[1])
    except:
        center = 9999
    return center


def calibrationsofi(imgex, stdex, stdref, outputfile, _interactive):
    # print "LOGX:: Entering `calibrationsofi` method/function in
    # %(__file__)s" % globals()
    import pyfits
    import os
    import ntt
    import sys
    import numpy as np

    dataex, hdrex = pyfits.getdata(imgex, 0, header=True)  # read object
    _exptimeimg = ntt.util.readkey3(hdrex, 'exptime')
    _shiftex = ntt.util.readkey3(hdrex, 'shift')
    crvalex = ntt.util.readkey3(hdrex, 'CRVAL1')
    cdex = ntt.util.readkey3(hdrex, 'CD1_1')
    xxex = np.arange(len(dataex[0][0]))
    aaex = crvalex + (xxex) * cdex

    dataref, hdrref = pyfits.getdata(
        stdref, 0, header=True)  # read refrence standard
    crvalref = pyfits.open(stdref)[0].header.get('CRVAL1')
    cdref = pyfits.open(stdref)[0].header.get('CD1_1')
    yyref = pyfits.open(stdref)[0].data
    xxref = np.arange(len(yyref))
    aaref = crvalref + (xxref) * cdref
    yyrefinterp = np.interp(aaex, aaref, yyref)

    datastd, hdrstd = pyfits.getdata(
        stdex, 0, header=True)  # read tellurich standard
    _exptimestd = ntt.util.readkey3(hdrstd, 'exptime')
    _magstd = ntt.util.readkey3(hdrstd, 'magstd')
    _shiftstd = ntt.util.readkey3(hdrstd, 'shift')
    crvalstd = ntt.util.readkey3(hdrstd, 'CRVAL1')
    cdstd = ntt.util.readkey3(hdrstd, 'CD1_1')
    xxstd = np.arange(len(datastd[0][0]))

    if not _shiftex:
        _shiftex = 0  # or _shiftstd!=_shiftex:
    if not _shiftstd:
        _shiftstd = _shiftex  # or _shiftstd!=_shiftex:
    # print '\n### check in wavelengh not performed for tellurich standard\n
    # or std and sn shifts are different\n'
    answ = 'yes'
    while answ == 'yes':
        answ = raw_input('\n### which shift you want to use for the standard ? \n - same as the supernova (' +
                         str(_shiftex) + ')[s] \n' + ' - not shift [n] \n - shift in Ang [' + str(_shiftstd) +
                         ']\n shift [' + str(_shiftstd) + '] ')
        if not answ:
            answ = str(_shiftstd)
    if answ.lower() not in ['s', 'n']:
        try:
            answ = str(float(answ))
        except:
            print '\n### warning: value not valid'
            answ = 'yes'
    if answ.lower() == 'n':
        aastd = crvalstd + (xxstd) * cdstd
    elif answ.lower() == 's':
        aastd = crvalex + (xxstd) * cdstd
    else:
        aastd = crvalstd + float(answ) + (xxstd) * cdstd

    #   scale standard observed to vega or sun magnitude (to be checked)
    #  tellurich stars sun-type are calibrated to the I sun mag
    #  tellurich stars vega-type are calibrated to the V vega mag
    if 'sun' in stdref:
        yystd_sc = (datastd[0][0] / _exptimestd) / \
            ((10 ** ((-27.47) / 2.5)) / (10 ** ((_magstd) / 2.5)))
    elif 'vega' in stdref:
        yystd_sc = (datastd[0][0] / _exptimestd) / \
            ((10 ** ((0.0023) / 2.5)) / (10 ** ((_magstd) / 2.5)))
    else:
        sys.exit('problem with type of star (not sun, not vega) ')
    yystdinterp = np.interp(aaex, aastd, yystd_sc)
    #     calibrate extracted spectrum with tellurich standard
    ntt.util.delete(outputfile)
    dataex[0][0] = (dataex[0][0] / _exptimeimg) / (yystdinterp / yyrefinterp)
    dataex[1][0] = (dataex[1][0] / _exptimeimg) / (yystdinterp / yyrefinterp)
    dataex[2][0] = (dataex[2][0] / _exptimeimg) / (yystdinterp / yyrefinterp)
    dataex[3][0] = (dataex[3][0] / _exptimeimg) / (yystdinterp / yyrefinterp)
    pyfits.writeto(outputfile, np.float32(dataex), hdrex)
    ntt.util.updateheader(outputfile, 0, {'SENSFUN': [stdex, 'tell stand frame'],
                                          'BUNIT': ['erg/cm2/s/Angstrom', 'Physical unit of array values']})
    # get sensitivity teluric function
    senstelluric = 'senstel_' + ntt.util.readkey3(hdrstd, 'date-night') + '_' + \
                   ntt.util.readkey3(hdrstd, 'grism') + '_' + \
        ntt.util.readkey3(hdrstd, 'filter')
    print senstelluric, stdex
    senstelluric = ntt.name_duplicate(stdex, senstelluric, '')
    if os.path.isfile(senstelluric):
        ntt.util.delete(senstelluric)
    yyrefinterp2 = np.interp(aastd, aaref, yyref)
    datasens = 2.5 * np.log10(np.array((yystd_sc / yyrefinterp2) / cdstd))
    # yys=(10**(yys1/2.5))*cds                  # from sens iraf in sens flux
    # yys1=2.5*log10(yys/cds)                   # from sens flux to sens iraf
    pyfits.writeto(senstelluric, np.float32(datasens), hdrstd)
    from pyraf import iraf

    matching = [s for s in hdrex.keys() if "TRACE" in s]
    for imcmb in matching:
        aaa = iraf.hedit(senstelluric, imcmb, delete='yes',
                         update='yes', verify='no', Stdout=1)
    ntt.util.updateheader(
        senstelluric, 0, {'TRACE1': [stdex, 'Originating file']})
    return outputfile, senstelluric


def sofispec1Dredu(files, _interactive, _ext_trace, _dispersionline, _automaticex, _verbose=False):
    # print "LOGX:: Entering `sofispec1Dredu` method/function in %(__file__)s"
    # % globals()
    import re
    import string
    import sys
    import os
    os.environ["PYRAF_BETA_STATUS"] = "1"
    import ntt
    import pyfits
    import numpy as np
    import datetime
    import pylab as pl
    from pyraf import iraf

    dv = ntt.dvex()
    now = datetime.datetime.now()
    datenow = now.strftime('20%y%m%d%H%M')
    MJDtoday = 55927 + (datetime.date.today() - datetime.date(2012, 01, 01)).days
    scal = np.pi / 180.
    hdr0 = ntt.util.readhdr(re.sub('\n', '', files[0]))
    _gain = ntt.util.readkey3(hdr0, 'gain')
    _rdnoise = ntt.util.readkey3(hdr0, 'ron')
    std_sun, rastd_sun, decstd_sun, magstd_sun = ntt.util.readstandard(
        'standard_sofi_sun.txt')
    std_vega, rastd_vega, decstd_vega, magstd_vega = ntt.util.readstandard(
        'standard_sofi_vega.txt')
    std_phot, rastd_phot, decstd_phot, magstd_phot = ntt.util.readstandard(
        'standard_sofi_phot.txt')
    outputfile = []
    objectlist, RA, DEC = {}, {}, {}
    for img in files:
        img = re.sub('\n', '', img)
        hdr = ntt.util.readhdr(img)
        _ra = ntt.util.readkey3(hdr, 'RA')
        _dec = ntt.util.readkey3(hdr, 'DEC')
        _grism = ntt.util.readkey3(hdr, 'grism')
        _filter = ntt.util.readkey3(hdr, 'filter')
        _slit = ntt.util.readkey3(hdr, 'slit')
        cc_sun = np.arccos(np.sin(_dec * scal) * np.sin(decstd_sun * scal) + np.cos(_dec * scal) *
                           np.cos(decstd_sun * scal) * np.cos((_ra - rastd_sun) * scal)) * ((180 / np.pi) * 3600)
        cc_vega = np.arccos(np.sin(_dec * scal) * np.sin(decstd_vega * scal) + np.cos(_dec * scal) *
                            np.cos(decstd_vega * scal) * np.cos((_ra - rastd_vega) * scal)) * ((180 / np.pi) * 3600)
        cc_phot = np.arccos(np.sin(_dec * scal) * np.sin(decstd_phot * scal) + np.cos(_dec * scal) *
                            np.cos(decstd_phot * scal) * np.cos((_ra - rastd_phot) * scal)) * ((180 / np.pi) * 3600)
        if min(cc_sun) < 100:
            _type = 'sun'
        elif min(cc_phot) < 100:
            _type = 'stdp'
        elif min(cc_vega) < 100:
            _type = 'vega'
        else:
            _type = 'obj'
        if min(cc_phot) < 100:
            if _verbose:
                print img, 'phot', str(min(cc_phot)), str(std_phot[np.argmin(cc_phot)])
            ntt.util.updateheader(img, 0, {'stdname': [std_phot[np.argmin(cc_phot)], ''],
                                           'magstd': [float(magstd_phot[np.argmin(cc_phot)]), '']})
        # ntt.util.updateheader(img,0,{'magstd':[float(magstd_phot[argmin(cc_phot)]),'']})
        elif min(cc_sun) < 100:
            if _verbose:
                print img, 'sun', str(min(cc_sun)), str(std_sun[np.argmin(cc_sun)])
            ntt.util.updateheader(img, 0, {'stdname': [std_sun[np.argmin(cc_sun)], ''],
                                           'magstd': [float(magstd_sun[np.argmin(cc_sun)]), '']})
        # ntt.util.updateheader(img,0,{'magstd':[float(magstd_sun[argmin(cc_sun)]),'']})
        elif min(cc_vega) < 100:
            if _verbose:
                print img, 'vega', str(min(cc_vega)), str(std_vega[np.argmin(cc_vega)])
            ntt.util.updateheader(img, 0, {'stdname': [std_vega[np.argmin(cc_vega)], ''],
                                           'magstd': [float(magstd_vega[np.argmin(cc_vega)]), '']})
        # ntt.util.updateheader(img,0,{'magstd':[float(magstd_vega[argmin(cc_vega)]),'']})
        else:
            if _verbose:
                print img, 'object'

        _OBID = (ntt.util.readkey3(hdr, 'esoid'))
        if _type not in objectlist:
            objectlist[_type] = {}
        if _grism not in objectlist[_type]:
            objectlist[_type][_grism] = {}
        if _OBID not in objectlist[_type][_grism]:
            objectlist[_type][_grism][_OBID] = []
        objectlist[_type][_grism][_OBID].append(img)

    if 'stdp' not in objectlist:
        print '###  warning: not photometric standard'
    else:
        print '### photometric standard in the list of object'
    if 'sun' not in objectlist:
        print '### warning: not telluric G standard (sun type)'
    else:
        print '### telluric G standard (sun type) in the list of object'
    if 'vega' not in objectlist:
        print '### warning: not telluric A standard (vega type)'
    else:
        print '### telluric A standard (vega type) in the list of object'

    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.specred(_doprint=0)
    iraf.immatch(_doprint=0)
    iraf.imutil(_doprint=0)
    toforget = ['specred.apall', 'specred.transform']
    for t in toforget:
        iraf.unlearn(t)
    iraf.specred.apall.readnoi = _rdnoise
    iraf.specred.apall.gain = _gain
    iraf.specred.dispaxi = 2
    for _type in objectlist:
        for setup in objectlist[_type]:
            for _ID in objectlist[_type][setup]:
                listmerge = objectlist[_type][setup][_ID]
                listmerge = ntt.sortbyJD(listmerge)
                _object = ntt.util.readkey3(
                    ntt.util.readhdr(listmerge[0]), 'object')
                if string.count(_object, '/') or string.count(_object, '.') or string.count(_object, ' '):
                    nameobj = string.split(_object, '/')[0]
                    nameobj = string.split(nameobj, ' ')[0]
                    nameobj = string.split(nameobj, '.')[0]
                else:
                    nameobj = _object
                _date = ntt.util.readkey3(
                    ntt.util.readhdr(listmerge[0]), 'date-night')
                outputimage = nameobj + '_' + _date + \
                    '_' + setup + '_merge_' + str(MJDtoday)
                outputimage = ntt.util.name_duplicate(
                    listmerge[0], outputimage, '')
                print '### setup= ', setup, ' name field= ', nameobj, ' merge image= ', outputimage, '\n'
                if os.path.isfile(outputimage) and _interactive:
                    answ = raw_input(
                        'combine frame of dithered spectra already created. Do you want to make it again [[y]/n] ? ')
                    if not answ:
                        answ = 'y'
                else:
                    answ = 'y'
                if answ in ['Yes', 'y', 'Y', 'yes']:
                    if _interactive:
                        automaticmerge = raw_input(
                            '\n### Do you want to try to find the dither bethween frames automatically [[y]/n]')
                        if not automaticmerge:
                            automaticmerge = 'yes'
                        elif automaticmerge.lower() in ['y', 'yes']:
                            automaticmerge = 'yes'
                        else:
                            automaticmerge = 'no'
                    else:
                        automaticmerge = 'yes'
                    if automaticmerge == 'yes':
                        offset = 0
                        offsetvec = []
                        _center0 = ntt.sofispec1Ddef.findaperture(
                            listmerge[0], False)
                        _offset0 = ntt.util.readkey3(
                            ntt.util.readhdr(listmerge[0]), 'xcum')
                        print '\n### Try to merge spectra considering their offset along x axes .......'
                        f = open('_offset', 'w')
                        for img in listmerge:
                            _center = ntt.sofispec1Ddef.findaperture(
                                img, False)
                            _center2 = (
                                float(_center) + (float(_offset0) - float(_center0))) * (-1)
                            _offset = (-1) * \
                                ntt.util.readkey3(
                                    ntt.util.readhdr(img), 'xcum')
                            if abs(_center2 - _offset) >= 20:
                                automaticmerge = 'no'
                                break
                            else:
                                offset3 = _center2
                            offsetvec.append(offset3)
                            line = str(offset3) + '   0\n'
                            f.write(line)
                        f.close()
                    if automaticmerge == 'yes':
                        print '### automatic merge .......... done'
                    else:
                        print '\n### warning: try identification of spectra position in interactive way '
                        offset = 0
                        offsetvec = []
                        _z1, _z2, goon = ntt.util.display_image(
                            listmerge[0], 1, '', '', False)
                        print '\n### find aperture on first frame and use it as reference position of ' \
                              'the spectra (mark with ' + '"' + 'm' + '"' + ')'
                        _center0 = ntt.sofispec1Ddef.findaperture(
                            listmerge[0], True)
                        _offset0 = ntt.util.readkey3(
                            ntt.util.readhdr(listmerge[0]), 'xcum')
                        print '\n### find the aperture on all the spectra frames (mark with ' + '"' + 'm' + '"' + ')'
                        f = open('_offset', 'w')
                        for img in listmerge:
                            print '\n### ', img
                            _z1, _z2, goon = ntt.util.display_image(
                                img, 1, '', '', False)
                            _center = ntt.sofispec1Ddef.findaperture(img, True)
                            _center2 = (
                                float(_center) + (float(_offset0) - float(_center0))) * (-1)
                            _offset = (-1) * \
                                ntt.util.readkey3(
                                    ntt.util.readhdr(img), 'xcum')
                            print '\n### position from  dither header: ' + str(_offset)
                            print '### position identified interactively: ' + str(_center2)
                            offset3 = raw_input(
                                '\n### which is the right position [' + str(_center2) + '] ?')
                            if not offset3:
                                offset3 = _center2
                            offsetvec.append(offset3)
                            line = str(offset3) + '   0\n'
                            f.write(line)
                        f.close()
                    start = int(max(offsetvec) - min(offsetvec))
                    f = open('_goodlist', 'w')
                    print listmerge
                    for img in listmerge:
                        f.write(img + '\n')
                    f.close()
                    ntt.util.delete(outputimage)
                    ntt.util.delete('_output.fits')
                    yy1 = pyfits.open(listmerge[0])[0].data[:, 10]
                    iraf.immatch.imcombine('@_goodlist', '_output', combine='sum', reject='none', offset='_offset',
                                           masktyp='', rdnoise=_rdnoise, gain=_gain, zero='mode', Stdout=1)
                    iraf.imutil.imcopy(
                        '_output[' + str(start) + ':1024,*]', output=outputimage, verbose='no')
                    print outputimage
                    print len(listmerge)
                    hdr1 = ntt.util.readhdr(outputimage)
                    ntt.util.updateheader(outputimage, 0,
                                          {'SINGLEXP': [False, 'TRUE if resulting from single exposure'],
                                           'M_EPOCH': [False, 'TRUE if resulting from multiple epochs'],
                                           'EXPTIME': [ntt.util.readkey3(hdr1, 'EXPTIME') * len(listmerge),
                                                       'Total integration time per pixel (s)'],
                                           'TEXPTIME': [float(ntt.util.readkey3(hdr1, 'TEXPTIME')) * len(listmerge),
                                                        'Total integration time of all exposures (s)'],
                                           'APERTURE': [2.778e-4 * float(re.sub('long_slit_', '',
                                                                                ntt.util.readkey3(hdr1, 'slit'))),
                                                        '[deg] Aperture diameter'],
                                           'NOFFSETS': [2, 'Number of offset positions'],
                                           'NUSTEP': [0, 'Number of microstep positions'],
                                           'NJITTER': [int(ntt.util.readkey3(hdr1, 'NCOMBINE') / 2),
                                                       'Number of jitter positions']})
                    hdr = ntt.util.readhdr(outputimage)
                    matching = [s for s in hdr.keys() if "IMCMB" in s]
                    for imcmb in matching:
                        aaa = iraf.hedit(outputimage, imcmb, delete='yes', update='yes',
                                         verify='no', Stdout=1)
                    if 'SKYSUB' in hdr.keys():
                        aaa = iraf.hedit(outputimage, 'SKYSUB', delete='yes', update='yes',
                                         verify='no', Stdout=1)

                    mjdend = []
                    mjdstart = []
                    num = 0
                    for img in listmerge:
                        num = num + 1
                        hdrm = ntt.util.readhdr(img)
                        ntt.util.updateheader(outputimage, 0,
                                              {'PROV' + str(num):
                                               [ntt.util.readkey3(
                                                   hdrm, 'ARCFILE'), 'Originating file'],
                                               'TRACE' + str(num): [img, 'Originating file']})
                        mjdend.append(ntt.util.readkey3(hdrm, 'MJD-END'))
                        mjdstart.append(ntt.util.readkey3(hdrm, 'MJD-OBS'))
                    _dateobs = ntt.util.readkey3(ntt.util.readhdr(
                        listmerge[np.argmin(mjdstart)]), 'DATE-OBS')

                    _telapse = (max(mjdend) - min(mjdstart)) * \
                        60. * 60 * 24.  # *86400
                    _tmid = (max(mjdend) + min(mjdstart)) / 2

                    _title = str(_tmid)[0:9] + ' ' + str(ntt.util.readkey3(hdr, 'object')) + ' ' + str(
                        ntt.util.readkey3(hdr, 'grism')) + ' ' + \
                        str(ntt.util.readkey3(hdr, 'filter')) + \
                        ' ' + str(ntt.util.readkey3(hdr, 'slit'))
                    ntt.util.updateheader(outputimage, 0,
                                          {'MJD-OBS': [min(mjdstart), 'MJD start'],
                                           'MJD-END': [max(mjdend), 'MJD end'],
                                           'TELAPSE': [_telapse, 'Total elapsed time [days]'],
                                           'TMID': [_tmid, '[d] MJD mid exposure'],
                                           'TITLE': [_title, 'Dataset title'],
                                           'DATE-OBS': [_dateobs, 'Date of observation']})
                    # missing: merge airmass
                else:
                    print '\n### skip making again combined spectrum'
                objectlist[_type][setup][_ID] = [outputimage]
                print '\n### setup= ', setup, ' name field= ', nameobj, ' merge image= ', outputimage, '\n'
                if outputimage not in outputfile:
                    outputfile.append(outputimage)
                ntt.util.updateheader(outputimage, 0, {'FILETYPE': [
                                      42116, 'combine 2D spectra frame']})

    if _verbose:
        if 'obj' in objectlist:
            print objectlist['obj']
        if 'stdp' in objectlist:
            print objectlist['stdp']
        if 'sun' in objectlist:
            print objectlist['sun']
        if 'vega' in objectlist:
            print objectlist['vega']

    if 'obj' not in objectlist.keys():
        sys.exit('\n### error: no objects in the list')

    sens = {}
    print '\n############################################\n### extract the spectra  '
    # print objectlist
    for setup in objectlist['obj']:
        reduced = []
        for _ID in objectlist['obj'][setup]:
            for img in objectlist['obj'][setup][_ID]:
                hdr = ntt.util.readhdr(img)
                print '\n### next object\n ', img, ntt.util.readkey3(hdr, 'object')
                _grism = ntt.util.readkey3(hdr, 'grism')
                _exptimeimg = ntt.util.readkey3(hdr, 'exptime')
                _JDimg = ntt.util.readkey3(hdr, 'JD')

                imgex = ntt.util.extractspectrum(img, dv, _ext_trace, _dispersionline, _interactive, 'obj',
                                                 automaticex=_automaticex)
                if imgex not in outputfile:
                    outputfile.append(imgex)
                ntt.util.updateheader(imgex, 0, {'FILETYPE': [42107, 'extracted 1D wave calib'],
                                                 'PRODCATG': ['SCIENCE.' + ntt.util.readkey3(hdr, 'tech').upper(),
                                                              'Data product category']})
                hdr = ntt.util.readhdr(imgex)
                matching = [s for s in hdr.keys() if "TRACE" in s]
                for imcmb in matching:
                    aaa = iraf.hedit(imgex, imcmb, delete='yes',
                                     update='yes', verify='no', Stdout=1)
                ntt.util.updateheader(
                    imgex, 0, {'TRACE1': [img, 'Originating file']})

                if os.path.isfile('database/ap' + re.sub('_ex.fits', '', imgex)):
                    if 'database/ap' + re.sub('_ex.fits', '', imgex) not in outputfile:
                        outputfile.append(
                            'database/ap' + re.sub('_ex.fits', '', imgex))

                ###########################   telluric standard   #############
                if 'sun' in objectlist and setup in objectlist['sun']:
                    _type = 'sun'
                elif 'vega' in objectlist and setup in objectlist['vega']:
                    _type = 'vega'
                else:
                    _type = 'none'
                if _type in ['sun', 'vega']:
                    stdref = ntt.__path__[
                        0] + '/standard/fits/' + str(_type) + '.fits'
                    stdvec, airmassvec, JDvec = [], [], []
                    for _ID in objectlist[_type][setup]:
                        for std in objectlist[_type][setup][_ID]:
                            _airmassstd = ntt.util.readkey3(
                                ntt.util.readhdr(std), 'airmass')
                            _JDstd = ntt.util.readkey3(
                                ntt.util.readhdr(std), 'JD')
                            JDvec.append(abs(_JDstd - _JDimg))
                            stdvec.append(std)
                            airmassvec.append(_airmassstd)
                    stdtelluric = stdvec[np.argmin(JDvec)]
                    _exptimestd = ntt.util.readkey3(
                        ntt.util.readhdr(stdtelluric), 'exptime')
                    _magstd = ntt.util.readkey3(
                        ntt.util.readhdr(stdtelluric), 'magstd')
                    print '\n\n ##### closer standard for telluric corrections  #### \n\n'
                    print stdtelluric, airmassvec[np.argmin(JDvec)]
                    stdtelluric_ex = ntt.util.extractspectrum(stdtelluric, dv, False, False, _interactive, 'std',
                                                              automaticex=_automaticex)
                    if stdtelluric_ex not in outputfile:
                        outputfile.append(stdtelluric_ex)
                    ntt.util.updateheader(stdtelluric_ex, 0, {'FILETYPE': [
                                          42107, 'extracted 1D wave calib ']})
                    ntt.util.updateheader(stdtelluric_ex, 0, {'PRODCATG': [
                        'SCIENCE.' + ntt.util.readkey3(
                            ntt.util.readhdr(stdtelluric_ex), 'tech').upper(), 'Data product category']})

                    hdr = ntt.util.readhdr(stdtelluric_ex)
                    matching = [s for s in hdr.keys() if "TRACE" in s]
                    for imcmb in matching:
                        aaa = iraf.hedit(
                            stdtelluric_ex, imcmb, delete='yes', update='yes', verify='no', Stdout=1)
                    ntt.util.updateheader(stdtelluric_ex, 0, {'TRACE1': [
                                          stdtelluric, 'Originating file']})
                    ###########################################################
                    #               SN tellurich calibration
                    imgf = re.sub('_ex.fits', '_f.fits', imgex)
                    imgf, senstelluric = ntt.sofispec1Ddef.calibrationsofi(imgex, stdtelluric_ex, stdref, imgf,
                                                                           _interactive)
                    if imgf not in outputfile:
                        outputfile.append(imgf)
                    if senstelluric not in outputfile:
                        outputfile.append(senstelluric)
                    ntt.util.updateheader(imgf, 0, {'FILETYPE': [42208, '1D wave calib, tell cor.'],
                                                    #                                                    'SNR': [ntt.util.StoN(imgf, 50),
                                                    'SNR': [ntt.util.StoN2(imgf, False),
                                                            'Average signal to noise ratio per pixel'],
                                                    'TRACE1': [imgex, 'Originating file'],
                                                    'ASSON1': [re.sub('_f.fits', '_2df.fits', imgf),
                                                               'Name of associated file'],
                                                    'ASSOC1': ['ANCILLARY.2DSPECTRUM', 'Category of associated file']})
                    ###########################################################
                    imgd = ntt.efoscspec1Ddef.fluxcalib2d(
                        img, senstelluric)  # flux calibration 2d images
                    ntt.util.updateheader(
                        imgd, 0, {'FILETYPE': [42209, '2D wavelength and flux calibrated spectrum']})
                    iraf.hedit(imgd, 'PRODCATG', delete='yes',
                               update='yes', verify='no')
                    hdrd = ntt.util.readhdr(imgd)
                    matching = [s for s in hdrd.keys() if "TRACE" in s]
                    for imcmb in matching:
                        aaa = iraf.hedit(
                            imgd, imcmb, delete='yes', update='yes', verify='no', Stdout=1)
                    ntt.util.updateheader(
                        imgd, 0, {'TRACE1': [img, 'Originating file']})
                    if imgd not in outputfile:
                        outputfile.append(imgd)
                ###############################################################
                if 'stdp' in objectlist and setup in objectlist['stdp']:
                    print '\n #####  photometric calibration   ######\n '
                    standardfile = []
                    for _ID in objectlist['stdp'][setup]:
                        for stdp in objectlist['stdp'][setup][_ID]:
                            stdp_ex = ntt.util.extractspectrum(stdp, dv, False, _dispersionline, _interactive, 'std',
                                                               automaticex=_automaticex)
                            standardfile.append(stdp_ex)
                            if stdp_ex not in outputfile:
                                outputfile.append(stdp_ex)
                            ntt.util.updateheader(stdp_ex, 0, {
                                'FILETYPE': [42107, 'extracted 1D wave calib'],
                                'TRACE1': [stdp_ex, 'Originating file'],
                                'PRODCATG': ['SCIENCE.' + ntt.util.readkey3(ntt.util.readhdr(stdp_ex), 'tech').upper(),
                                             'Data product category']})
                    print '\n### ', standardfile, ' \n'
                    if len(standardfile) >= 2:
                        standardfile0 = raw_input(
                            'which one do you want to use [' + str(standardfile[0]) + '] ? ')
                        if not standardfile0:
                            standardfile0 = standardfile[0]
                    else:
                        standardfile0 = standardfile[0]
                    print standardfile0
                    stdpf = re.sub('_ex.fits', '_f.fits', standardfile0)
                    stdpf, senstelluric2 = ntt.sofispec1Ddef.calibrationsofi(standardfile0, stdtelluric_ex, stdref,
                                                                             stdpf, _interactive)
                    if stdpf not in outputfile:
                        outputfile.append(stdpf)
                    ntt.util.updateheader(stdpf, 0, {'FILETYPE': [42208, '1D wave calib, tell cor'],
                                                     'TRACE1': [stdp, 'Originating file']})
                    stdname = ntt.util.readkey3(
                        ntt.util.readhdr(standardfile0), 'stdname')
                    standardfile = ntt.__path__[
                        0] + '/standard/flux/' + stdname
                    xx, yy = ntt.util.ReadAscii2(standardfile)

                    crval1 = pyfits.open(stdpf)[0].header.get('CRVAL1')
                    cd1 = pyfits.open(stdpf)[0].header.get('CD1_1')
                    datastdpf, hdrstdpf = pyfits.getdata(stdpf, 0, header=True)
                    xx1 = np.arange(len(datastdpf[0][0]))
                    aa1 = crval1 + (xx1) * cd1
                    yystd = np.interp(aa1, xx, yy)
                    rcut = np.compress(
                        ((aa1 < 13000) | (aa1 > 15150)) & ((11700 < aa1) | (aa1 < 11000)) & (aa1 > 10000) &
                        ((aa1 < 17800) | (aa1 > 19600)) & (aa1 < 24000), datastdpf[0][0] / yystd)
                    aa11 = np.compress(
                        ((aa1 < 13000) | (aa1 > 15150)) & ((11700 < aa1) | (aa1 < 11000)) & (aa1 > 10000) &
                        ((aa1 < 17800) | (aa1 > 19600)) & (aa1 < 24000), aa1)
                    yy1clean = np.interp(aa1, aa11, rcut)
                    aa1 = np.array(aa1)
                    yy1clean = np.array(yy1clean)
                    A = np.ones((len(rcut), 2), dtype=float)
                    A[:, 0] = aa11
                    result = np.linalg.lstsq(A, rcut)  # result=[zero,slope]
                    p = [result[0][1], result[0][0]]
                    yfit = ntt.util.pval(aa1, p)

                    pl.clf()
                    pl.ion()
                    pl.plot(aa1, datastdpf[0][0] / yystd,
                            color='red', label='std')
                    pl.plot(aa1, yfit, color='blue', label='fit')
                    pl.legend(numpoints=1, markerscale=1.5)
                    #   sens function sofi spectra
                    outputsens = 'sens_' + stdpf
                    ntt.util.delete(outputsens)
                    datastdpf[0][0] = yfit
                    pyfits.writeto(outputsens, np.float32(datastdpf), hdrstdpf)
                    #################
                    imgsc = re.sub('_ex.fits', '_sc.fits', imgex)
                    ntt.util.delete(imgsc)
                    crval2 = pyfits.open(imgf)[0].header.get('CRVAL1')
                    cd2 = pyfits.open(imgf)[0].header.get('CD1_1')
                    dataf, hdrf = pyfits.getdata(imgf, 0, header=True)
                    xx2 = np.arange(len(dataf[0][0]))
                    aa2 = crval2 + (xx2) * cd2
                    yyscale = np.interp(aa2, aa1, yfit)

                    dataf[0][0] = dataf[0][0] / yyscale
                    dataf[1][0] = dataf[1][0] / yyscale
                    dataf[2][0] = dataf[2][0] / yyscale
                    dataf[3][0] = dataf[3][0] / yyscale

                    pyfits.writeto(imgsc, np.float32(dataf), hdrf)
                    ntt.util.updateheader(imgsc, 0, {'SENSPHOT': [outputsens, 'sens used to flux cal'],
                                                     'FILETYPE': [42208, '1D wave,flux calib, tell cor'],
                                                     'TRACE1': [imgf, 'Originating file']})
                    #                ntt.util.updateheader(imgsc,0,{'FILETYPE':[42208,'1D wave,flux calib, tell cor']})
                    #                ntt.util.updateheader(imgsc,0,{'TRACE1':[imgf,'']})
                    print '\n### flux calibrated spectrum= ', imgf, ' with the standard= ', stdpf
                    if imgsc not in outputfile:
                        outputfile.append(imgsc)
                else:
                    print '\n### photometric calibrated not performed \n'

    print '\n### adding keywords for phase 3 ....... '
    reduceddata = ntt.util.rangedata(outputfile)
    f = open('logfile_spec1d_' + str(reduceddata) +
             '_' + str(datenow) + '.raw.list', 'w')
    for img in outputfile:
        if str(img)[-5:] == '.fits':
            hdr = ntt.util.readhdr(img)
            # added for DR2
            if 'NCOMBINE' in hdr:
                _ncomb = ntt.util.readkey3(hdr, 'NCOMBINE')
            else:
                _ncomb = 1.0

            _effron = 12. * \
                (1 / np.sqrt(ntt.util.readkey3(hdr, 'ndit') * _ncomb)) * \
                np.sqrt(np.pi / 2)

            try:
                ntt.util.phase3header(img)  # phase 3 definitions
                ntt.util.updateheader(img, 0, {'quality': ['Final', ''],
                                               'EFFRON': [_effron, 'Effective readout noise per output (e-)']})
                f.write(ntt.util.readkey3(
                    ntt.util.readhdr(img), 'arcfile') + '\n')
            except:
                print 'Warning: ' + img + ' is not a fits file'
    f.close()
    return outputfile, 'logfile_spec1d_' + str(reduceddata) + '_' + str(datenow) + '.raw.list'
# #############################################################################################
