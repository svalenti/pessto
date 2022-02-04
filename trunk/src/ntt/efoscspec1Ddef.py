def telluric_atmo(imgstd):

    import numpy as np
    import ntt
    from pyraf import iraf

    try:        import pyfits
    except:     from astropy.io import fits as pyfits

    iraf.images(_doprint=0, Stdout=0)
    iraf.noao(_doprint=0, Stdout=0)
    iraf.twodspec(_doprint=0, Stdout=0)
    iraf.longslit(_doprint=0, Stdout=0)
    iraf.onedspec(_doprint=0, Stdout=0)
    toforget = ['imfilter.gauss', 'specred.apall', 'longslit.identify', 'longslit.reidentify', 'specred.standard',
                'onedspec.wspectext']
    for t in toforget:
        iraf.unlearn(t)

    _grism = ntt.util.readkey3(ntt.util.readhdr(imgstd), 'grism')
    imgout = 'invers_atmo_' + imgstd
    ntt.util.delete(imgout)
    iraf.set(direc=ntt.__path__[0] + '/')
    _cursor = 'direc$standard/ident/cursor_sky_0'
    iraf.noao.onedspec.bplot(imgstd, cursor=_cursor,
                             spec2=imgstd, new_ima=imgout, overwri='yes')
    xxstd, ffstd = ntt.util.readspectrum(imgout)
    if _grism in ['Gr13', 'Gr16']:
        llo2 = np.compress((np.array(xxstd) >= 7550) & (
            np.array(xxstd) <= 7750), np.array(xxstd))
        llh2o = np.compress((np.array(xxstd) >= 7100) & (
            np.array(xxstd) <= 7500), np.array(xxstd))
        ffo2 = np.compress((np.array(xxstd) >= 7550) & (
            np.array(xxstd) <= 7750), np.array(ffstd))
        ffh2o = np.compress((np.array(xxstd) >= 7100) & (
            np.array(xxstd) <= 7500), np.array(ffstd))
    elif _grism in ['Gr11', 'Gr18', 'Gr20']:
        llo2 = np.compress((np.array(xxstd) >= 6830) & (
            np.array(xxstd) <= 7100), np.array(xxstd))
        llh2o = np.compress((np.array(xxstd) >= 7100) & (
            np.array(xxstd) <= 7500), np.array(xxstd))
        ffo2 = np.compress((np.array(xxstd) >= 6830) & (
            np.array(xxstd) <= 7100), np.array(ffstd))
        ffh2o = np.compress((np.array(xxstd) >= 7100) & (
            np.array(xxstd) <= 7500), np.array(ffstd))
    if _grism in ['Gr13', 'Gr16', 'Gr11', 'Gr18', 'Gr20']:
        _skyfileh2o = 'direc$standard/ident/ATLAS_H2O.fits'
        _skyfileo2 = 'direc$standard/ident/ATLAS_O2.fits'
        atlas_smooto2 = '_atlas_smoot_o2.fits'
        atlas_smooth2o = '_atlas_smoot_h2o.fits'
        _sigma = 200
        ntt.util.delete(atlas_smooto2)
        ntt.util.delete(atlas_smooth2o)
        iraf.imfilter.gauss(_skyfileh2o, output=atlas_smooth2o, sigma=_sigma)
        iraf.imfilter.gauss(_skyfileo2, output=atlas_smooto2, sigma=_sigma)
        llskyh2o, ffskyh2o = ntt.util.readspectrum(atlas_smooth2o)
        llskyo2, ffskyo2 = ntt.util.readspectrum(atlas_smooto2)
        ffskyo2cut = np.interp(llo2, llskyo2, ffskyo2)
        ffskyh2ocut = np.interp(llh2o, llskyh2o, ffskyh2o)
        _scaleh2o = []
        integral_h2o = []
        for i in range(1, 21):
            j = 0.6 + i * 0.04
            _ffskyh2ocut = list((np.array(ffskyh2ocut) * j) + 1 - j)
            diff_h2o = abs(_ffskyh2ocut - ffh2o)
            integraleh2o = np.trapz(diff_h2o, llh2o)
            integral_h2o.append(integraleh2o)
            _scaleh2o.append(j)
        _scaleo2 = []
        integral_o2 = []
        for i in range(1, 21):
            j = 0.6 + i * 0.04
            _ffskyo2cut = list((np.array(ffskyo2cut) * j) + 1 - j)
            diff_o2 = abs(_ffskyo2cut - ffo2)
            integraleo2 = np.trapz(diff_o2, llo2)
            integral_o2.append(integraleo2)
            _scaleo2.append(j)
        sh2o = _scaleh2o[np.argmin(integral_h2o)]
        so2 = _scaleo2[np.argmin(integral_o2)]
        telluric_features = ((np.array(ffskyh2o) * sh2o) +
                             1 - sh2o) + ((np.array(ffskyo2) * so2) + 1 - so2) - 1
        telluric_features = np.array([1] + list(telluric_features) + [1])
        llskyo2 = np.array([1000] + list(llskyo2) + [15000])
        telluric_features_cut = np.interp(xxstd, llskyo2, telluric_features)
        _imgout = 'atmo_' + imgstd

        data1, hdr = pyfits.getdata(imgstd, 0, header=True)
        data1[0] = np.array(telluric_features_cut)
        data1[1] = data1[1] / data1[1]
        data1[2] = data1[2] / data1[2]
        data1[3] = data1[3] / data1[3]
        ntt.util.delete(_imgout)
        pyfits.writeto(_imgout, np.float32(data1), hdr)
        ntt.util.delete(atlas_smooto2)
        ntt.util.delete(atlas_smooth2o)
        ntt.util.delete(imgout)
    else:
        _imgout = ''
        print('### telluric correction with model not possible ')
    return _imgout


def fluxcalib2d(img2d, sensfun):  # flux calibrate 2d images

    try:        import pyfits
    except:     from astropy.io import fits as pyfits

    import re
    import string
    import numpy as np
    import ntt

    data2d, hdr2d = pyfits.getdata(img2d, 0, header=True)
    xxd = np.arange(len(data2d[:, 0]))
    # aad=crvald+(xxd)*cdd
    # crvald=readkey3(hdr2d,'CRVAL2')
    # cdd=readkey3(hdr2d,'CD2_2')
    crvald = pyfits.open(img2d)[0].header.get('CRVAL2')
    cdd = pyfits.open(img2d)[0].header.get('CD2_2')
    _exptime = ntt.util.readkey3(ntt.util.readhdr(img2d), 'exptime')
    _airmass = ntt.util.readkey3(ntt.util.readhdr(img2d), 'airmass')
    #  read sensfunction and interpole pixel of 2D image
    yys = pyfits.open(sensfun)[0].data
    crvals = pyfits.open(sensfun)[0].header.get('CRVAL1')
    cds = pyfits.open(sensfun)[0].header.get('CD1_1')
    yys = (10 ** (yys / 2.5)) * cds  # from sens iraf in sens flux
    xxs = np.arange(len(yys))
    aasens = crvals + (xxs) * cds
    xxs2 = (aasens - crvald) / cdd
    aasens2 = np.interp(xxd, xxs2, yys)
    #  read atmosferic function and interpole pixel of 2D image
    aae, yye = ntt.util.ReadAscii2(
        ntt.__path__[0] + '/standard/extinction/lasilla2.txt')
    aae, yye = np.array(aae, float), np.array(yye, float)
    xxe = (aae - crvald) / cdd
    atm_xx = np.interp(xxd, xxe, yye)
    aircorr = 10 ** (0.4 * np.array(atm_xx) * _airmass)
    img2df = re.sub('.fits', '_2df.fits', img2d)
    for i in range(len(data2d[0, :])):
        data2d[:, i] = ((np.array(data2d[:, i] / _exptime) *
                         np.array(aircorr)) / aasens2) * 1e20
    ntt.util.delete(img2df)
    pyfits.writeto(img2df, np.float32(data2d), hdr2d)
    ntt.util.updateheader(
        img2df, 0, {'SENSFUN': [sensfun.split('/')[-1], '']})
    ntt.util.updateheader(img2df, 0, {
                          'BUNIT': ['10^20 erg/cm2/s/Angstrom', 'Physical unit of array values']})
    return img2df


def checkwavestd(imgex, _interactive):

    import ntt
    import numpy as np

    try:        import pyfits
    except:     from astropy.io import fits as pyfits

    print('\n### Warning: check in wavelenght with sky lines not performed\n')
    if _interactive in ['yes', 'YES', 'Yes', 'Y', 'y']:
        answ = input(
            '\n### Do you want to check the wavelengh calibration with tellurich lines [[y]/n]? ')
        if not answ:
            answ = 'y'
    else:
        answ = 'y'
    if answ in ['y', 'yes']:
        print('\n### check wavelength calibration with tellurich lines \n')
	# sky
        _skyfile = ntt.__path__[0] + '/standard/ident/sky_new_0.fits'
        skyff = 1 - (pyfits.open(_skyfile)[0].data)
        crval1 = pyfits.open(_skyfile)[0].header.get('CRVAL1')
        cd1 = pyfits.open(_skyfile)[0].header.get('CD1_1')
        skyxx = np.arange(len(skyff))
        skyaa = crval1 + (skyxx) * cd1

	# object
        atmofile = ntt.efoscspec1Ddef.atmofile(imgex, 'atmo2_' + imgex)
        atmoff = 1 - (pyfits.open(atmofile)[0].data[0][0])
        crval1 = pyfits.open(atmofile)[0].header.get('CRVAL1')
        cd1 = pyfits.open(atmofile)[0].header.get('CD1_1')
        atmoxx = np.arange(len(atmoff))
        atmoaa = crval1 + (atmoxx) * cd1

        if 'Gr18' in imgex.split('_'):
	        shift = ntt.efoscspec2Ddef.checkwavelength_arc(
                atmoaa, atmoff, skyaa, skyff, 5500, 6800)  # for plotting purposes only.
        else:
	        shift = ntt.efoscspec2Ddef.checkwavelength_arc(
                atmoaa, atmoff, skyaa, skyff, 6800, 7800)
    else:
        shift = 0
    zro = pyfits.open(imgex)[0].header.get('CRVAL1')
    if _interactive in ['yes', 'YES', 'Yes', 'Y', 'y']:
        answ = input(
            '\n### do you want to correct the wavelengh calibration with this shift: ' + str(shift) + ' [[y]/n] ? ')
        if not answ:
            answ = 'y'
        if answ.lower() in ['y', 'yes']:
            ntt.util.updateheader(imgex, 0, {'CRVAL1': [zro + int(shift), '']})
            ntt.util.updateheader(imgex, 0, {'shift': [float(shift), '']})
    else:
        ntt.util.updateheader(imgex, 0, {'CRVAL1': [zro + int(shift), '']})
        ntt.util.updateheader(imgex, 0, {'shift': [float(shift), '']})


# ###################################

def atmofile(imgstd, imgout=''):

    from pyraf import iraf
    import os
    import ntt

    iraf.noao(_doprint=0, Stdout=0)
    iraf.onedspec(_doprint=0, Stdout=0)
    iraf.set(direc=ntt.__path__[0] + '/')
    _cursor = 'direc$standard/ident/cursor_sky_0'
    if not imgout:
        imgout = 'atmo_' + imgstd
    os.system('rm -rf ' + imgout)
    iraf.noao.onedspec.bplot(imgstd, cursor=_cursor,
                             spec2=imgstd, new_ima=imgout, overwri='yes')
    return imgout


def sensfunction(standardfile, _function, _order, _interactive):

    import re
    import os
    import sys
    import ntt
    import datetime

    try:       import pyfits  #  added later
    except:    from astropy.io import fits as pyfits

    from pyraf import iraf
    import numpy as np

    MJDtoday = 55927 + (datetime.date.today() - datetime.date(2012, 1, 1)).days
    iraf.noao(_doprint=0, Stdout=0)
    iraf.imred(_doprint=0, Stdout=0)
    iraf.specred(_doprint=0, Stdout=0)
    toforget = ['specred.scopy', 'specred.sensfunc', 'specred.standard']
    for t in toforget:
        iraf.unlearn(t)
    iraf.specred.scopy.format = 'multispec'
    iraf.specred.verbose = 'no'
    hdrs = ntt.util.readhdr(standardfile)
    try:
        _outputsens = 'sens_' + str(ntt.util.readkey3(hdrs, 'date-night')) + '_' + \
                      str(ntt.util.readkey3(hdrs, 'grism')) + '_' + str(ntt.util.readkey3(hdrs, 'filter')) + '_' + \
                      re.sub('.dat', '', ntt.util.readkey3(
                          hdrs, 'stdname')) + '_' + str(MJDtoday)
    except:
        sys.exit('Error: missing header -stdname- in standard ' +
                 str(standardfile) + '  ')

    _outputsens = ntt.util.name_duplicate(standardfile, _outputsens, '')
    if os.path.isfile(_outputsens):
        if _interactive.lower() != 'yes':
            ntt.util.delete(_outputsens)
        else:
            answ = input(
                'sensitivity function already computed, do you want to do it again [[y]/n] ? ')
            if not answ:
                answ = 'y'
            if answ.lower() in ['y', 'yes']:
                ntt.util.delete(_outputsens)

    if not os.path.isfile(_outputsens):
        iraf.set(direc=ntt.__path__[0] + '/')
        _caldir = 'direc$standard/MAB/'
        _extinctdir = 'direc$standard/extinction/'
        _observatory = 'lasilla'
        _extinction = 'lasilla2.txt'
        refstar = 'm' + \
            re.sub('.dat', '', pyfits.open(standardfile)
                   [0].header.get('stdname'))
        _airmass = ntt.util.readkey3(hdrs, 'airmass')
        _exptime = ntt.util.readkey3(hdrs, 'exptime')
        _outputstd = 'std_' + str(ntt.util.readkey3(hdrs, 'grism')) + '_' + \
                     str(ntt.util.readkey3(hdrs, 'filter')) + '.fits'
        ntt.util.delete(_outputstd)
        ntt.util.delete(_outputsens)
        iraf.specred.standard(input=standardfile, output=_outputstd, extinct=_extinctdir + _extinction,
                              caldir=_caldir, observa=_observatory, star_nam=refstar, airmass=_airmass,
                              exptime=_exptime, interac=_interactive)
        iraf.specred.sensfunc(standard=_outputstd, sensitiv=_outputsens, extinct=_extinctdir + _extinction,
                              ignorea='yes', observa=_observatory, functio=_function, order=_order,
                              interac=_interactive)

        data, hdr = pyfits.getdata(standardfile, 0, header=True)  # added later
        data1, hdr1 = pyfits.getdata(
            _outputsens, 0, header=True)  # added later
        ntt.util.delete(_outputsens)  # added later
        pyfits.writeto(_outputsens, np.float32(data1), hdr)  # added later
    return _outputsens


def efoscspec1Dredu(files, _interactive, _ext_trace, _dispersionline, liststandard, listatmo0, _automaticex,
                    _verbose=False):

    import ntt

    try:        import pyfits
    except:     from astropy.io import fits as pyfits

    import re
    import string
    import sys
    import os
    import numpy as np

    os.environ["PYRAF_BETA_STATUS"] = "1"
    _extinctdir = 'direc$standard/extinction/'
    _extinction = 'lasilla2.txt'
    _observatory = 'lasilla'
    import datetime

    now = datetime.datetime.now()
    datenow = now.strftime('20%y%m%d%H%M')
    MJDtoday = 55927 + (datetime.date.today() - datetime.date(2012, 1, 1)).days
    dv = ntt.dvex()
    scal = np.pi / 180.
    _gain = ntt.util.readkey3(ntt.util.readhdr(
        re.sub('\n', '', files[0])), 'gain')
    _rdnoise = ntt.util.readkey3(
        ntt.util.readhdr(re.sub('\n', '', files[0])), 'ron')
    std, rastd, decstd, magstd = ntt.util.readstandard(
        'standard_efosc_mab.txt')
    objectlist = {}
    for img in files:
        hdr = ntt.util.readhdr(img)
        img = re.sub('\n', '', img)
        ntt.util.correctcard(img)
        _ra = ntt.util.readkey3(hdr, 'RA')
        _dec = ntt.util.readkey3(hdr, 'DEC')
        _object = ntt.util.readkey3(hdr, 'object')
        _grism = ntt.util.readkey3(hdr, 'grism')
        _filter = ntt.util.readkey3(hdr, 'filter')
        _slit = ntt.util.readkey3(hdr, 'slit')
        dd = np.arccos(np.sin(_dec * scal) * np.sin(decstd * scal) + np.cos(_dec * scal) *
                       np.cos(decstd * scal) * np.cos((_ra - rastd) * scal)) * ((180 / np.pi) * 3600)
        if min(dd) < 100:
            _type = 'stdsens'
        else:
            _type = 'obj'
        if min(dd) < 100:
            ntt.util.updateheader(
                img, 0, {'stdname': [std[np.argmin(dd)], '']})
            ntt.util.updateheader(
                img, 0, {'magstd': [float(magstd[np.argmin(dd)]), '']})

        if _type not in objectlist:
            objectlist[_type] = {}
        if (_grism, _filter, _slit) not in objectlist[_type]:
            objectlist[_type][_grism, _filter, _slit] = [img]
        else:
            objectlist[_type][_grism, _filter, _slit].append(img)

    from pyraf import iraf
    iraf.noao(_doprint=0, Stdout=0)
    iraf.imred(_doprint=0, Stdout=0)
    iraf.specred(_doprint=0, Stdout=0)
    iraf.imutil(_doprint=0, Stdout=0)
    toforget = ['imutil.imcopy', 'specred.sarith', 'specred.standard']
    for t in toforget:
        iraf.unlearn(t)
    iraf.specred.verbose = 'no'
    iraf.specred.dispaxi = 2
    iraf.set(direc=ntt.__path__[0] + '/')
    sens = {}
    outputfile = []
    if 'obj' in objectlist.keys():
        tpe = 'obj'
    elif 'stdsens' in objectlist.keys():
        tpe = 'stdsens'
    else:
        sys.exit('error: no objects and no standards in the list')

    for setup in objectlist[tpe]:
        extracted = []
        listatmo = []
        if setup not in sens:
            sens[setup] = []
        if tpe == 'obj':
            print('\n### setup= ', setup, '\n### objects= ', objectlist['obj'][setup], '\n')
            for img in objectlist['obj'][setup]:
                #              hdr=readhdr(img)
                print('\n\n### next object= ', img, ' ', ntt.util.readkey3(ntt.util.readhdr(img), 'object'), '\n')
                if os.path.isfile(re.sub('.fits', '_ex.fits', img)):
                    if ntt.util.readkey3(ntt.util.readhdr(re.sub('.fits', '_ex.fits', img)), 'quality') == 'Rapid':
                        ntt.util.delete(re.sub('.fits', '_ex.fits', img))
                imgex = ntt.util.extractspectrum(img, dv, _ext_trace, _dispersionline, _interactive, 'obj',
                                                 automaticex=_automaticex)
                if not os.path.isfile(imgex):
                    sys.exit('### error, extraction not computed')
                if not ntt.util.readkey3(ntt.util.readhdr(imgex), 'shift') and \
                        ntt.util.readkey3(ntt.util.readhdr(imgex), 'shift') != 0.0:
                    # if not readkey3(readhdr(imgex),'shift'):
                    ntt.efoscspec1Ddef.checkwavestd(imgex, _interactive)
                extracted.append(imgex)
                if imgex not in outputfile:
                    outputfile.append(imgex)
                ntt.util.updateheader(
                    imgex, 0, {'FILETYPE': [22107, 'extracted 1D spectrum ']})
                ntt.util.updateheader(imgex, 0, {
                    'PRODCATG': ['SCIENCE.' +
                                 ntt.util.readkey3(ntt.util.readhdr(imgex), 'tech').upper(), 'Data product category']})
                ntt.util.updateheader(
                    imgex, 0, {'TRACE1': [img, 'Originating file']})
                if os.path.isfile('database/ap' + re.sub('_ex.fits', '', imgex)):
                    if 'database/ap' + re.sub('_ex.fits', '', imgex) not in outputfile:
                        outputfile.append(
                            'database/ap' + re.sub('_ex.fits', '', imgex))
            print('\n### all object with this setup extracted\n')
        if liststandard:
            standardlist = liststandard
            _type = 'stdfromdreducer'
        else:
            try:
                standardlist = objectlist['stdsens'][setup]
                _type = 'stdsens'
            except:
                standardlist = ''
                _type = ''
        if _type == 'stdfromdreducer' and len(extracted) >= 1:
            _outputsens2 = ntt.util.searchsens(extracted[0], standardlist)[0]
            print('\n### using standard from reducer ' + str(_outputsens2))
        elif _type not in ['stdsens', 'stdfromdreducer'] and len(extracted) >= 1:
            _outputsens2 = ntt.util.searchsens(extracted[0], '')[0]
            os.system('cp ' + _outputsens2 + ' .')
            _outputsens2 = _outputsens2.split('/')[-1]
            print('\n### no standard in the list, using standard from archive')
        else:
            for simg in standardlist:
                print('\n###  standard for setup ' + \
                      str(setup) + ' = ', simg, ' ', ntt.util.readkey3(
                          ntt.util.readhdr(simg), 'object'), '\n')
                simgex = ntt.util.extractspectrum(
                    simg, dv, False, False, _interactive, 'std', automaticex=_automaticex)
                ntt.util.updateheader(
                    simgex, 0, {'FILETYPE': [22107, 'extracted 1D spectrum']})
                ntt.util.updateheader(simgex, 0, {
                    'PRODCATG': [
                        'SCIENCE.' + ntt.util.readkey3(ntt.util.readhdr(simgex), 'tech').upper(), 'Data product category']})
                ntt.util.updateheader(
                    simgex, 0, {'TRACE1': [simg, 'Originating file']})
                if not ntt.util.readkey3(ntt.util.readhdr(simgex), 'shift') and \
                        ntt.util.readkey3(ntt.util.readhdr(simgex), 'shift') != 0.0:
                    #                if not readkey3(readhdr(simgex),'shift'):
                    ntt.efoscspec1Ddef.checkwavestd(simgex, _interactive)
                atmofile = ntt.efoscspec1Ddef.telluric_atmo(
                    simgex)  # atmo file2
                ntt.util.updateheader(
                    atmofile, 0, {'TRACE1': [simgex, 'Originating file']})
                ntt.util.updateheader(
                    atmofile, 0, {'FILETYPE': [21211, 'telluric correction 1D spectrum ']})
                if tpe != 'obj' and atmofile not in outputfile:
                    outputfile.append(atmofile)
                if not listatmo0:
                    listatmo.append(atmofile)
                sens[setup].append(simgex)
                if simgex not in outputfile:
                    outputfile.append(simgex)
                if setup[0] == 'Gr13' and setup[1] == 'Free':
                    if os.path.isfile(re.sub('Free', 'GG495', simg)):
                        print('\n### extract standard frame with blocking filter to correct for second order contamination\n')
                        simg2 = re.sub('Free', 'GG495', simg)
                        simgex2 = ntt.util.extractspectrum(simg2, dv, False, False, _interactive, 'std',
                                                           automaticex=_automaticex)
                        ntt.util.updateheader(
                            simgex2, 0, {'FILETYPE': [22107, 'extracted 1D spectrum']})
                        ntt.util.updateheader(simgex2, 0, {
                            'PRODCATG': ['SCIENCE.' +
                                         ntt.util.readkey3(
                                             ntt.util.readhdr(simgex2), 'tech').upper(), 'Data product category']})
                        if not ntt.util.readkey3(ntt.util.readhdr(simgex2), 'shift') and \
                                ntt.util.readkey3(ntt.util.readhdr(simgex2), 'shift') != 0.0:
                            # if not readkey3(readhdr(simgex2),'shift'):
                            ntt.efoscspec1Ddef.checkwavestd(
                                simgex2, _interactive)
                        ntt.util.updateheader(
                            simgex2, 0, {'TRACE1': [simg2, 'Originating file']})
            print('\n### standard available: ', sens[setup])
            if tpe == 'obj':
                if len(sens[setup]) > 1:
                    goon = 'no'
                    while goon != 'yes':
                        stdused = input(
                            '\n### more than one standard for this setup, which one do you want to use [' + sens[setup][
                                0] + '] ?')
                        if not stdused:
                            stdused = sens[setup][0]
                        if os.path.isfile(stdused):
                            goon = 'yes'
                else:
                    stdused = sens[setup][0]
                stdvec = [stdused]
            else:
                stdvec = sens[setup]
            for stdused in stdvec:
                stdusedclean = re.sub('_ex', '_clean', stdused)
                ntt.util.delete(stdusedclean)
                iraf.specred.sarith(
                    input1=stdused, op='/', input2=atmofile, output=stdusedclean, format='multispec')
                _outputsens2 = ntt.efoscspec1Ddef.sensfunction(
                    stdusedclean, 'spline3', 16, _interactive)
                ntt.util.updateheader(_outputsens2, 0, {'FILETYPE': [
                                      21212, 'sensitivity function']})
                ntt.util.updateheader(
                    _outputsens2, 0, {'TRACE1': [stdused, 'Originating file']})

                if setup[0] == 'Gr13' and setup[1] == 'Free':
                    if os.path.isfile(re.sub('Free', 'GG495', stdused)):
                        print('\n### compute sensitivity function of grim 13 with blocking filter ' \
                              'to correct for second order contamination \n')
                        stdused2 = re.sub('Free', 'GG495', stdused)
                        if not ntt.util.readkey3(ntt.util.readhdr(stdused2), 'STDNAME'):
                            ntt.util.updateheader(stdused2, 0, {
                                'STDNAME': [ntt.util.readkey3(ntt.util.readhdr(stdused), 'STDNAME'), '']})
                        atmofile2 = ntt.efoscspec1Ddef.telluric_atmo(
                            stdused2)  # atmo file2
                        stdusedclean2 = re.sub('_ex', '_clean', stdused2)
                        ntt.util.delete(stdusedclean2)
                        iraf.specred.sarith(input1=stdused2, op='/', input2=atmofile2, output=stdusedclean2,
                                            format='multispec')
                        _outputsens3 = ntt.efoscspec1Ddef.sensfunction(
                            stdusedclean2, 'spline3', 16, _interactive)
                        ntt.util.updateheader(_outputsens3, 0, {'FILETYPE': [
                                              21212, 'sensitivity function']})
                        ntt.util.updateheader(
                            _outputsens3, 0, {'TRACE1': [stdused2, 'Originating file']})
                        _outputsens2 = correctsens(_outputsens2, _outputsens3)

                if _outputsens2 not in outputfile:
                    outputfile.append(_outputsens2)
        if _outputsens2 and tpe == 'obj':
            ####################################################
            for img in objectlist['obj'][setup]:  # flux calibrate 2d images
                imgd = fluxcalib2d(img, _outputsens2)
                ntt.util.updateheader(
                    imgd, 0, {'FILETYPE': [22209, '2D wavelength and flux calibrated spectrum ']})
                ntt.util.updateheader(
                    imgd, 0, {'TRACE1': [img, 'Originating files']})
                iraf.hedit(imgd, 'PRODCATG', delete='yes',
                           update='yes', verify='no')
                if imgd not in outputfile:
                    outputfile.append(imgd)
            ####################################################
            #    flux calib in the standard way
            if not listatmo and listatmo0:
                listatmo = listatmo0[:]
            for _imgex in extracted:
                _airmass = ntt.util.readkey3(
                    ntt.util.readhdr(_imgex), 'airmass')
                _exptime = ntt.util.readkey3(
                    ntt.util.readhdr(_imgex), 'exptime')
                _imgf = re.sub('_ex.fits', '_f.fits', _imgex)
                ntt.util.delete(_imgf)
                qqq = iraf.specred.calibrate(input=_imgex, output=_imgf, sensiti=_outputsens2, extinct='yes',
                                             flux='yes',
                                             extinction=_extinctdir + _extinction, observatory=_observatory,
                                             airmass=_airmass, ignorea='yes', exptime=_exptime, fnu='no')
                hedvec = {'SENSFUN': [_outputsens2, ''],
                          'FILETYPE': [22208, '1D wavelength and flux calibrated spectrum', ''],
                          #                     'SNR':[ntt.util.StoN(_imgf,50),'Average signal to noise ratio per pixel'],
                          'SNR': [ntt.util.StoN2(_imgf, False), 'Average signal to noise ratio per pixel'],
                          'BUNIT': ['erg/cm2/s/Angstrom', 'Physical unit of array values'],
                          'TRACE1': [_imgex, 'Originating file'],
                          'ASSON1': [re.sub('_f.fits', '_2df.fits', _imgf), 'Name of associated file'],
                          'ASSOC1': ['ANCILLARY.2DSPECTRUM', 'Category of associated file']}
                ntt.util.updateheader(_imgf, 0, hedvec)
                if _imgf not in outputfile:
                    outputfile.append(_imgf)
                if listatmo:
                    atmofile = ntt.util.searcharc(_imgex, listatmo)[0]
                    if atmofile:
                        _imge = re.sub('_f.fits', '_e.fits', _imgf)
                        ntt.util.delete(_imge)
                        iraf.specred.sarith(input1=_imgf, op='/', input2=atmofile, output=_imge, w1='INDEF', w2='INDEF',
                                            format='multispec')
                        try:
                            iraf.imutil.imcopy(
                                input=_imgf + '[*,1,2]', output=_imge + '[*,1,2]', verbose='no')
                        except:
                            pass
                        try:
                            iraf.imutil.imcopy(
                                input=_imgf + '[*,1,3]', output=_imge + '[*,1,3]', verbose='no')
                        except:
                            pass
                        try:
                            iraf.imutil.imcopy(
                                input=_imgf + '[*,1,4]', output=_imge + '[*,1,4]', verbose='no')
                        except:
                            pass
                        if _imge not in outputfile:
                            outputfile.append(_imge)
                        ntt.util.updateheader(
                            _imge, 0, {'FILETYPE': [22210, '1D, wave, flux calib, telluric corr.']})
                        if atmofile not in outputfile:
                            outputfile.append(atmofile)
                        ntt.util.updateheader(
                            _imge, 0, {'ATMOFILE': [atmofile, '']})
                        ntt.util.updateheader(
                            _imge, 0, {'TRACE1': [_imgf, 'Originating file']})
                        imgin = _imge
                    else:
                        imgin = _imgf
                else:
                    imgin = _imgf
                imgasci = re.sub('.fits', '.asci', imgin)

                ntt.util.delete(imgasci)
                iraf.onedspec(_doprint=0, Stdout=0)
                iraf.onedspec.wspectext(
                    imgin + '[*,1,1]', imgasci, header='no')
                if imgasci not in outputfile:
                    outputfile.append(imgasci)

    print('\n### adding keywords for phase 3 ....... ')
    for img in outputfile:
        if str(img)[-5:] == '.fits':
            try:
                ntt.util.phase3header(img)  # phase 3 definitions
                ntt.util.updateheader(img, 0, {'quality': ['Final', '']})
            except:
                print('Warning: ' + img + ' is not a fits file')
            try:
                if int(re.sub('\.', '', str(pyfits.__version__))[:2]) <= 30:
                    aa = 'HIERARCH '
                else:
                    aa = ''
            except:
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

    print(outputfile)
    reduceddata = ntt.rangedata(outputfile)
    f = open('logfile_spec1d_' + str(reduceddata) +
             '_' + str(datenow) + '.raw.list', 'w')
    for img in outputfile:
        try:
            f.write(ntt.util.readkey3(ntt.util.readhdr(img), 'arcfile') + '\n')
        except:
            pass
    f.close()
    return outputfile, 'logfile_spec1d_' + str(reduceddata) + '_' + str(datenow) + '.raw.list'


##########################################################################

def correctsens(img1, img2):

    import os
    import re

    try:        import pyfits
    except:     from astropy.io import fits as pyfits

    import numpy as np

    # read spectrum 1
    yys = pyfits.open(img1)[0].data
    crvals = pyfits.open(img1)[0].header.get('CRVAL1')
    cds = pyfits.open(img1)[0].header.get('CD1_1')
    xxs = np.arange(len(yys))
    aasens = crvals + (xxs) * cds
    #    read spectrum 2
    yysgg = pyfits.open(img2)[0].data
    crvalsgg = pyfits.open(img2)[0].header.get('CRVAL1')
    cdsgg = pyfits.open(img2)[0].header.get('CD1_1')
    xxsgg = np.arange(len(yysgg))
    aasensgg = crvalsgg + (xxsgg) * cdsgg

    # scale Gr13GG to Gr13
    yysgg2 = yysgg * np.interp(5500, aasens, yys) / \
        np.interp(5500, aasensgg, yysgg)
    yysgg2cut = np.compress((aasensgg > 5500), yysgg2)
    aasensggcut = np.compress((aasensgg > 5500), aasensgg)
    yysgg2cutblu = np.compress(aasens <= 5500, yys)  # blue part
    if aasens[-1] > aasensgg[-1]:
        aasenscut = np.compress((aasens > 5500) & (
            aasens < aasensgg[-1]), aasens)
        yysgg2cutred = np.interp(
            aasenscut, aasensggcut, yysgg2cut)  # medium part
        yysggred2 = np.compress((aasens >= aasensggcut[-1]), yys) + (
            yysgg2cutred[-1] - np.compress((aasens >= aasensggcut[-1]), yys)[0])  # red part
        finale = np.array(list(yysgg2cutblu) +
                          list(yysgg2cutred) + list(yysggred2))
    else:
        aasenscut = np.compress((aasens > 5500), aasens)
        yysgg2cutred = np.interp(
            aasenscut, aasensggcut, yysgg2cut)  # medium part
        finale = np.array(list(yysgg2cutblu) + list(yysgg2cutred))
    data, hdr = pyfits.getdata(img1, 0, header=True)
    data = finale
    imgout = re.sub('.fits', '_2ord.fits', img1)
    if os.path.isfile(imgout):
        os.remove(imgout)
    pyfits.writeto(imgout, np.float32(data), hdr)
    return imgout
##########################################################################
