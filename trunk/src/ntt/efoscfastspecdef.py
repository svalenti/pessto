
def efoscfastredu(imglist, _listsens, _listarc, _ext_trace, _dispersionline, _cosmic, _interactive):
    import string
    import os
    import re
    import sys
    os.environ["PYRAF_BETA_STATUS"] = "1"
    import pyfits
    from ntt.util import readhdr, readkey3
    import ntt
    import numpy as np
    dv = ntt.dvex()
    scal = np.pi / 180.
    if not _interactive:
        _interactive = False
        _inter = 'NO'
    else:
        _inter = 'YES'
    from pyraf import iraf

    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.ccdred(_doprint=0)
    iraf.twodspec(_doprint=0)
    iraf.longslit(_doprint=0)
    iraf.onedspec(_doprint=0)
    iraf.specred(_doprint=0)
    toforget = ['ccdproc', 'imcopy', 'specred.apall', 'longslit.identify', 'longslit.reidentify', 'specred.standard',
                'longslit.fitcoords', 'onedspec.wspectext']
    for t in toforget:
        iraf.unlearn(t)
    iraf.ccdred.verbose = 'no'  # not print steps
    iraf.specred.verbose = 'no'  # not print steps
    iraf.ccdproc.darkcor = 'no'
    iraf.ccdproc.fixpix = 'no'
    iraf.ccdproc.flatcor = 'no'
    iraf.ccdproc.zerocor = 'no'
    iraf.ccdproc.ccdtype = ''
    _gain = ntt.util.readkey3(ntt.util.readhdr(imglist[0]), 'gain')
    _ron = ntt.util.readkey3(ntt.util.readhdr(imglist[0]), 'ron')
    iraf.specred.apall.readnoi = _ron
    iraf.specred.apall.gain = _gain
    iraf.specred.dispaxi = 2
    iraf.longslit.dispaxi = 2
    iraf.longslit.mode = 'h'
    iraf.specred.mode = 'h'
    iraf.noao.mode = 'h'
    iraf.ccdred.instrument = "ccddb$kpno/camera.dat"
    iraf.set(direc=ntt.__path__[0] + '/')
    for img in imglist:
        hdr = ntt.util.readhdr(img)
        _tech = ntt.util.readkey3(hdr, 'tech')
        if _tech != 'SPECTRUM': sys.exit('error: ' + str(img) + ' is not a spectrum ')
        print '\n####  image name = ' + img + '\n'
        _grism0 = readkey3(hdr, 'grism')
        _filter0 = readkey3(hdr, 'filter')
        _slit0 = readkey3(hdr, 'slit')
        _object0 = readkey3(hdr, 'object')
        _date0 = readkey3(hdr, 'date-night')
        setup = (_grism0, _filter0, _slit0)
        _biassec0 = '[3:1010,1026:1029]'
        if _grism0 == 'Gr16':
            _trimsec0 = '[100:950,1:950]'
        elif _grism0 == 'Gr13':
            if _filter0 == 'Free':
                _trimsec0 = '[100:950,1:1015]'
            elif _filter0 == 'GG495':
                _trimsec0 = '[100:950,208:1015]'
            elif _filter0 == 'OG530':
                _trimsec0 = '[100:950,300:1015]'
        elif _grism0 == 'Gr11':
            _trimsec0 = '[100:950,5:1015]'
        else:
            _trimsec0 = '[100:950,5:1015]'
        _object0 = re.sub(' ', '', _object0)
        _object0 = re.sub('/', '_', _object0)
        nameout0 = 't' + str(_object0) + '_' + str(_date0)
        for _set in setup:
            nameout0 = nameout0 + '_' + _set
        nameout0 = ntt.util.name_duplicate(img, nameout0, '')
        timg = nameout0
        if os.path.isfile(timg): os.system('rm -rf ' + timg)
        iraf.imcopy(img, output=timg)
        iraf.ccdproc(timg, output='', overscan='no', trim='yes', zerocor="no", flatcor="no", readaxi='column',
                     trimsec=str(_trimsec0), biassec=_biassec0, Stdout=1)
        img = timg
        if _listarc:
            arcfile = ntt.util.searcharc(img, _listarc)[0]
        else:
            arcfile = ''
        if not arcfile:
            arcfile = ntt.util.searcharc(img, '')[0]
        else:
            iraf.ccdproc(arcfile, output='t' + arcfile, overscan='no', trim='yes', zerocor="no", flatcor="no",
                         readaxi='column', trimsec=str(_trimsec0), biassec=str(_biassec0), Stdout=1)
            arcfile = 't' + arcfile

        if _cosmic:
            # print cosmic rays rejection
            ntt.cosmics.lacos(img, output='', gain=_gain, readn=_ron, xorder=9, yorder=9, sigclip=4.5, sigfrac=0.5,
                              objlim=1, verbose=True, interactive=False)
            print '\n### cosmic rays rejections ........ done '

        if not arcfile:
            print '\n### warning no arcfile \n exit '
        else:
            arcref = ntt.util.searcharc(img, '')[0]
            if arcfile[0] == '/':
                os.system('cp ' + arcfile + ' ' + string.split(arcfile, '/')[-1])
                arcfile = string.split(arcfile, '/')[-1]
            arcref = string.split(arcref, '/')[-1]
            if arcref:
                os.system('cp ' + arcref + ' .')
                arcref = string.split(arcref, '/')[-1]
                if not os.path.isdir('database/'):   os.mkdir('database/')
                if os.path.isfile(ntt.util.searcharc(img, '')[1] + '/database/id' + re.sub('.fits', '', arcref)):
                    os.system('cp ' + ntt.util.searcharc(img, '')[1] + '/database/id' + re.sub('.fits', '',
                                                                                               arcref) + ' database/')
                iraf.longslit.reidentify(referenc=arcref, images=arcfile, interac=_inter, section='column 10',
                                         coordli='direc$standard/ident/Lines_HgCdHeNeAr600.dat', overrid='yes', step=0,
                                         newaps='no', nsum=5, nlost=2, mode='h', verbose='no')
            else:
                iraf.longslit.identify(images=arcfile, section='column 10',
                                       coordli='direc$standard/ident/Lines_HgCdHeNeAr600.dat', nsum=10, fwidth=7,
                                       order=3, mode='h')
            iraf.longslit.reident(referenc=arcfile, images=arcfile, interac='NO', section='column 10',
                                  coordli='direc$standard/ident/Lines_HgCdHeNeAr600.dat', overrid='yes', step=10,
                                  newaps='yes', nsum=5, nlost=2, mode='h', verbose='no')
            qqq = iraf.longslit.fitcoords(images=re.sub('.fits', '', arcfile), fitname=re.sub('.fits', '', arcfile),
                                          interac='no', combine='yes', databas='database',
                                          function='legendre', yorder=4, logfile='logfile', plotfil='', mode='h')
            iraf.specred.transform(input=img, output=img, minput='', fitnames=re.sub('.fits', '', arcfile),
                                   databas='database',
                                   x1='INDEF', x2='INDEF', y1='INDEF', y2='INDEF', flux='yes', mode='h',
                                   logfile='logfile')
            # ######################  check wavelength calibration ###################
            _skyfile = ntt.__path__[0] + '/standard/ident/sky_' + setup[0] + '_' + setup[1] + '.fits'
            shift = ntt.efoscspec2Ddef.skyfrom2d(img, _skyfile)
            print '\n###     check in wavelengh performed ...... spectrum shifted of  ' + str(shift) + ' Angstrom \n'
            zro = pyfits.open(img)[0].header.get('CRVAL2')
            ntt.util.updateheader(img, 0, {'CRVAL2': [zro + int(shift), '']})
            std, rastd, decstd, magstd = ntt.util.readstandard('standard_efosc_mab.txt')
            hdrt = readhdr(img)
            _ra = readkey3(hdrt, 'RA')
            _dec = readkey3(hdrt, 'DEC')
            _object = readkey3(hdrt, 'object')
            dd = np.arccos(np.sin(_dec * scal) * np.sin(decstd * scal) + np.cos(_dec * scal) *
                           np.cos(decstd * scal) * np.cos((_ra - rastd) * scal)) * ((180 / np.pi) * 3600)
            if min(dd) < 100:
                _type = 'stdsens'
                ntt.util.updateheader(img, 0, {'stdname': [std[np.argmin(dd)], '']})
                ntt.util.updateheader(img, 0, {'magstd': [float(magstd[np.argmin(dd)]), '']})
            else:
                _type = 'obj'
            print '\n###      EXTRACTION USING IRAF TASK APALL \n'
            result = []
            if _type == 'obj':
                imgex = ntt.util.extractspectrum(img, dv, _ext_trace, _dispersionline, _interactive, _type)
                ntt.util.updateheader(imgex, 0, {'FILETYPE': [22107, 'extracted 1D spectrum ']})
                ntt.util.updateheader(imgex, 0, {
                'PRODCATG': ['SCIENCE.' + readkey3(readhdr(imgex), 'tech').upper(), 'Data product category']})
                ntt.util.updateheader(imgex, 0, {'TRACE1': [img, '']})
                result.append(imgex)
                if _listsens:
                    sensfile = ntt.util.searchsens(img, _listsens)[0]
                else:
                    sensfile = ''
                if not sensfile:      sensfile = ntt.util.searchsens(img, '')[0]
                if sensfile:
                    imgf = re.sub('.fits', '_f.fits', img)
                    _extinctdir = 'direc$standard/extinction/'
                    _extinction = 'extinction_lasilla.dat'
                    _observatory = 'lasilla'
                    _exptime = readkey3(hdrt, 'exptime')
                    _airmass = readkey3(hdrt, 'airmass')
                    ntt.util.delete(imgf)
                    iraf.specred.calibrate(input=imgex, output=imgf, sensiti=sensfile, extinct='yes',
                                           flux='yes', ignorea='yes', extinction=_extinctdir + _extinction,
                                           observatory=_observatory, airmass=_airmass, exptime=_exptime,
                                           fnu='no')
                    hedvec = {'SENSFUN': [string.split(sensfile, '/')[-1], 'sensitivity function'],
                              'FILETYPE': [22208, '1D wavelength and flux calibrated spectrum '],
                              'SNR': [ntt.util.StoN2(imgf, False), 'Average S/N ratio'],
                              'BUNIT': ['erg/cm2/s/A', 'Flux Calibration Units'], 'TRACE1': [imgex, '']}
                    ntt.util.updateheader(imgf, 0, hedvec)
                    imgout = imgf
                    imgd = ntt.efoscspec1Ddef.fluxcalib2d(img, sensfile)
                    ntt.util.updateheader(imgd, 0, {'FILETYPE': [22209, '2D wavelength and flux calibrated spectrum ']})
                    ntt.util.updateheader(imgd, 0, {'TRACE1': [img, '']})
                    imgasci = re.sub('.fits', '.asci', imgout)
                    ntt.util.delete(imgasci)
                    iraf.onedspec.wspectext(imgout + '[*,1,1]', imgasci, header='no')
                    result = result + [imgout, imgd, imgasci]
            else:
                imgex = ntt.util.extractspectrum(img, dv, _ext_trace, _dispersionline, _interactive, 'std')
                imgout = ntt.efoscspec1Ddef.sensfunction(imgex, 'spline3', 6, _inter)
                result = result + [imgout]

    for img in result:
        if img[-5:] == '.fits':
            ntt.util.phase3header(img)  # phase 3 definitions
            ntt.util.airmass(img)  # phase 3 definitions
            ntt.util.updateheader(img, 0, {'quality': ['Rapid', 'Final or Rapid reduction']})
    return result
