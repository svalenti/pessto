def makefringingmask(listimg, _output, _interactive, _combine='average', _rejection='avsigclip'):
    # print "LOGX:: Entering `makefringingmask` method/function in
    # %(__file__)s" % globals()
    import ntt
    from ntt.util import readhdr, readkey3, delete, updateheader
    import glob
    import os
    import sys
    import re
    import string
    from pyraf import iraf
    iraf.noao(_doprint=0)
    iraf.immatch(_doprint=0)
    iraf.imutil(_doprint=0)
    iraf.nproto(_doprint=0)
    iraf.proto(_doprint=0)
    toforget = ['nproto.objmasks', 'proto.fixpix']
    for t in toforget:
        iraf.unlearn(t)

    if _interactive == True:
        listimg2 = []
        for img in listimg:
            _exptime = readkey3(readhdr(img), 'exptime')
            if float(_exptime) >= 10:
                answ = 'xxx'
                while answ.lower() not in ['y', 'n', 's', 'a']:
                    iraf.display(img, frame=1, fill='yes')
                    answ = raw_input(
                        'use this image (yes,no,stop (not more images),all) [[y]/n/s/a] ? ')
                    if not answ:
                        answ = 'y'
                    if answ.lower() == 'y':
                        listimg2.append(img)
                    elif answ.lower() == 'a':
                        listimg2 = listimg[:]
                if answ.lower() in ['a', 's']:
                    break
        listimg = listimg2[:]

    iraf.nproto.objmasks1.fitxord = 1
    iraf.nproto.objmasks1.fityord = 1
    hdr0 = readhdr(listimg[0])
    _date = readkey3(hdr0, 'date-obs')
    _filter = readkey3(hdr0, 'filter')
    _exptime = readkey3(hdr0, 'exptime')
    _instrume = readkey3(hdr0, 'instrume')
    _ron = readkey3(hdr0, 'ron')
    _gain = readkey3(hdr0, 'gain')
    badpixelmask = 'bad_pixel_mask.pl'
    if not os.path.isfile(badpixelmask):
        os.system('cp ' + ntt.__path__[0] + '/archive/' + _instrume +
                  '/badpixels/badpixel_20100210.pl ' + badpixelmask)
    ff = open('_listmask', 'w')
    hh = open('_listobgz', 'w')
    for img in listimg:
        _exptime = readkey3(readhdr(img), 'exptime')
        hh.write('z_' + img + '\n')
        ff.write('mask_' + img + '\n')
        delete('mask_' + img)
        aaa = iraf.hedit(img, delete='yes', field='OBJMASK',
                         up='yes', verify='no', Stdout=1)
        aaa = iraf.hedit(img, delete='yes', field='BPM',
                         up='yes', verify='no', Stdout=1)
        delete('z_' + img)
        iraf.imutil.imexpr(expr='(a - median(a))/' + str(_exptime),
                           a=img, output='z_' + img, verbose='no')
        ntt.util.updateheader('z_' + img, 0, {'EXPTIME': [1, '']})
    ff.close()
    hh.close()
    if not _output:
        _output = 'fringing_' + str(_date) + '_' + str(_filter) + '.fits'
    delete(_output)
    print ' making mask for each frame .......'
    ccc = iraf.nproto.objmasks(images='@_listobgz', objmasks='@_listmask', omtype='boolean',
                               blksize=-16, convolv='block 3 3', hsigma=5, lsigma=3, minpix=10, ngrow=2, agrow=4., Stdout=1)
    print 'combining all frames, masking the objects .....'
    iraf.imcombine('@_listobgz', output=_output, masktyp='!OBJMASK', maskval=0, combine=_combine, reject=_rejection,
                   scale='none', statsec='[100:800,100:800]', rdnoise='', gain='', nlow=1, nhigh=1, logfile='imcombinelog')

    ntt.util.phase3header(_output)
    ntt.util.updateheader(
        _output, 0, {'BUNIT': ['ADU', 'pixel units(ADU,electrons)']})
    ntt.util.updateheader(_output, 0, {'FILETYPE': [11231, 'fringing frame']})
    return _output

##########################################################################


def makefringing(listimg, _output, _xorder, _yorder, _interactive, combine='average', rejection='avsigclip'):
    # print "LOGX:: Entering `makefringing` method/function in %(__file__)s" %
    # globals()
    import os
    import string
    import re
    import ntt
    from ntt.util import readhdr, readkey3, delete, updateheader
    from pyraf import iraf
    from pyfits import open as popen
    from numpy import median, where, mean
    iraf.noao(_doprint=0)
    iraf.nproto(_doprint=0)
    if _interactive == True:
        listimg2 = []
        for img in listimg:
            _exptime = readkey3(readhdr(img), 'exptime')
            if float(_exptime) >= 100:
                answ = 'xxx'
                while answ.lower() not in ['y', 'n', 's', 'a']:
                    iraf.display(img, frame=1, fill='yes')
                    answ = raw_input(
                        'use this image (yes,no,stop (not more images),all) [[y]/n/s/a] ? ')
                    if not answ:
                        answ = 'y'
                    if answ.lower() == 'y':
                        listimg2.append(img)
                    elif answ.lower() == 'a':
                        listimg2 = listimg[:]
                if answ.lower() in ['a', 's']:
                    break
        listimg = listimg2[:]

    iraf.nproto.objmasks1.fitxord = 1
    iraf.nproto.objmasks1.fityord = 1
    listmask = []
    ff = open('_listmask', 'w')
    gg = open('_listobg', 'w')
    hh = open('_listobgz', 'w')
    for img in listimg:
        _exptime = readkey3(readhdr(img), 'exptime')
        gg.write(img + '\n')
        hh.write('z_' + img + '\n')
        ff.write('mask_' + img + '\n')
        listmask.append('mask_' + img)
        os.system('rm -rf mask_' + img)
        os.system('rm -rf t_' + img)
        os.system('rm -rf z_' + img)

        iraf.imsurfit(img, 't_' + img, xorder=_xorder,
                      yorder=_yorder, type_ou='residual', regions='all')
        iraf.imarith('t_' + img, '/', _exptime, 'z_' + img)
        os.system('rm -rf t_' + img)
    gg.close()
    ff.close()
    hh.close()

    os.system('rm ' + _output)
    iraf.objmasks(images='@_listobgz', objmasks='@_listmask', omtype='boolean',
                  blksize=-16, convolv='block 3 3', hsigma=5, lsigma=3, minpix=10, ngrow=2, agrow=4.)
    iraf.imcombine('@_listobgz', output=_output, masktyp='!OBJMASK', maskval=0, combine=_combine, reject=_rejection,
                   scale='none', statsec='[100:800,100:800]', offsets='', rdnoise='', gain='', nlow=1, nhigh=1, logfile='imcombinelog')

    for img in listimg:
        os.system('rm -rf mask_' + img)
        os.system('rm -rf z_' + img)

    os.system('rm _listmask')
    os.system('rm _listobg')
    iraf.display(_output, frame=2, fill='yes')

##########################################################################
