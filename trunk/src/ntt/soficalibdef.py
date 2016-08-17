def makeflat(lista):
    # print "LOGX:: Entering `makeflat` method/function in %(__file__)s" %
    # globals()
    flat = ''
    import datetime
    import glob
    import os
    import ntt
    from ntt.util import readhdr, readkey3, delete, name_duplicate, updateheader, correctcard
    from pyraf import iraf
    iraf.images(_doprint=0)
    iraf.imutil(_doprint=0)
    iraf.imgeom(_doprint=0)
    # iraf.blkavg(_doprint=0)
    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.generic(_doprint=0)
    toforget = ['imgeom.blkavg', 'imutil.imarith',
                'immatch.imcombine', 'noao.imred']
    for t in toforget:
        iraf.unlearn(t)
    import datetime
    MJDtoday = 55927 + (datetime.date.today() - datetime.date(2012, 01, 01)).days
    _date = readkey3(readhdr(lista[0]), 'date-night')
    _filter = readkey3(readhdr(lista[0]), 'filter')
    output = name_duplicate(
        lista[3], 'flat_' + str(_date) + '_' + str(_filter) + '_' + str(MJDtoday), '')
    if os.path.isfile(output):
        answ = raw_input('file already prooduced, do again [y/[n]] ? ')
        if not answ:
            answ = 'n'
    else:
        answ = 'y'
    if answ in ['yes', 'y', 'YES', 'Y', 'Yes']:
        delete("temp_off.fits,temp_off_mask.fits,temp_on_mask.fits,temp_on.fits")
        iraf.image.immatch.imcombine(
            lista[0] + ',' + lista[7], output="temp_off.fits")
        iraf.image.immatch.imcombine(
            lista[1] + ',' + lista[6], output="temp_off_mask.fits")
        iraf.image.immatch.imcombine(
            lista[2] + ',' + lista[5], output="temp_on_mask.fits")
        iraf.image.immatch.imcombine(
            lista[3] + ',' + lista[4], output="temp_on.fits")
        #   create the bias correction for the flat-on according to the
        #   Lidman technique0
        delete("temp_onA.fits,temp_onC.fits,temp_onB.fits,temp_onAC.fits,temp_onACB.fits,temp_onACB_2D.fits")
        delete("temp_on_bias.fits")
        iraf.imgeom.blkavg(
            input="temp_on.fits[500:600,*]", output="temp_onA.fits", option="average", b1=101, b2=1)
        iraf.imgeom.blkavg(
            input="temp_on_mask.fits[500:600,*]", output="temp_onC.fits", option="average", b1=101, b2=1)
        iraf.imgeom.blkavg(
            input="temp_on_mask.fits[50:150,*]", output="temp_onB.fits", option="average", b1=101, b2=1)
        iraf.imutil.imarith("temp_onA.fits", "-",
                            "temp_onC.fits", "temp_onAC.fits")
        iraf.imutil.imarith("temp_onAC.fits", "+",
                            "temp_onB.fits", "temp_onACB.fits")
        iraf.imgeom.blkrep(input="temp_onACB.fits",
                           output="temp_onACB_2D.fits", b1=1024, b2=1)
        iraf.imutil.imarith("temp_on.fits", "-",
                            "temp_onACB_2D.fits", "temp_on_bias.fits")
    #   same as above for the flat-off
        delete("temp_offA.fits,temp_offC.fits,temp_offB.fits,temp_offAC.fits,temp_offACB.fits,temp_offACB_2D.fits")
        delete("temp_off_bias.fits")
        iraf.imgeom.blkavg(
            input="temp_off.fits[500:600,*]", output="temp_offA.fits", option="average", b1=101, b2=1)
        iraf.imgeom.blkavg(
            input="temp_off_mask.fits[500:600,*]", output="temp_offC.fits", option="average", b1=101, b2=1)
        iraf.imgeom.blkavg(
            input="temp_off_mask.fits[50:150,*]", output="temp_offB.fits", option="average", b1=101, b2=1)
        iraf.imutil.imarith("temp_offA.fits", "-",
                            "temp_offC.fits", "temp_offAC.fits")
        iraf.imutil.imarith("temp_offAC.fits", "+",
                            "temp_offB.fits", "temp_offACB.fits")
        iraf.imgeom.blkrep(input="temp_offACB.fits",
                           output="temp_offACB_2D.fits", b1=1024, b2=1)
        iraf.imutil.imarith("temp_off.fits", "-",
                            "temp_offACB_2D.fits", "temp_off_bias.fits")
        #   create the corrected flat-field
        #    output=name_duplicate("temp_on_bias.fits",'flat_'+str(_date)+'_'+str(_filter)+'_'+str(MJDtoday),'')
        output = name_duplicate(
            lista[3], 'flat_' + str(_date) + '_' + str(_filter) + '_' + str(MJDtoday), '')
    #    print lista[0],'flat_'+str(_date)+'_'+str(_filter)+'_'+str(MJDtoday)
        delete(output)
        iraf.imutil.imarith("temp_on_bias.fits", "-",
                            "temp_off_bias.fits", output)
        iraf.noao.imred.generic.normalize(output)  # normalize the flat-field
        correctcard(output)
        delete("temp_on*.fits")  # delete the temporary images
        delete("temp_off*.fits")
        print 'flat -> ' + str(output)
    else:
        print 'skip redoing the flat'
    return output


def makeillumination(lista, flatfield):  # ,outputfile,illum_frame):
    # print "LOGX:: Entering `makeillumination` method/function in
    # %(__file__)s" % globals()
    import os
    import glob
    import string
    import re
    import pyfits
    import ntt
    from ntt.util import readhdr, readkey3, delete, display_image, defsex,  name_duplicate, correctcard
    from numpy import compress, array, argmax, argmin, min, argsort, float32
    import datetime
    MJDtoday = 55927 + (datetime.date.today() - datetime.date(2012, 01, 01)).days
    _date = readkey3(readhdr(lista[0]), 'date-night')
    _filter = readkey3(readhdr(lista[0]), 'filter')
    illum_frame = name_duplicate(
        lista[0], 'illum_' + _date + '_' + _filter + '_' + str(MJDtoday), '')
    from pyraf import iraf
    iraf.images(_doprint=0)
    iraf.imutil(_doprint=0)
    iraf.utilities(_doprint=0)
    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.ccdred(_doprint=0)
    iraf.digiphot(_doprint=0)
    iraf.daophot(_doprint=0)
    iraf.generic(_doprint=0)
    toforget = ['digiphot.daophot', 'imutil.imarith',
                'image', 'utilities.surfit']
    for t in toforget:
        iraf.unlearn(t)
    n = len(lista)
#   start loop to read image names from the input file
    lista1 = []
    iraf.ccdred.verbose = 'no'
    ff = open('templist.lst', 'w')
    for i in range(0, len(lista)):
        ff.write('C' + lista[i] + '\n')
        delete('C' + lista[i])
        delete('C' + re.sub('.fits', '_sub.fits', lista[i]))
        ntt.sofiphotredudef.crosstalk(lista[i], 'C' + lista[i])
        iraf.noao.imred.ccdred.ccdproc('C' + lista[i], output='', overscan="no", trim="yes", ccdtype='', darkcor='no', fixpix='no', zerocor="no", flatcor='yes',
                                       illumco='no', trimsec='[1:1024,1:1007]', biassec='', flat=flatfield, illum='')
        correctcard('C' + lista[i])
        lista1.append('C' + lista[i])
    ff.close()
    print '\n### prereducing STD frames to compute illumination correction ........'
    lista2, skyfile = ntt.sofiphotredudef.skysub(lista1, readkey3(
        readhdr(lista1[0]), 'ron'), readkey3(readhdr(lista1[0]), 'gain'), True)
    lista2 = ntt.sofiphotredudef.sortbyJD(lista2)
    print '\n### use x on the star and q  to continue....'
    display_image(lista2[0], 2, '', '', False)
    delete('tmpone.coo')
    iraf.image.tv.imexamine(lista2[0], 2, logfile='tmpone.coo',
                            keeplog='yes', xformat='', yformat='', wcs='logical')
    iraf.tvmark(2, 'tmpone.coo', mark="circle", number='yes',
                label='no', radii=8, nxoffse=5, nyoffse=5, color=204, txsize=2)
    xycoo = iraf.proto.fields('tmpone.coo', '1,2', Stdout=1)
    x0, y0 = string.split(xycoo[0])
    x0 = float(x0)
    y0 = float(y0)
    xcum0 = readkey3(readhdr(lista2[0]), 'xcum')
    ycum0 = readkey3(readhdr(lista2[0]), 'ycum')
    iraf.digiphot(_doprint=0)
    iraf.daophot(_doprint=0)
    iraf.noao.digiphot.daophot.datapars.datamin = -1000
    iraf.noao.digiphot.daophot.datapars.datamax = 60000
    iraf.noao.digiphot.daophot.daopars.function = 'gauss'
    iraf.noao.digiphot.daophot.photpars.zmag = 0
    namesex = defsex('default.sex')
    for i in range(0, len(lista2)):
        j = i + 1
        xcum = readkey3(readhdr(lista2[i]), 'xcum')
        ycum = readkey3(readhdr(lista2[i]), 'ycum')
        xx = x0 - xcum0 + xcum
        yy = y0 - ycum0 + ycum
        # sex objects
        os.system('sex ' + lista2[i] + ' -c ' + namesex + '>  _logsex')
        delete('_logsex')
        xpix = iraf.proto.fields('detections.cat', fields='2', Stdout=1)
        ypix = iraf.proto.fields('detections.cat', fields='3', Stdout=1)
        cm = iraf.proto.fields('detections.cat', fields='4', Stdout=1)
        cm = compress((array(xpix) != ''), array(cm, float))
        ypix = compress((array(xpix) != ''), array(ypix, float))
        xpix = compress((array(xpix) != ''), array(xpix, float))
        if len(xpix) > 300:
            num = 300
        else:
            num = len(xpix) - 1
        xpix = xpix[argsort(cm)][0:num]
        ypix = ypix[argsort(cm)][0:num]
        distance = (ypix - yy)**2 + (xpix - xx)**2
        xx1, yy1 = xpix[argmin(distance)], ypix[argmin(distance)]
        f = open('tmpone.coo', 'w')
        f.write(str(xx1) + ' ' + str(yy1) + '\n')
        f.close()
        display_image(lista2[i], 1, '', '', False)
        iraf.tvmark(1, 'tmpone.coo', mark="circle", number='yes',
                    label='no', radii=8, nxoffse=5, nyoffse=5, color=204, txsize=2)
        answ = 'n'
        while answ != 'y':
            answ = raw_input('selected the right one [[y]/n] ?')
            if not answ:
                answ = 'y'
            if answ in ['y', 'YES', 'yes', 'Y']:
                print lista2[i]
                delete('pippo.' + str(j) + '.mag')
                gggg = iraf.digiphot.daophot.phot(
                    lista2[i], "tmpone.coo", output="pippo." + str(j) + ".mag", verify='no', interac='no', Stdout=1)
                try:
                    float(string.split(gggg[0])[3])
                    answ = 'y'
                except:
                    print '\n### warning'
                    answ = 'n'
            else:
                print '\n### select the std star'
                display_image(lista2[i], 1, '', '', False)
                iraf.image.tv.imexamine(lista2[
                                        i], 1, logfile='tmpone.coo', keeplog='yes', xformat='', yformat='', wcs='logical')
                xycoo = iraf.proto.fields('tmpone.coo', '1,2', Stdout=1)
                x2, y2 = string.split(xycoo[0])
                f = open('tmpone.coo', 'w')
                f.write(str(x2) + ' ' + str(y2) + '\n')
                f.close()
                delete('pippo.' + str(j) + '.mag')
                print '###### new selection ' + str(x2), str(y2)
                gggg = iraf.digiphot.daophot.phot(
                    lista2[i], "tmpone.coo", output='pippo.' + str(j) + '.mag', verify='no', interac='no', Stdout=1)
                try:
                    float(string.split(gggg[0])[3])
                    answ = 'y'
                except:
                    print '\n### warning'
                    answ = 'n'

    os.system('ls pippo.*.mag > tempmag.lst')
    tmptbl0 = iraf.txdump(textfile="@tempmag.lst",
                          fields="XCENTER,YCENTER,FLUX", expr='yes', Stdout=1)
    ff = open('magnitudini', 'w')
    for i in tmptbl0:
        ff.write(i + '\n')
    ff.close()
#   delete the temporary images and files
    delete("temp*.fits")
    delete('temp*.lst')
    delete(illum_frame)
    print '\n### fitting the illumination surface...'
    aaa = iraf.utilities.surfit('magnitudini', image=illum_frame, function="polynomial",
                                xorder=2, yorder=2, xterms="full", ncols=1024, nlines=1024, Stdout=1)
    iraf.noao.imred.generic.normalize(illum_frame)
    correctcard(lista[0])
    data, hdr = pyfits.getdata(illum_frame, 0, header=True)
    data0, hdr0 = pyfits.getdata(lista[0], 0, header=True)
    delete(illum_frame)
    pyfits.writeto(illum_frame, float32(data), hdr0)
    flatfield0 = string.split(flatfield, '/')[-1]
    ntt.util.updateheader(
        illum_frame, 0, {'MKILLUM': [flatfield0, 'flat field']})
    display_image(illum_frame, 1, '', '', False)
    for i in range(0, len(lista)):  # in lista:
        img = lista[i]
        delete('pippo.' + str(i) + '.mag')
        delete('C' + img)
        delete('C' + re.sub('.fits', '_sky.fits', img))
#    delete('C*.fits.mag.1')
#    iraf.hedit(illum_frame,'MKILLUM','Illum. corr. created '+flatfield,add='yes',update='yes',verify='no')
    return illum_frame

###############################################################################
    #  select files


def doflatsofi(flats, _doflat, illum, _output):
    # print "LOGX:: Entering `doflatsofi` method/function in %(__file__)s" %
    # globals()
    import ntt
    from ntt.util import display_image, delete, searchflat, name_duplicate
    from pyraf import iraf
    import glob
    import string
    onofflimit = {'J': 1000, 'H': 1000, 'Ks': 5000}
    masklimit = {'J': {'ON': 1000, 'OFF': 30}, 'H': {
        'ON': 1000, 'OFF': 40}, 'Ks': {'ON': 1000, 'OFF': 1000}}
    if flats and _doflat:
        listaflat = []
        for _filter in flats:
            for ID in flats[_filter]:
                images = flats[_filter][ID]
                if len(images) == 8:
                    mflat = makeflat(images)
                    listaflat.append(mflat)
                    display_image(mflat, 1, '', '', False)
                    raw_input('go on ')
                elif len(images) != 8:  # % 8 == 0:
                    print '\n###  to compute a flat field you need a sequence of 8 calibration files in the following orders:'
                    print 'OFF  OFFMASK  ONMASK  ON  ON   ONMASK   OFFMASK     OFF\n'
                    print len(images), _filter, ID
                    tipo = ['OFF', 'OFFMASK', 'ONMASK', 'ON',
                            'ON', 'ONMASK', 'OFFMASK', 'OFF']
                    listtmp = []
                    ii = 0
                    nn = 0
                    for img in images:
                        onoffvalue = float(string.split(iraf.imstat(
                            img + '[500:600,900:1000]', Stdout=1)[1])[2])
                        maskvalue = float(string.split(iraf.imstat(
                            img + '[100:200,900:1000]', Stdout=1)[1])[2])
                        if onoffvalue >= onofflimit[_filter]:
                            onoff = 'ON'
                        else:
                            onoff = 'OFF'
                        if maskvalue >= masklimit[_filter][onoff]:
                            mask = 'none'
                        else:
                            mask = 'MASK'
#                        display_image(img,1,'','',False)
                        print onoff, mask, onoffvalue, maskvalue, img, tipo[nn]
                    for img in images:
                        onoffvalue = float(string.split(iraf.imstat(
                            img + '[500:600,900:1000]', Stdout=1)[1])[2])
                        maskvalue = float(string.split(iraf.imstat(
                            img + '[100:200,900:1000]', Stdout=1)[1])[2])
                        if onoffvalue >= onofflimit[_filter]:
                            onoff = 'ON'
                        else:
                            onoff = 'OFF'
                        if maskvalue >= masklimit[_filter][onoff]:
                            mask = 'none'
                        else:
                            mask = 'MASK'
                        display_image(img, 1, '', '', False)
                        print onoff, mask, onoffvalue, maskvalue, img, tipo[nn]
                        answ = raw_input('ok  [[y]/n/r/s] ? ')
                        if not answ:
                            answ = 'y'
                        if answ == 'y':
                            listtmp.append(img)
                            ii = ii + 1
                            nn = nn + 1
                            if len(listtmp) == 8:
                                print '### number images selected: ', str(len(listtmp))
                                mflat = ntt.soficalibdef.makeflat(listtmp)
                                listaflat.append(mflat)
                                display_image(mflat, 1, '', '', False)
                                nn = 0
                                listtmp = []
                        elif answ == 'r':
                            listtmp = []
                            ii = 0
                            nn = 0
                        elif answ == 's':
                            break
                        else:
                            print len(images), _filter, ID
                        print '### number images selected: ', str(len(listtmp))
    else:
        listaflat = flats
    listaillum = []
    if illum:
        for _filter in illum:
            for ID in illum[_filter]:
                images = illum[_filter][ID]
                flatfield = searchflat(images[0], listaflat)[0]
                if flatfield:
                    illum_frame = ntt.soficalibdef.makeillumination(
                        images, flatfield)
                    listaillum.append(illum_frame)
                else:
                    print 'flat field not found'

        for img in listaflat:
            try:
                ntt.util.phase3header(img)  # phase 3 definitions
                ntt.util.updateheader(
                    img, 0, {'BUNIT': ['ADU', 'pixel units(ADU,electrons)']})
                ntt.util.updateheader(
                    img, 0, {'FILETYPE': [31202, 'flat field']})
            except:
                print '\n### problems with phase 3 definitions'

        for img in listaillum:
            try:
                ntt.util.phase3header(img)  # phase 3 definitions
                ntt.util.updateheader(
                    img, 0, {'BUNIT': ['ADU', 'pixel units(ADU,electrons)']})
                ntt.util.updateheader(
                    img, 0, {'FILETYPE': [31213, 'illum corr  frames']})
            except:
                print '\n### problems with phase 3 definitions'
    return listaflat, listaillum

###########################################################################
