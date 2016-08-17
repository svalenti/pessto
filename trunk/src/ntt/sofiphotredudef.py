def pesstocombine(imglist, _combine, _rejection, outputimage):
    # print "LOGX:: Entering `pesstocombine` method/function in %(__file__)s"
    # % globals()
    import ntt
    from pyraf import iraf
    from numpy import compress, array, round, median, std, isnan, sqrt, argmin, argsort
    import string
    import os
    import sys
    import re

    max_sep = 20
    _ron = ntt.util.readkey3(ntt.util.readhdr(imglist[0]), 'ron')
    _gain = ntt.util.readkey3(ntt.util.readhdr(imglist[0]), 'gain')
    namesex = ntt.util.defsex('default.sex')
    os.system('sex ' + imglist[0] + ' -c ' + namesex + ' > _logsex')
    ntt.util.delete(namesex)
    ntt.util.delete('_logsex')
    a = iraf.proto.fields('detections.cat', fields='2,3', Stdout=1)
    a = compress(array(a) != '', array(a))
    if len(a) >= 50:
        am = array(iraf.proto.fields(
            'detections.cat', fields='4', Stdout=1), float)
        a = a[argsort(am)][0:50]
    ff = open('detection0.txt', 'w')
    for i in a:
        ff.write(i + '\n')
    ff.close()
    x, y = [], []
    for i in a:
        x.append(float(string.split(i)[0]))
        y.append(float(string.split(i)[1]))
    x, y = array(x), array(y)
    if ntt.util.readkey3(ntt.util.readhdr(imglist[0]), 'instrume') == 'sofi':
        ref_cumoffsetx0 = ntt.util.readkey3(
            ntt.util.readhdr(imglist[0]), 'xcum')
        ref_cumoffsety0 = ntt.util.readkey3(
            ntt.util.readhdr(imglist[0]), 'ycum')
    else:
        ref_cumoffsetx0 = 0
        ref_cumoffsety0 = 0
    goodlist, _offset = [], []
    for img in imglist:
        if img:
            hdr = ntt.util.readhdr(img)
            if ntt.util.readkey3(hdr, 'instrume') == 'sofi':
                cumoffsetx = ntt.util.readkey3(hdr, 'xcum')
                cumoffsety = ntt.util.readkey3(hdr, 'ycum')
            else:
                cumoffsetx = 0
                cumoffsety = 0
            x1 = x + (float(cumoffsetx) - float(ref_cumoffsetx0))
            y1 = y + (float(cumoffsety) - float(ref_cumoffsety0))
            namesex = ntt.util.defsex('default.sex')
            os.system('sex ' + img + ' -c ' + namesex + ' > _logsex')
            ntt.util.delete(namesex)
            b = iraf.proto.fields('detections.cat', fields='2,3', Stdout=1)
            b = compress(array(b) != '', array(b))
            if len(b) >= 50:
                bm = array(iraf.proto.fields(
                    'detections.cat', fields='4', Stdout=1), float)
                b = b[argsort(bm)][0:50]
            ff = open('detection1.txt', 'w')
            for i in b:
                ff.write(i + '\n')
            ff.close()
            x2, y2 = [], []
            for i in b:
                x2.append(float(string.split(i)[0]))
                y2.append(float(string.split(i)[1]))
            xdist, ydist, iref = [], [], []
            for i in range(len(x2)):
                dist = sqrt((x2[i] - x1) ** 2 + (y2[i] - y1) ** 2)
                idist = argmin(dist)
                if dist[idist] < max_sep:
                    xdist.append(x2[i] - x1[idist])
                    ydist.append(y2[i] - y1[idist])
                    iref.append(idist)

            xoff, xstd = round(median(xdist), 2), round(std(xdist), 2)
            yoff, ystd = round(median(ydist), 2), round(std(ydist), 2)
            _xdist, _ydist = array(xdist), array(ydist)
            __xdist = compress((abs(_xdist - xoff) < 3 * xstd)
                               & (abs(_ydist - yoff) < 3 * ystd), _xdist)
            __ydist = compress((abs(_xdist - xoff) < 3 * xstd)
                               & (abs(_ydist - yoff) < 3 * ystd), _ydist)
            xoff, xstd = round(median(__xdist), 2), round(std(__xdist), 2)
            yoff, ystd = round(median(__ydist), 2), round(std(__ydist), 2)
            if isnan(xoff):
                xoff = 0
            if isnan(yoff):
                yoff = 0
            x3 = x + (float(cumoffsetx) - float(ref_cumoffsetx0)) + float(xoff)
            y3 = y + (float(cumoffsety) - float(ref_cumoffsety0)) + float(yoff)
            if len(xdist) <= 5:
                print '\n### Warning: less than 5 stars, check  ? !!! \n '
            print '### measure offsets between frames ', str(xoff), str(yoff), str(
                round(cumoffsetx + float(xoff), 3)), str(round(cumoffsety + float(yoff), 3))
            if ntt.util.readkey3(hdr, 'instrume') == 'sofi':
                ntt.util.updateheader(img, 0, {'HIERARCH ESO SEQ CUMOFFSETX': [
                                      round(cumoffsetx + float(xoff), 3), '']})
                ntt.util.updateheader(img, 0, {'HIERARCH ESO SEQ CUMOFFSETY': [
                                      round(cumoffsety + float(yoff), 3), '']})
            goodlist.append(img)
            xoffnew = (float(cumoffsetx) -
                       float(ref_cumoffsetx0)) + float(xoff)
            yoffnew = (float(cumoffsety) -
                       float(ref_cumoffsety0)) + float(yoff)
            _offset.append(str((-1) * (xoffnew)) + '   ' +
                           str((-1) * (yoffnew)) + '\n')
    f = open('_offset', 'w')
    for i in _offset:
        f.write(i)
    f.close()
    f = open('_goodlist', 'w')
    for i in goodlist:
        f.write(i + '\n')
    f.close()
    ntt.util.delete('imcombinelog')
    ntt.util.delete(outputimage)
    iraf.images.immatch.imcombine('@_goodlist', outputimage, combine=_combine, reject=_rejection, offset='_offset',
                                  masktyp='', rdnoise=_ron, gain=_gain, zero='mode', logfile='imcombinelog')
    ntt.util.updateheader(
        outputimage, 0, {'M_EPOCH': [False, 'TRUE if resulting from multiple epochs']})
    ntt.util.updateheader(outputimage, 0, {'SINGLEXP': [
                          False, 'TRUE if resulting from single exposure']})
    ntt.util.updateheader(outputimage, 0, {'TEXPTIME': [ntt.util.readkey3(hdr, 'TEXPTIME') * len(goodlist),
                                                        'Total integration time of all exposures (s)']})
    ntt.util.updateheader(outputimage, 0, {'EXPTIME': [ntt.util.readkey3(hdr, 'EXPTIME') * len(goodlist),
                                                       'Total integration time of all exposures (s)']})
    ntt.util.updateheader(outputimage, 0, {'GAIN': [(2. / 3.) * len(goodlist) * ntt.util.readkey3(hdr, 'gain'),
                                                    'Total integration time of all exposures (s)']})

    num = 0
    for im in goodlist:
        num = num + 1
        ntt.util.updateheader(outputimage, 0, {'PROV' + str(num): [ntt.util.readkey3(ntt.util.readhdr(im), 'ARCFILE'),
                                                                   'Originating file']})
        # replace with im if we submit single images
        ntt.util.updateheader(
            outputimage, 0, {'TRACE' + str(num): [im, 'Originating file']})
    hdr = ntt.util.readhdr(outputimage)
    matching = [s for s in hdr.keys() if "IMCMB" in s]
    for imcmb in matching:
        aaa = iraf.hedit(outputimage, imcmb, delete='yes',
                         update='yes', verify='no', Stdout=1)
    if 'SKYSUB' in hdr.keys():
        aaa = iraf.hedit(outputimage, 'SKYSUB', delete='yes',
                         update='yes', verify='no', Stdout=1)

    ntt.util.delete(
        '_offset,_goodlist,detection*.txt,detections.cat,detection_sex*')
    return outputimage


# #######################################################################################

def registersofi(imglist):
    # print "LOGX:: Entering `registersofi` method/function in %(__file__)s" %
    # globals()
    import ntt
    from ntt.util import readhdr, readkey3
    from pyraf import iraf
    from numpy import compress, array, round, median, std, isnan, sqrt, argmin, sort, argsort
    import string
    import os
    import sys
    import re

    max_sep = 20
    hdr0 = ntt.util.readhdr(imglist[0])
    namesex = ntt.util.defsex('default.sex')
    os.system('sex ' + imglist[0] + ' -c ' + namesex + ' > _logsex')
    ntt.util.delete(namesex)
    ntt.util.delete('_logsex')
    a = iraf.proto.fields('detections.cat', fields='2,3', Stdout=1)
    a = compress(array(a) != '', array(a))
    if len(a) >= 50:
        am = array(iraf.proto.fields(
            'detections.cat', fields='4', Stdout=1), float)
        a = a[argsort(am)][0:50]
    ff = open('detection0.txt', 'w')
    for i in a:
        ff.write(i + '\n')
    ff.close()
    x, y = [], []
    for i in a:
        x.append(float(string.split(i)[0]))
        y.append(float(string.split(i)[1]))
    x, y = array(x), array(y)

    if ntt.util.readkey3(hdr0, 'instrume') == 'sofi':
        ref_cumoffsetx0 = readkey3(hdr0, 'xcum')
        ref_cumoffsety0 = readkey3(hdr0, 'ycum')
    else:
        ref_cumoffsetx0 = 0
        ref_cumoffsety0 = 0
    goodlist, _offset = [], []
    for img in imglist:
        if img:
            hdr = readhdr(img)
            if readkey3(hdr, 'instrume') == 'sofi':
                cumoffsetx = readkey3(hdr, 'xcum')
                cumoffsety = readkey3(hdr, 'ycum')
            else:
                cumoffsetx = 0
                cumoffsety = 0
            x1 = x + (float(cumoffsetx) - float(ref_cumoffsetx0))
            y1 = y + (float(cumoffsety) - float(ref_cumoffsety0))
            namesex = ntt.util.defsex('default.sex')
            os.system('sex ' + img + ' -c ' + namesex + ' > _logsex')
            ntt.util.delete(namesex)
            b = iraf.proto.fields('detections.cat', fields='2,3', Stdout=1)
            b = compress(array(b) != '', array(b))
            if len(b) >= 50:
                bm = array(iraf.proto.fields(
                    'detections.cat', fields='4', Stdout=1), float)
                b = b[argsort(bm)][0:50]
            ff = open('detection1.txt', 'w')
            for i in b:
                ff.write(i + '\n')
            ff.close()
            x2, y2 = [], []
            for i in b:
                x2.append(float(string.split(i)[0]))
                y2.append(float(string.split(i)[1]))
            xdist, ydist, iref = [], [], []
            for i in range(len(x2)):
                dist = sqrt((x2[i] - x1) ** 2 + (y2[i] - y1) ** 2)
                idist = argmin(dist)
                if dist[idist] < max_sep:
                    xdist.append(x2[i] - x1[idist])
                    ydist.append(y2[i] - y1[idist])
                    iref.append(idist)

            xoff, xstd = round(median(xdist), 2), round(std(xdist), 2)
            yoff, ystd = round(median(ydist), 2), round(std(ydist), 2)
            _xdist, _ydist = array(xdist), array(ydist)
            __xdist = compress((abs(_xdist - xoff) < 3 * xstd)
                               & (abs(_ydist - yoff) < 3 * ystd), _xdist)
            __ydist = compress((abs(_xdist - xoff) < 3 * xstd)
                               & (abs(_ydist - yoff) < 3 * ystd), _ydist)
            xoff, xstd = round(median(__xdist), 2), round(std(__xdist), 2)
            yoff, ystd = round(median(__ydist), 2), round(std(__ydist), 2)
            if isnan(xoff):
                xoff = 0
            if isnan(yoff):
                yoff = 0

            x3 = x + (float(cumoffsetx) - float(ref_cumoffsetx0)) + float(xoff)
            y3 = y + (float(cumoffsety) - float(ref_cumoffsety0)) + float(yoff)
            if len(xdist) >= 5:
                print '### measure offsets between frames ', str(xoff), str(yoff), str(
                    round(cumoffsetx + float(xoff), 3)), str(round(cumoffsety + float(yoff), 3))
                ntt.util.updateheader(img, 0,
                                      {'CRPIX1': [readkey3(hdr, 'CRPIX1') + round(cumoffsetx + float(xoff), 3), '']})
                ntt.util.updateheader(img, 0,
                                      {'CRPIX2': [readkey3(hdr, 'CRPIX2') + round(cumoffsety + float(yoff), 3), '']})
            else:
                ntt.util.updateheader(
                    img, 0, {'CRPIX1': [readkey3(hdr, 'CRPIX1') + round(cumoffsetx, 3), '']})
                ntt.util.updateheader(
                    img, 0, {'CRPIX2': [readkey3(hdr, 'CRPIX2') + round(cumoffsety, 3), '']})

            ntt.util.updateheader(
                img, 0, {'CRVAL1': [readkey3(hdr0, 'CRVAL1'), '']})
            ntt.util.updateheader(
                img, 0, {'CRVAL2': [readkey3(hdr0, 'CRVAL2'), '']})


def pesstocombine2(imglist, _combine, outputimage):
    # print "LOGX:: Entering `pesstocombine2` method/function in %(__file__)s"
    # % globals()
    import ntt
    from ntt.util import readhdr, readkey3
    from pyraf import iraf
    from numpy import min, max, argmin, float32
    import pyfits
    import string
    import os
    import sys
    import re

    hdr0 = ntt.util.readhdr(imglist[0])
    _ron = ntt.util.readkey3(hdr0, 'ron')
    _gain = ntt.util.readkey3(hdr0, 'gain')
    inputimg = string.join(imglist, ",")
    nameswarp = ntt.util.defswarp(
        'default.swarp', outputimage, _combine, gain=_gain)
    os.system('swarp  ' + str(inputimg) + ' > _logsex')
    ntt.util.delete(nameswarp)
    ntt.util.delete('_logsex')
    data, hdr = pyfits.getdata(outputimage, 0, header=True)
    ntt.util.delete(outputimage)
    pyfits.writeto(outputimage, float32(data), hdr0)
    listheader = ['NAXIS1', 'NAXIS2', 'CRVAL1',
                  'CRVAL2', 'CRPIX1', 'CRPIX2', 'CD2_2', 'CD1_1']
    hedvec = {}
    for hed in listheader:
        hedvec[hed] = [readkey3(hdr, hed), '']
    hedvec['M_EPOCH'] = [False, 'TRUE if resulting from multiple epochs']
    hedvec['SINGLEXP'] = [False, 'TRUE if resulting from single exposure']
    hedvec['TEXPTIME'] = [readkey3(
        hdr0, 'TEXPTIME') * len(imglist), 'Total integration time of all exposures (s)']
    hedvec['EXPTIME'] = [readkey3(
        hdr0, 'EXPTIME') * len(imglist), 'Total integration time of all exposures (s)']
    hedvec['GAIN'] = [(2. / 3.) * len(imglist) * readkey3(hdr0,
                                                          'gain'), 'Total integration time of all exposures (s)']
    #     hedvec['NCOMBINE']:[len(imglist),'Number of raw science data']
    ntt.util.updateheader(outputimage, 0, hedvec)
    iraf.images(_doprint=0)
    num = 0
    #######
    mjdend = []
    mjdstart = []
    for im in imglist:
        num = num + 1
        ntt.util.updateheader(outputimage, 0,
                              {'PROV' + str(num): [readkey3(readhdr(im), 'ARCFILE'), 'Originating file']})
        # replace with im if we submit single images
        ntt.util.updateheader(
            outputimage, 0, {'TRACE' + str(num): [im, 'Originating file']})
        mjdend.append(readkey3(readhdr(im), 'MJD-END'))
        mjdstart.append(readkey3(readhdr(im), 'MJD-OBS'))
    _telapse = (max(mjdend) - min(mjdstart)) * 60. * 60 * 24.
    _tmid = (max(mjdend) + min(mjdstart)) / 2
#    _tmid = (mjdend+float(readkey3(hdr,'MJD-OBS')))/2
    _dateobs = readkey3(readhdr(imglist[argmin(mjdstart)]), 'DATE-OBS')

    header0 = {'DATE-OBS': [_dateobs, 'Date the observation was started (UTC)'],
               'MJD-OBS': [min(mjdstart), 'Start of observations (days)'],
               'MJD-END': [max(mjdend), 'End of observations (days)'],
               'TELAPSE': [_telapse, 'Total elapsed time [days]'], 'TMID': [_tmid, '[d] MJD mid exposure'],
               'ASSON1': [re.sub('.fits', '.weight.fits', outputimage), 'Name of associated file'],
               'ASSOC1': ['ANCILLARY.WEIGHTMAP', 'Associated weigh map image']}
    ntt.util.updateheader(outputimage, 0, header0)

    if 'SKYSUB' in hdr0.keys():
        aaa = iraf.hedit(outputimage, 'SKYSUB', delete='yes',
                         update='yes', verify='no', Stdout=1)

    data, hdr = pyfits.getdata(outputimage, 0, header=True)
    dataw, hdrw = pyfits.getdata(
        re.sub('.fits', '.weight.fits', outputimage), 0, header=True)
    ntt.util.delete(re.sub('.fits', '.weight.fits', outputimage))
    if 'PRODCATG' in hdr:
        hdr.pop('PRODCATG')
    if 'PC1_2' in hdr:
        hdr.pop('PC1_2')
    if 'PC2_1' in hdr:
        hdr.pop('PC2_1')
    pyfits.writeto(re.sub('.fits', '.weight.fits',
                          outputimage), float32(dataw), hdr)

    ntt.util.delete('detection*.txt,detections.cat,detection_sex*')
    return outputimage, re.sub('.fits', '.weight.fits', outputimage)


##########################################################################

def sortbyJD(lista):
    # print "LOGX:: Entering `sortbyJD` method/function in %(__file__)s" %
    # globals()
    from pyfits import open as popen
    from numpy import array, argsort

    JDlist = []
    for img in lista:
        JDlist.append(popen(img)[0].header.get('MJD-OBS'))
    lista = array(lista)
    JDlist = array(JDlist)
    inds = JDlist.argsort()
    sortedlista = lista[inds]
    return list(sortedlista)


def crosstalk(inname, outname):
    # print "LOGX:: Entering `crosstalk` method/function in %(__file__)s" %
    # globals()
    import ntt
    from pyraf import iraf

    iraf.images(_doprint=0)
    iraf.imutil(_doprint=0)
    iraf.imgeom(_doprint=0)
    iraf.ctio(_doprint=0)
    toforget = ['imgeom.blkavg', 'imutil.imcopy',
                'imutil.imarith', 'ctio.imcreate']
    for t in toforget:
        iraf.unlearn(t)

    ntt.util.delete(
        "temp_sum.fits,temp_alpha.fits,temp_corr_A.fits,temp_corr_B.fits,temp_corr.fits," + outname)
    iraf.imgeom.blkavg(input=inname, output="temp_sum.fits",
                       option="sum", b1=1024, b2=1)
    iraf.imutil.imarith("temp_sum.fits", "*", 1.4e-5,
                        "temp_alpha.fits", verbose='no')
    iraf.imgeom.blkrep("temp_alpha.fits", "temp_corr_A.fits", 1024, 1)
    iraf.ctio.imcreate("temp_corr_B.fits", naxis=2,
                       naxis1=1024, naxis2=1024, header='new')
    #  calculate the second part of the correction
    iraf.imutil.imcopy(
        "temp_corr_A.fits[*,1:512]", "temp_corr_B.fits[*,513:1024]", verbose='no')
    iraf.imutil.imcopy(
        "temp_corr_A.fits[*,513:1024]", "temp_corr_B.fits[*,1:512]", verbose='no')
    iraf.imutil.imarith("temp_corr_A.fits", "+",
                        "temp_corr_B.fits", "temp_corr.fits", verbose='no')
    #  combine the 2 corrections
    iraf.imutil.imarith(inname, "-", "temp_corr.fits",
                        outname, verbose='no')  # apply the correction
    #   delete the temporary images and files
    ntt.util.delete(
        "temp_sum.fits,temp_alpha.fits,temp_corr_A.fits,temp_corr_B.fits,temp_corr.fits")


def searchill(flat, listill):
    # print "LOGX:: Entering `searchill` method/function in %(__file__)s" %
    # globals()
    from ntt.util import readhdr, readkey3
    import ntt
    import glob
    import string

    hdrf = readhdr(flat)
    _instrume = readkey3(hdrf, 'instrume')
    filter0 = readkey3(hdrf, 'filter')
    if not listill:
        directory = ntt.__path__[0] + '/archive/' + \
            str(_instrume) + '/illumination/' + filter0
        listill = glob.glob(directory + '/*fits')
    else:
        directory = ''
    illfile = ''
    if listill:
        for ill in listill:
            if string.split(flat, '/')[-1] in readkey3(readhdr(ill), 'MKILLUM'):
                illfile = ill
    return illfile, directory


def skysub(lista, _ron, _gain, _interactive, regi='crreject'):
    # print "LOGX:: Entering `skysub` method/function in %(__file__)s" %
    # globals()
    import re
    import string
    import os
    import ntt
    from ntt.util import readkey3, readhdr
    from pyraf import iraf
    from pyfits import open as popen
    from numpy import mean

    iraf.nproto(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.ccdred(_doprint=0)
    iraf.images(_doprint=0)
    toforget = ['imutil.imarith', 'ccdred.ccdproc', 'nproto.objmasks', 'ccdred.flatcombine', 'imutil.hedit',
                'immatch.imcombine']
    for t in toforget:
        iraf.unlearn(t)

    iraf.nproto.objmasks1.fitxord = 1
    iraf.nproto.objmasks1.fityord = 1
    iraf.ccdproc.darkcor = 'no'
    iraf.ccdproc.fixpix = 'no'
    iraf.ccdproc.flatcor = 'no'
    iraf.ccdproc.zerocor = 'no'
    iraf.ccdproc.ccdtype = ''
    iraf.ccdproc.overscan = 'no'
    iraf.ccdred.instrument = "/dev/null"

    f = open('tmplist0', 'w')
    #     ll=open('tmplist_fastmask','w')
    gg = open('tmplist_mask', 'w')
    dd = open('tmplist_fastsub', 'w')
    ee = open('tmplist_sky', 'w')
    for im in lista:
        f.write(im + '\n')
        gg.write('mask_' + im + '\n')
        ntt.util.delete('mask_' + im)
        dd.write('fastsub_' + re.sub('.fits', '', im) + '\n')
        ntt.util.delete('fastsub_' + im)
        #          ll.write('fastmask_'+re.sub('.fits','',im)+'\n')
        #          delete('fastmask_'+im)
        ee.write(re.sub('.fits', '_sky.fits', im) + '\n')
        ntt.util.delete(re.sub('.fits', '_sky.fits', im))
    f.close()
    gg.close()
    #     ll.close()
    dd.close()
    ee.close()
    ntt.util.delete('fastsky_' + lista[0])
    ntt.util.delete('sky_' + lista[0])
    listaout = []
    #################################
    #     ccc=iraf.nproto.objmasks(images='@tmplist0', objmasks='@tmplist_fastmask', omtype= 'boolean',\
    #                        blksize=-16, convolv='block 3 3',hsigma=5, lsigma=3, minpix= 10, ngrow= 2, agrow= 4., logfile='',Stdout=1)
    #     for im in lista:  ntt.util.updateheader(im,0,{'OBJMASK':[readkey3(readhdr('fastmask_'+im),'OBJMASK'),'mask']})
    #     iraf.image.immatch.imcombine('@tmplist0',output='fastsky_'+lista[0],masktyp='!OBJMASK',maskval=0,combine='median', reject='none',
    #                                  #scale='mode',statsec='[100:800,100:800]',
    #                                  offsets='', rdnoise=_ron,gain=_gain, nlow=1, nhigh=1, logfile='')
    #     for im in lista:  aaa=iraf.hedit(im,'OBJMASK',delete='yes',update='yes',verify='no',Stdout=1)
    iraf.ccdred.flatcombine('@tmplist0', output='fastsky_' + lista[0], rdnoise=_ron, gain=_gain, ccdtype='',
                            combine='average', reject='avsigclip')
    #####################
    iraf.imarith('@tmplist0', '-', 'fastsky_' +
                 lista[0], result='@tmplist_fastsub', verbose='no')
    ntt.util.delete('objmasklog')
    #  create object mask for each frame
    ccc = iraf.nproto.objmasks(images='@tmplist_fastsub', objmasks='@tmplist_mask', omtype='boolean',
                               blksize=-16, convolv='block 3 3', hsigma=5, lsigma=4, minpix=10, ngrow=2, agrow=4.,
                               # blksize=-16, convolv='block 3 3', hsigma=10,
                               # lsigma=10, minpix=10, ngrow=2, agrow=4.,
                               logfile='', Stdout=1)
    ntt.util.delete('imcombinelog')

    # create sky masking the objects
    for im in lista:
        ntt.util.updateheader(
            im, 0, {'OBJMASK': [readkey3(readhdr('fastsub_' + im), 'OBJMASK'), 'mask']})
    iraf.image.immatch.imcombine('@tmplist0', output='sky_' + lista[0], masktyp='!OBJMASK', maskval=0, combine='median',
                                 reject=regi,
                                 scale='mode', statsec='[100:800,100:800]', offsets='', rdnoise=_ron, gain=_gain,
                                 nlow=1, nhigh=1, logfile='imcombinelog')

    for im in lista:
        aaa = iraf.hedit(im, 'OBJMASK', delete='yes',
                         update='yes', verify='no', Stdout=1)

    hdr = readhdr('sky_' + lista[0])
    matching = [s for s in hdr.keys() if "IMCMB" in s]
    for imcmb in matching:
        aaa = iraf.hedit(
            'sky_' + lista[0], imcmb, delete='yes', update='yes', verify='no', Stdout=1)

    skyfile = ['sky_' + lista[0]]
    ntt.util.updateheader(
        'sky_' + lista[0], 0, {'FILETYPE': [31117, 'sky image']})

    iraf.imutil.imarith('@tmplist0', '-', 'sky_' +
                        lista[0], result='@tmplist_sky', verbose='no')
    num = 0
    for im in lista:
        num = num + 1
        hedvec = {'skysub': ['sky_' + lista[0], 'sky image subtracted'],
                  'FILETYPE': [32215, 'pre-reduced image sky subtracted'],
                  'TRACE1': [im, ''], 'MBKG': [mean(popen('sky_' + lista[0])[0].data), 'background level']}
        ntt.util.updateheader(re.sub('.fits', '_sky.fits', im), 0, hedvec)
        ntt.util.updateheader('sky_' + lista[0], 0,
                              {'PROV' + str(num): [readkey3(readhdr(im), 'ARCFILE'), 'Originating file']})
        ntt.util.updateheader(
            'sky_' + lista[0], 0, {'TRACE' + str(num): [im, 'Originating file']})

    ntt.util.delete('tmplist0,tmplist_s,tmplist_mask,tmplist_sky')
    for img in lista:
        #          delete('fastsky_'+img)
        #          delete('fastsub_'+img)
        #          delete('mask_'+img)                    # deleting masks
        listaout.append(re.sub('.fits', '_sky.fits', img))
    return listaout, skyfile


def skysuboff(listaon, listaoff, _ron, _gain, _interactive, namesky, regi='crreject'):
    # print "LOGX:: Entering `skysuboff` method/function in %(__file__)s" %
    # globals()
    import re
    import string
    import os
    import ntt
    from ntt.util import readkey3, readhdr
    from pyraf import iraf
    from pyfits import open as popen
    from numpy import mean

    iraf.nproto(_doprint=0)
    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.images(_doprint=0)
    iraf.ccdred(_doprint=0)
    toforget = ['imutil.imarith', 'ccdred.ccdproc', 'nproto.objmasks', 'ccdred.flatcombine', 'imutil.hedit',
                'immatch.imcombine']
    for t in toforget:
        iraf.unlearn(t)

    iraf.nproto.objmasks1.fitxord = 1
    iraf.nproto.objmasks1.fityord = 1
    iraf.ccdproc.darkcor = 'no'
    iraf.ccdproc.fixpix = 'no'
    iraf.ccdproc.flatcor = 'no'
    iraf.ccdproc.zerocor = 'no'
    iraf.ccdproc.ccdtype = ''
    iraf.ccdproc.overscan = 'no'
    iraf.ccdred.instrument = "/dev/null"

    listaout = []
    gg = open('tmplist_on', 'w')
    hh = open('tmplist_sky', 'w')
    for im in listaon:
        gg.write(im + '\n')
        hh.write(re.sub('.fits', '_sky.fits', im) + '\n')  # 'sky_'+im+'\n')
        ntt.util.delete(re.sub('.fits', '_sky.fits', im))  # 'sky_'+im)
        listaout.append(re.sub('.fits', '_sky.fits', im))  # 'sky_'+im)
    gg.close()
    hh.close()
    f = open('tmplist_mask', 'w')
    gg = open('tmplist_off', 'w')
    hh = open('tmplist_s', 'w')
    for im in listaoff:
        gg.write(im + '\n')
        f.write('mask_' + im + '\n')
        ntt.util.delete('mask_' + im)
        hh.write('fastsub_' + re.sub('.fits', '', im) + '\n')
        ntt.util.delete('fastsub_' + im)
    f.close()
    gg.close()
    hh.close()
    ntt.util.delete('fastskyoff.fits')

    iraf.flatcombine('@tmplist_off', output='fastskyoff.fits', rdnoise=_ron, gain=_gain, ccdtype='', combine='average',
                     reject='avsigclip')
    iraf.imarith('@tmplist_off', '-', 'fastskyoff.fits',
                 result='@tmplist_s', verbose='no')
    ntt.util.delete('logobjmask')
    ccc = iraf.nproto.objmasks(images='@tmplist_s', objmasks='@tmplist_mask', omtype='boolean',
                               blksize=-16, convolv='block 3 3', hsigma=4, lsigma=3, minpix=10, ngrow=2, agrow=4.,
                               logfile='', Stdout=1)
    ntt.util.delete('imcombinelog')
    ntt.util.delete(namesky)

    for im in listaoff:
        ntt.util.updateheader(im, 0, {'OBJMASK': [ntt.util.readkey3(
            ntt.util.readhdr('fastsub_' + im), 'OBJMASK'), 'mask']})
    iraf.images.immatch.imcombine('@tmplist_off', output=namesky, masktyp='!OBJMASK', maskval=0, combine='median',
                                  reject=regi, nlow=1, nhigh=2,
                                  scale='mode', rdnoise=_ron, gain=_gain, offsets='', logfile='imcombinelog')
    for im in listaoff:
        aaa = iraf.hedit(im, 'OBJMASK', delete='yes',
                         update='yes', verify='no', Stdout=1)

    skyfile = [namesky]
    ntt.util.updateheader(namesky, 0, {'FILETYPE': [31117, 'sky image']})

    hdr = readhdr(namesky)
    matching = [s for s in hdr.keys() if "IMCMB" in s]
    for imcmb in matching:
        aaa = iraf.hedit(namesky, imcmb, delete='yes',
                         update='yes', verify='no', Stdout=1)

    iraf.imarith('@tmplist_on', '-', namesky,
                 result='@tmplist_sky', verbose='no')
    for im in listaon:
        hedvec = {'skysub': [namesky, 'sky image subtracted'], 'FILETYPE': [32215, 'pre-reduced image sky subtracted'],
                  'TRACE1': [im, ''], 'MBKG': [mean(popen(namesky)[0].data), 'background level']}
        ntt.util.updateheader(re.sub('.fits', '_sky.fits', im), 0, hedvec)
    ntt.util.delete(
        'tmplist_on,tmplist_off,tmplist_sky,tmplist_s,tmplist_mask')
    num = 0
    for im in listaoff:
        num = num + 1
        ntt.util.updateheader(namesky, 0, {
                              'PROV' + str(num): [ntt.util.readkey3(ntt.util.readhdr(im), 'ARCFILE'), 'Originating file']})
        ntt.util.updateheader(namesky, 0, {'TRACE' + str(num): [im, '']})
        ntt.util.delete('fastsub_' + im)
        ntt.util.delete('mask_' + im)
    return listaout, skyfile


##########################################################################

def sofireduction(imglist, listill, listflat, _docross, _doflat, _doill, _interactive, _regi='crreject', _verbose=False,
                  method='iraf'):
    # print "LOGX:: entering sofireduction\n"

    from numpy import pi, cos, sin, arccos, array, argmin, min, isnan, sqrt
    from ntt.util import readhdr, readkey3, delete, name_duplicate, updateheader, correctcard
    import time
    import ntt
    import datetime
    import string
    import os
    import re
    import pickle

    scal = pi / 180.
    from pyraf import iraf

    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.ccdred(_doprint=0)

    toforget = ['ccdred.ccdproc', 'imutil.imarith']
    for t in toforget:
        iraf.unlearn(t)

    iraf.ccdproc.darkcor = 'no'
    iraf.ccdproc.fixpix = 'no'
    iraf.ccdproc.flatcor = 'no'
    iraf.ccdproc.zerocor = 'no'
    iraf.ccdproc.ccdtype = ''
    iraf.ccdred.instrument = "/dev/null"

    if _verbose:
        iraf.ccdred.verbose = 'yes'
    else:
        iraf.ccdred.verbose = 'no'
    now = datetime.datetime.now()
    datenow = now.strftime('20%y%m%d%H%M')
    MJDtoday = 55927 + (datetime.date.today() - datetime.date(2012, 01, 01)).days
    imglist = ntt.sofiphotredudef.sortbyJD(imglist)

    objects = {}
    for img in imglist:
        hdr = readhdr(img)
        _object = readkey3(hdr, 'object')
        _filter = readkey3(hdr, 'filter')
        _type = ''
        if _filter.lower() in ['grf', 'gbf']:
            _type = 'spectrum'
        if not _type:
            if _object.lower() == 'flat':
                _type = 'flat'
            elif _object.lower() == 'dark':
                _type = 'dark'
            elif _object.lower() == 'sky':
                _type = 'sky'
        if not _type:
            if str(readkey3(hdr, 'xcum')) != 'None' and str(readkey3(hdr, 'ycum')) != 'None':
                if _filter not in objects:
                    objects[_filter] = [img]
                else:
                    objects[_filter].append(img)
    fieldlist = {}
    OBID = {}
    outputobject = []
    for _setup in objects:
        for obj in objects[_setup]:
            hdro = readhdr(obj)
            _object_name = readkey3(hdro, 'object')
            _ra = readkey3(hdro, 'RA')
            _dec = readkey3(hdro, 'DEC')
            _OBID = (readkey3(hdro, 'esoid'), _setup)
            if string.count(_object_name, '/') or string.count(_object_name, '.') or string.count(_object_name, ' '):
                nameobj = string.split(_object_name, '/')[0]
                nameobj = string.split(nameobj, ' ')[0]
                nameobj = string.split(nameobj, '.')[0]
            else:
                nameobj = _object_name
            if _OBID not in OBID:
                count = 1
                nameobj0 = nameobj + '_' + str(count)
                answ = 'yes'
                while answ == 'yes':
                    if nameobj0 in fieldlist:
                        count = count + 1
                        nameobj0 = nameobj + '_' + str(count)
                    else:
                        answ = 'no'
                fieldlist[nameobj0] = {}
                OBID[readkey3(hdro, 'esoid'), _setup] = nameobj0
            if _setup not in fieldlist[nameobj0]:
                fieldlist[nameobj0][_setup] = []
            fieldlist[nameobj0][_setup].append(obj)
    flat = {}
    if _doflat:
        if listflat:
            for imgflat in listflat:
                _filter = readkey3(readhdr(imgflat), 'filter')
                if _filter not in flat:
                    flat[_filter] = imgflat
        for _set in objects:
            if _set not in flat:
                flat[_set] = ntt.util.searchflat(objects[_set][0], '')[0]

    if _verbose:
        print flat
        print listill

    ill = {}
    if _doill and _doflat:
        for _set in flat:
            if listill:
                for imgill in listill:
                    print flat[_set], _set
                    print imgill
                    if string.split(flat[_set], '/')[-1] in readkey3(readhdr(imgill),
                                                                     'MKILLUM'):  # popen(imgill)[0].header.get('MKILLUM'):
                        ill[_set] = imgill
            else:
                ill[_set] = ntt.sofiphotredudef.searchill(flat[_set], '')[0]

    fieldlist2 = {}
    for _set in objects:
        for field in fieldlist:
            if _set in fieldlist[field]:
                if _interactive:
                    print '\n### next object'
                    for image in fieldlist[field][_set]:
                        try:
                            print '### ', image, str(readkey3(readhdr(image), 'xcum')), str(
                                readkey3(readhdr(image), 'ycum'))
                        except:
                            print '### ', image, ' dither  not define'
                    answ = raw_input('\n### do you want to reduce this object ' + str(field) + ' and filter ' + str(
                        _set) + ' [[y],n] ? ')
                    if not answ:
                        answ = 'y'
                else:
                    answ = 'y'
                if answ in ['YES', 'yes', 'y', 'Y']:
                    img = fieldlist[field][_set][0]
                    if _doill and _set in ill:
                        if str(ill[_set])[0] == '/':
                            _illum = string.split(ill[_set], '/')[-1]
                            ntt.util.delete(_illum)
                            iraf.images.imutil.imcopy(
                                ill[_set], output=_illum, verbose='no')
                        else:
                            _illum = ill[_set]
                    else:
                        _illum = ''
                    if _doflat and flat[_set]:
                        if str(flat[_set])[0] == '/':
                            _masterflat = string.split(flat[_set], '/')[-1]
                            os.system('rm -rf  ' + _masterflat)
                            iraf.images.imutil.imcopy(
                                flat[_set], output=_masterflat, verbose='no')
                        else:
                            _masterflat = flat[_set]
                    else:
                        _masterflat = ''
                    for image in fieldlist[field][_set]:
                        _object_name = readkey3(readhdr(image), 'object')
                        _date = readkey3(readhdr(image), 'date-night')
                        if string.count(_object_name, '/') or string.count(_object_name, '.') or string.count(
                                _object_name, ' '):
                            nameobj = string.split(_object_name, '/')[0]
                            nameobj = string.split(nameobj, ' ')[0]
                            nameobj = string.split(nameobj, '.')[0]
                        else:
                            nameobj = _object_name
                        if field not in fieldlist2:
                            fieldlist2[field] = {}
                        if _set not in fieldlist2[field]:
                            fieldlist2[field][_set] = []
                        nameobj2 = nameobj + '_' + str(_date)
                        nameobj2 = nameobj2 + '_' + \
                            str(_set) + '_' + str(MJDtoday)
                        nameobjnew = name_duplicate(image, nameobj2, '')
                        # print '\n### ',image, nameobjnew,_masterflat,_illum
                        if _docross:
                            ntt.util.delete('C' + nameobjnew)
                            ntt.sofiphotredudef.crosstalk(
                                image, 'C' + nameobjnew)
                            correctcard('C' + nameobjnew)
                            ntt.util.updateheader(
                                'C' + nameobjnew, 0, {'CROSSTAL': ['True', '']})
                            print '\n### image corrected for cross talk   ...... done '
                        else:
                            iraf.images.imutil.imcopy(
                                image, 'C' + nameobjnew, verbose='no')
                            correctcard('C' + nameobjnew)
                            ntt.util.updateheader(
                                'C' + nameobjnew, 0, {'CROSSTAL': ['False', '']})
                        if _doill and _illum:
                            _illumco = 'yes'
                            print '### image corrected for illumination correction   ...... done '
                        else:
                            _illumco = 'no'
                        if _doflat and _masterflat:
                            _flatcor = 'yes'
                            print '### image corrected for flat field   ...... done '
                        else:
                            _flatcor = 'no'
                        ntt.util.delete(nameobjnew)
                        iraf.noao.imred.ccdred.ccdproc('C' + nameobjnew, output=nameobjnew, overscan="no", trim="yes",
                                                       zerocor="no", flatcor=_flatcor,
                                                       illumco=_illumco, trimsec='[1:1024,1:1007]', biassec='',
                                                       flat=_masterflat, illum=_illum, Stdout=1)
                        correctcard(nameobjnew)
                        print '### input= ' + str(image)
                        print '### output= ' + str(nameobjnew)
                        if nameobjnew not in outputobject:
                            outputobject.append(nameobjnew)
                        hdrn = readhdr(nameobjnew)
                        #  changed for DR2   2014-05-18
                        mjdend = float(readkey3(hdrn, 'MJD-OBS')) + (float(readkey3(hdrn, 'ndit')) * float(
                            readkey3(hdrn, 'dit')) + 1.8) / (60. * 60. * 24.)
                        texp = float(readkey3(hdrn, 'dit')) * \
                            float(readkey3(hdrn, 'ndit'))
                        strtexp = time.strftime('%H:%M:%S', time.gmtime(texp))

                        hedvec = {'DIT': [readkey3(hdrn, 'dit'), 'Integration Time'],
                                  'NDIT': [readkey3(hdrn, 'ndit'), 'Number of sub-integrations'],
                                  'FILETYPE': [32104, 'pre-reduced image'],
                                  'M_EPOCH': [False, 'TRUE if resulting from multiple epochs'],
                                  'SINGLEXP': [True, 'TRUE if resulting from single exposure'],
                                  'TEXPTIME': [texp, 'Total integ. time of all exposure (s)'],
                                  'EXPTIME': [texp, 'Total integ. time of all exposure (s) ' + strtexp],
                                  'PROV1': [readkey3(hdrn, 'ARCFILE'), 'Originating file'],
                                  'TRACE1': [readkey3(hdrn, 'ARCFILE'), 'Originating file'],
                                  'MJD-END': [mjdend, 'End of observations (days)']}
                        ntt.util.updateheader(nameobjnew, 0, hedvec)
                        ntt.util.airmass(nameobjnew)  # phase 3 definitions
                        if _doill and _illum:
                            updateheader(nameobjnew, 0, {'ILLUMCOR': [
                                         _illum, 'illumination correction']})
                        if _doflat and _masterflat:
                            updateheader(nameobjnew, 0, {'FLATCOR': [
                                         _masterflat, 'flat correction']})
                        ntt.util.delete('C' + nameobjnew)
                        fieldlist2[field][_set].append(nameobjnew)
    for i in fieldlist2:
        for fil in fieldlist2[i]:
            print '\n### next set of images ' + str(i), str(fil)
            lista = fieldlist2[i][fil]
            lista = ntt.sofiphotredudef.sortbyJD(lista)
            if len(lista) <= 3:
                print '\n### warning: less than 4 images'
                answ = raw_input('Are you sure you want to go on ? [[y]/n]')
                if not answ:
                    answ = 'y'
            else:
                answ = 'y'
            if answ in ['Yes', 'YES', 'yes', 'Y', 'y']:
                hdr0 = readhdr(lista[0])
                _ra0 = readkey3(hdr0, 'RA')
                _dec0 = readkey3(hdr0, 'DEC')
                _ron = readkey3(hdr0, 'ron')
                _gain = readkey3(hdr0, 'gain')
                lista2 = ''
                distance, _ycum, _xcum = [], [], []  # check the offset   in-between images
                for img in lista:
                    hdr = readhdr(img)
                    _ra = readkey3(hdr, 'RA')
                    _dec = readkey3(hdr, 'DEC')
                    _xcum.append(readkey3(hdr, 'xcum'))
                    _ycum.append(readkey3(hdr, 'ycum'))
            else:
                listaout = ''
            if max(max(array(_xcum)), max(array(_ycum))) <= 400 and max(max(array(_xcum)), max(array(_ycum))) >= 1:
                ncombine = 1  # number to count all source raw images that are used
                print '\n### Dithering on source'
                if _interactive:
                    print '### ', str(_xcum), str(_ycum)
                    try:
                        from pylab import plot, ion, show, clf

                        ion()
                        clf()
                        plot(_xcum, _ycum, 'o')
                    except:
                        pass
                    asxx = 'yes'
                    while asxx == 'yes':
                        asxx = 'no'
                        num = raw_input(
                            '\n### How many positions has this mask [' + str(len(_xcum)) + '] ? ')
                        if not num:
                            num = int(len(_xcum))
                        else:
                            if num.isdigit():
                                num = int(num)
                            else:
                                asxx = 'yes'
                                print '\n### Warning: value not valid, try again.....\n'
                else:
                    num = len(lista)
                    if num in [8, 12, 16, 20, 24, 28, 32, 36, 40]:
                        num = 4
                if len(lista) > num and readkey3(readhdr(lista[0]), 'nexp') % num == 0 and len(lista) % num == 0:
                    # popen(lista[0])[0].header.get('nexp') % num == 0 and
                    # len(lista) % num == 0:
                    print '\n### split lista in sample of ' + str(num) + ' images'
                    ii = 0
                    listatmp, listaout = [], []
                    skyfile = []
                    for img in lista:
                        listatmp.append(img)
                        ii = ii + 1
                        if ii == num:
                            listaout0, skyfile0 = ntt.sofiphotredudef.skysub(
                                listatmp, _ron, _gain, _interactive, _regi)
                            listaout = listaout + listaout0
                            skyfile = skyfile + skyfile0
                            ii = 0
                            listatmp = []
                else:
                    listaout, skyfile = ntt.sofiphotredudef.skysub(
                        lista, _ron, _gain, _interactive, _regi)
            elif max(max(array(_xcum)), max(array(_ycum))) <= 1:
                listaout = ''
                print '\n### images without dithering, probably these are aquisition images.'
            else:
                # changed for DR2
                ncombine = 1  # number to count all images that are used ON source
                print '\n## Warning: ON OFF'
                listaon = []
                listaoff = []
                _xcum = []
                _ycum = []
                for img in lista:
                    _xcum.append(readkey3(readhdr(img), 'xcum'))
                    _ycum.append(readkey3(readhdr(img), 'ycum'))
                    print '### ', str(img), str(_xcum[-1]), str(_ycum[-1])
                    if _interactive:
                        _z1, _z2, goon = ntt.util.display_image(
                            img, 1, '', '', False)
                        if abs(int(_ycum[-1])) > 400 or abs(int(_xcum[-1])) > 400:
                            answ0 = 2
                        else:
                            answ0 = 1
                        answ = raw_input(
                            'is this ON[1] or OFF[2] ? [' + str(answ0) + '] ')
                        if not answ:
                            answ = answ0
                        answ = int(answ)
                        if answ == 1:
                            listaon.append(img)
                        elif answ == 2:
                            listaoff.append(img)
                    else:
                        if abs(int(_ycum[-1])) > 400 or abs(int(_xcum[-1])) > 400:
                            listaoff.append(img)
                        else:
                            listaon.append(img)
                print '### image ON \n', listaon, '\n image OFF \n', listaoff
                if _interactive:
                    try:
                        from pylab import plot, ion, show, clf

                        ion()
                        clf()
                        plot(_xcum, _ycum, 'o')
                    #                             show()
                    except:
                        pass
                    asxx = 'yes'
                    while asxx == 'yes':
                        asxx = 'no'
                        num = raw_input(
                            '\n### How many positions has this mask [' + str(len(listaon)) + '] ? ')
                        if not num:
                            num = int(len(listaon))
                        else:
                            if num.isdigit():
                                num = int(num)
                            else:
                                asxx = 'yes'
                                print '\n### Warning: value not valid, try again.....\n'
                else:
                    num = len(listaon)
                hdron = readhdr(listaon[0])
                _object_name = readkey3(hdron, 'object')
                _date = readkey3(hdron, 'date-night')
                _setup = readkey3(hdron, 'filter')
                if string.count(_object_name, '/') or string.count(_object_name, '.') or string.count(_object_name,
                                                                                                      ' '):
                    nameobj = string.split(_object_name, '/')[0]
                    nameobj = string.split(nameobj, ' ')[0]
                    nameobj = string.split(nameobj, '.')[0]
                else:
                    nameobj = _object_name
                nameobj = nameobj + '_' + str(_date)
                nameobj = 'skyoff_' + nameobj + '_' + \
                    str(_setup) + '_' + str(MJDtoday)
                namesky = name_duplicate(listaoff[0], nameobj, '')
                if len(listaon) > num and len(listaoff) > num and readkey3(readhdr(listaon[0]), 'nexp') % num == 0 \
                        and len(listaon) % num == 0 and len(listaoff) % num == 0 and len(listaoff) == len(listaon):
                    #                          popen(lista[0])[0].header.get('nexp') % num == 0 \
                    nn = 0
                    mm = num
                    listaout = []
                    skyfile = []
                    for i in range(0, len(listaon) / num):
                        listaon1 = listaon[nn:mm]
                        listaoff1 = listaoff[nn:mm]
                        nn = nn + num
                        mm = mm + num
                        listaout0, skyfile0 = ntt.sofiphotredudef.skysuboff(listaon1, listaoff1, _ron, _gain,
                                                                            _interactive, namesky, _regi)
                        listaout = listaout + listaout0
                        skyfile = skyfile + skyfile0
                elif len(listaon) == num and len(listaoff) == num:
                    listaout, skyfile = ntt.sofiphotredudef.skysuboff(listaon, listaoff, _ron, _gain, _interactive,
                                                                      namesky, _regi)
                else:
                    print '\n### not enough images to do ON/OFF reduction'
                    print '\n### select manually the images on and off to use '
                    listaout = []
                    skyfile = []
                    goon = raw_input('stop this step ? [n/y]')
                    while goon not in ['yes', 'YES', 'Y', 'y']:
                        print 'liston'
                        kk = 0
                        for gg in listaon:
                            print gg, kk
                            kk = kk + 1
                        kk = 0
                        print 'listoff'
                        for gg in listaoff:
                            print gg, kk
                            kk = kk + 1

                        ddd = raw_input('select list ON [ 0,1,2,3 ] ')
                        listaon1 = []
                        for ii in string.split(ddd, ','):
                            listaon1.append(listaon[int(ii)])
                        print listaon1

                        ddd1 = raw_input('select list OFF [ 0,1,2 ] ')
                        listaoff1 = []
                        for ii in string.split(ddd1, ','):
                            listaoff1.append(listaoff[int(ii)])
                        print listaoff1

                        listaout0, skyfile0 = ntt.sofiphotredudef.skysuboff(listaon1, listaoff1, _ron, _gain,
                                                                            _interactive, namesky, _regi)
                        print listaout0, listaout
                        listaout = listaout + listaout0
                        skyfile = skyfile + skyfile0

                        goon = raw_input('stop this step ? [n/y]')

            if listaout:
                outputobject = outputobject + listaout
                outputobject = outputobject + skyfile
                hdrout = readhdr(listaout[0])
                _object_name = readkey3(hdrout, 'object')
                _filter = readkey3(hdrout, 'filter')
                _date = readkey3(hdrout, 'date-night')
                nameobj = _object_name + '_' + _date + '_' + \
                    _filter + '_merge' + '_' + str(MJDtoday)
                nameobjnew = name_duplicate(listaout[0], nameobj, '')
                _combine = 'median'
                #                   _rejection=_regi
                outputimage = 'merge.fits'
                ntt.util.delete('merge.fits')
                try:
                    ntt.sofiphotredudef.registersofi(listaout)
                    nameobjnew, nameobjnew1 = ntt.sofiphotredudef.pesstocombine2(
                        listaout, _combine, nameobjnew)
                    _xref, _yref = getreferencepixels(
                        listaout[0], nameobjnew)  # added for astrometry
                    hedvec0 = {'NCOMBINE': [len(listaout) * ncombine, 'Number of raw science data'],
                               'NOFFSETS': [int(num), 'Number of offset positions'],
                               'NUSTEP': [0, 'Number of microstep positions'],
                               'NJITTER': [int(float(readkey3(readhdr(nameobjnew), 'nexp')) / num),
                                           'Number of microstep positions'],
                               'NTCRPIX1': [float(_xref), 'reference x pixel of combined image'],
                               'NTCRPIX2': [float(_yref), 'reference y pixel of combined image']}
                    ntt.util.updateheader(nameobjnew1, 0, hedvec0)
                    ntt.util.updateheader(nameobjnew, 0, hedvec0)
                    ntt.util.updateheader(
                        nameobjnew, 0, {'FILETYPE': [32216, 'combined dithered sofi images']})
                    ntt.util.updateheader(
                        nameobjnew1, 0, {'FILETYPE': [31214, 'weight.mask']})
                    ntt.util.phase3header(nameobjnew)  # phase 3 definitions
                    ntt.util.phase3header(nameobjnew1)  # phase 3 definitions
                    if nameobjnew not in outputobject:
                        outputobject.append(nameobjnew)
                    if nameobjnew1 not in outputobject:
                        outputobject.append(nameobjnew1)
                except Exception, e:
                    print e

                try:
                    # print "LOGX:: running sextractor\n"
                    sexvec = ntt.efoscastrodef.sextractor(nameobjnew)
                    rmsx3, rmsy3, num3, fwhmgess, ellgess, ccc, rasys3, decsys3, mbkg3 = ntt.efoscastrodef.efoscastroloop(
                        [nameobjnew], '2mass', False, 40, 40, 20, 'rxyscale', 100, 30, sexvec, True, 10, method)
                    if rmsx3 > 2 or rmsy3 > 2:
                        rmsx3, rmsy3, num3, fwhmgess, ellgess, ccc, rasys3, decsys3, mbkg3 = ntt.efoscastrodef.efoscastroloop(
                            [nameobjnew], '2mass', False, int(20), int(20),
                            int(10), 'rxyscale', 100, 30, sexvec, True, 5, method)
                    astrostring = str(rmsx3) + ' ' + \
                        str(rmsy3) + ' ' + str(num3)
                    ntt.util.updateheader(
                        nameobjnew, 0, {'ASTROMET': [astrostring, 'rmsx rmsy nstars']})
                    print '\n### check astrometry: fine \n### rmsx rmsy nstars: ' + astrostring
                except Exception, e:
                    print e
                    rmsx3, rmsy3, num3, fwhmgess, ellgess, ccc, rasys3, decsys3, mbkg3 = '', '', '', '', '', '', '', '', ''
                    print '\n### problem with astrometry, do you have network ? '
                # print "LOGX:: PSF_FWHM = %(fwhmgess)s\n" % locals()
                if fwhmgess and fwhmgess < 99:
                    # print "LOGX:: PSF_FWHM < 99\n"
                    hedvec = {'PSF_FWHM': [fwhmgess, 'Spatial resolution (arcsec)'],
                              'ELLIPTIC': [ellgess, 'Average ellipticity of point sources'],
                              'CUNIT1': ['deg', 'unit of the coord. trans.'],
                              'CUNIT2': ['deg', 'unit of the coord. trans.'],
                              'CRDER1': [(1 / sqrt(2.)) * float(rmsx3) * (1. / 3600.), 'Random error (degree)'],
                              'CRDER2': [(1 / sqrt(2.)) * float(rmsy3) * (1. / 3600.), 'Random error (degree)'],
                              'CSYER1': [rasys3, 'Systematic error in  (RA_m - Ra_ref)'],
                              'CSYER2': [decsys3, 'Systematic error in (DEC_m - DEC_ref)']}
                    ntt.util.updateheader(nameobjnew, 0, hedvec)
                    result = ntt.efoscastrodef.zeropoint(
                        nameobjnew, '2mass', False, False)
                    if result:
                        if os.path.isfile(re.sub('.fits', '.ph', nameobjnew)):
                            if re.sub('.fits', '.ph', nameobjnew) not in outputobject:
                                outputobject.append(
                                    re.sub('.fits', '.ph', nameobjnew))
                        print '\n### zeropoint ..... done'
                        for ll in result:
                            valore = '%3.3s %6.6s %6.6s' % (
                                str(ll), str(result[ll][1]), str(result[ll][0]))
                            print '### ', valore
                            ntt.util.updateheader(
                                nameobjnew, 0, {'zp' + ll: [str(valore), '']})
                        valore = ''
                        ntt.util.updateheader(nameobjnew, 0,
                                              {'FLUXCAL': ['ABSOLUTE', 'Certifies the validity of PHOTZP']})
                    else:
                        ntt.util.updateheader(nameobjnew, 0,
                                              {'FLUXCAL': ['UNCALIBRATED', 'Certifies the validity of PHOTZP']})
                        ntt.util.updateheader(
                            nameobjnew, 0, {'PHOTSYS': ['NULL', 'Photometric system VEGA or AB']})
                        ntt.util.updateheader(
                            nameobjnew, 0, {'PHOTZP': [9999., 'MAG=-2.5*log(data)+PHOTZP']})
                else:
                    hedvec = {'PSF_FWHM': [9999., 'Spatial resolution (arcsec)'],
                              'ELLIPTIC': [9999., 'Average ellipticity of point sources'],
                              'CUNIT1': ['deg', 'unit of the coord. trans.'],
                              'CUNIT2': ['deg', 'unit of the coord. trans.'],
                              'CRDER1': [9999., 'Random error (degree)'],
                              'CRDER2': [9999., 'Random error (degree)'],
                              'CSYER1': [9999., 'Systematic error in  (RA_m - Ra_ref)'],
                              'CSYER2': [9999., 'Systematic error in (DEC_m - DEC_ref)'],
                              'FLUXCAL': ['UNCALIBRATED', 'Certifies the validity of PHOTZP'],
                              'PHOTSYS': ['NULL', 'Photometric system VEGA or AB'],
                              'PHOTZP': [9999., 'MAG=-2.5*log(data)+PHOTZP']}
                    ntt.util.updateheader(nameobjnew, 0, hedvec)
                if mbkg3:
                    if readkey3(readhdr(nameobjnew), 'FLUXCAL') == 'ABSOLUTE':
                        try:
                            ntt.util.updateheader(nameobjnew, 0, {
                                'ABMAGSAT': [float(mbkg3) + float(readkey3(readhdr(nameobjnew), 'PHOTZP')),
                                             'Saturation limit for point sources (AB mags)']})
                        except:
                            ntt.util.updateheader(nameobjnew, 0, {
                                'ABMAGSAT': [float(readkey3(readhdr(nameobjnew), 'PHOTZP')),
                                             'Saturation limit for point sources (AB mags)']})
                    else:
                        ntt.util.updateheader(nameobjnew, 0, {
                            'ABMAGSAT': [float(mbkg3), 'Saturation limit for point sources (AB mags)']})
                else:
                    ntt.util.updateheader(nameobjnew, 0,
                                          {'ABMAGSAT': [9999., 'Saturation limit for point sources (AB mags)']})
                print nameobjnew
                maglim = ntt.util.limmag(nameobjnew)
                if maglim:
                    ntt.util.updateheader(nameobjnew, 0,
                                          {'ABMAGLIM': [maglim, '5-sigma limiting AB magnitude for point sources']})
                else:
                    ntt.util.updateheader(nameobjnew, 0,
                                          {'ABMAGLIM': [9999., '5-sigma limiting AB magnitude for point sources']})

    reduceddata = ntt.util.rangedata(outputobject)
    print '\n### adding keiwords for phase 3 ....... '
    f = open('logfile_phot_' + str(reduceddata) +
             '_' + str(datenow) + '.raw.list', 'w')
    for img in outputobject:
        print img
        if str(img)[-5:] == '.fits':
            hdr = readhdr(img)
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
            #  added for DR2
            # print img
            # print readkey3(hdr,'NCOMBINE')
            # print readkey3(hdr,'ndit')
            if 'NCOMBINE' in hdr:
                _ncombine = readkey3(hdr, 'NCOMBINE')
            else:
                _ncombine = 1

            ntt.util.updateheader(
                img, 0, {'DETRON ': [12, 'Readout noise per output (e-)']})
            ntt.util.updateheader(img, 0, {'EFFRON': [12. * (1 / sqrt(readkey3(hdr, 'ndit') * _ncombine)) *
                                                      sqrt(pi / 2), 'Effective readout noise per output (e-)']})

            try:
                ntt.util.phase3header(img)  # phase 3 definitions
            except:
                print '\n### problems with phase 3 definitions'
            f.write(readkey3(hdr, 'arcfile') + '\n')
            ntt.util.updateheader(img, 0, {'quality': ['Final', 'final or rapid'],
                                           'BUNIT': ['ADU', 'Physical unit of array values']})
            if readkey3(hdr, 'tech'):
                ntt.util.updateheader(img, 0, {
                    'PRODCATG': ['SCIENCE.' + readkey3(hdr, 'tech').upper(), 'Data product category']})
            if not readkey3(hdr, 'pixscale'):
                ntt.util.updateheader(img, 0,
                                      {'pixscale': [0.288, 'pixel/scale (arcsec)']})
    f.close()
    return outputobject, 'logfile_phot_' + str(reduceddata) + '_' + str(datenow) + '.raw.list'


##########################################################################
def getreferencepixels(img, imgmerge):
    # print "LOGX:: entering getreferencepixels\n"
    # print "LOGX:: Entering `getreferencepixels` method/function in
    # %(__file__)s" % globals()
    from pyraf import iraf

    iraf.imcoords(_doprint=0)
    toforget = ['imcoords.wcsctran']
    for t in toforget:
        iraf.unlearn(t)
    import ntt
    import numpy as np

    ff = open('tmp.pix', 'w')
    ff.write('512  512')
    ff.close()
    ntt.util.delete('tmp.coo,tmp2.pix')
    iraf.wcsctran('tmp.pix', 'tmp.coo', img, inwcs='physical', outwcs='world', columns='1 2', formats='%10.6f %10.6f',
                  verbose='yes')
    iraf.wcsctran('tmp.coo', 'tmp2.pix', imgmerge, inwcs='world', units='degrees degrees', outwcs='logical',
                  columns='1 2', formats='%10.1f %10.1f', verbose='no')
    aa = iraf.proto.fields('tmp2.pix', fields='1', Stdout=1)
    bb = iraf.proto.fields('tmp2.pix', fields='2', Stdout=1)
    aa = np.compress(np.array(aa) != '', np.array(aa))
    bb = np.compress(np.array(bb) != '', np.array(bb))
    ntt.util.delete('tmp.*,tmp2.*')
    return aa[0], bb[0]

##########################################################################
