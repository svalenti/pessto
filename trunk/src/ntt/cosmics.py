#  cosmic correction
#  lacosmic  iraf modules rewritten in pyraf
#
#

import os
import numpy as np
import math

try:       from astropy.io import fits as pyfits
except:    import pyfits

# We define the laplacian kernel to be used
laplkernel = np.array([[0.0, -1.0, 0.0], [-1.0, 4.0, -1.0], [0.0, -1.0, 0.0]])
# Other kernels :
growkernel = np.ones((3, 3))

# FITS import - export


def fromfits(infilename, hdu=0, verbose=True):
    # print "LOGX:: Entering `fromfits` method/function in %(__file__)s" %
    # globals()
    """
        Reads a FITS file and returns a 2D numpy array of the data.
        Use hdu to specify which HDU you want (default = primary = 0)
        """

    pixelarray, hdr = pyfits.getdata(infilename, hdu, header=True)
    pixelarray = np.asarray(pixelarray).transpose()

    pixelarrayshape = pixelarray.shape
    if verbose:
        print("FITS import shape : (%i, %i)" % (pixelarrayshape[0], pixelarrayshape[1]))
        print("FITS file BITPIX : %s" % (hdr["BITPIX"]))
        print("Internal array type :", pixelarray.dtype.name)

    return pixelarray, hdr


def tofits(outfilename, pixelarray, hdr=None, verbose=True):
    # print "LOGX:: Entering `tofits` method/function in %(__file__)s" %
    # globals()
    """
        Takes a 2D numpy array and write it into a FITS file.
        If you specify a header (pyfits format, as returned by fromfits()) it will be used for the image.
        You can give me boolean numpy arrays, I will convert them into 8 bit integers.
        """
    pixelarrayshape = pixelarray.shape
    if verbose:
        print("FITS export shape : (%i, %i)" % (pixelarrayshape[0], pixelarrayshape[1]))

    if pixelarray.dtype.name == "bool":
        pixelarray = np.cast["uint8"](pixelarray)

    if os.path.isfile(outfilename):
        os.remove(outfilename)

    if hdr == None:  # then a minimal header will be created
        hdu = pyfits.PrimaryHDU(pixelarray.transpose())
    else:  # this if else is probably not needed but anyway ...
        hdu = pyfits.PrimaryHDU(pixelarray.transpose(), hdr)

    hdu.writeto(outfilename, output_verify='ignore')

    if verbose:
        print("Wrote %s" % outfilename)

###################################################


def lacos(_input0, output='clean.fits', outmask='mask.fits', gain=1.3, readn=9, xorder=9, yorder=9, sigclip=4.5, sigfrac=0.5, objlim=1, verbose=True, interactive=False):
    # print "LOGX:: Entering `lacos` method/function in %(__file__)s" %
    # globals()
    import ntt
    from ntt.util import delete
    import sys
    import re
    import os
    import string
    from pyraf import iraf
    import numpy as np

    oldoutput, galaxy, skymod, med5 = 'oldoutput.fits', 'galaxy.fits', 'skymod.fits', 'med5.fits'
    blk, lapla, med3, med7, sub5, sigima, finalsel = 'blk.fits', 'lapla.fits', 'med3.fits', 'med7.fits', 'sub5.fits', 'sigima.fits', 'finalsel.fits'
    deriv2, noise, sigmap, firstsel, starreject = 'deriv2.fits', 'noise.fits', 'sigmap.fits', 'firstsel.fits', 'starreject.fits'
    inputmask = 'inputmask.fits'
    # set some parameters in standard IRAF tasks
    iraf.convolve.bilinear = 'no'
    iraf.convolve.radsym = 'no'
    # create Laplacian kernel
    # laplkernel = np.array([[0.0, -1.0, 0.0], [-1.0, 4.0, -1.0], [0.0, -1.0, 0.0]])
    f = open('_kernel', 'w')
    f.write('0 -1 0;\n-1 4 -1;\n0 -1 0')
    f.close()
    # create growth kernel
    f = open('_gkernel', 'w')
    f.write('1 1 1;\n1 1 1;\n1 1 1')
    f.close()
    gkernel = np.array([[1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]])
    delete(galaxy)
    delete(skymod)
    delete(oldoutput)

    if not output:
        output = _input0
    else:
        delete(output)
        iraf.imcopy(_input0, output, verbose='no')

    delete('_xxx.fits,_yyy.fits')
    iraf.imcopy(_input0 + '[350:550,*]', '_xxx.fits', verbose='no')
    _input = '_xxx.fits'

    arrayinput, headerinput = ntt.cosmics.fromfits(_input, verbose=False)
    ntt.cosmics.tofits(outmask, np.float32(
        arrayinput - arrayinput), headerinput, verbose=False)

    # subtract object spectra if desired
    iraf.fit1d(_input, galaxy, "fit", axis=2, order=9, func="leg", low=4.,
               high=4., nav=1, inter='no', sample="*", niter=3, grow=0, cursor="")
    iraf.imarith(_input, "-", galaxy, oldoutput)
    # Subtract sky lines
    iraf.fit1d(oldoutput, skymod, "fit", axis=1, order=5, func="leg", low=4., high=4.,
               inter='no', sample="*", nav=1, niter=3, grow=0, cursor="")
    iraf.imarith(oldoutput, "-", skymod, oldoutput)

    arrayoldoutput, headeroldoutput = ntt.cosmics.fromfits(
        oldoutput, verbose=False)
    # add object spectra to sky model
    iraf.imarith(skymod, "+", galaxy, skymod)
    delete(med5)
    # add median of residuals to sky model
    iraf.median(oldoutput, med5, 5, 5, zlor='INDEF',
                zhir='INDEF', verbose='no')
#    m5 = ndimage.filters.median_filter(_inputarray, size=5, mode='mirror')
    iraf.imarith(skymod, "+", med5, med5)
    # take second-order derivative (Laplacian) of input image
    # kernel is convolved with subsampled image, in order to remove negative
    # pattern around high pixels
    delete(blk)
    delete(lapla)
    iraf.blkrep(oldoutput, blk, 2, 2)
    iraf.convolve(blk, lapla, '_kernel')
    iraf.imreplace(lapla, 0, upper=0, lower='INDEF')
    delete(deriv2)
    delete(noise)
    iraf.blkavg(lapla, deriv2, 2, 2, option="average")
    # create noise model
    iraf.imutil.imexpr(expr='sqrt(a*' + str(gain) + '+' + str(readn) +
                       '**2)/' + str(gain), a=med5, output=noise, verbose='no')
    iraf.imreplace(med5, 0.00001, upper=0, lower='INDEF')
    # divide Laplacian by noise model
    delete(sigmap)
    iraf.imutil.imexpr(expr='(a/b)/2', a=deriv2, b=noise,
                       output=sigmap, verbose='no')
    # removal of large structure (bright, extended objects)
    delete(med5)
    iraf.median(sigmap, med5, 5, 5, zlo='INDEF', zhi='INDEF', verbose='no')
    iraf.imarith(sigmap, "-", med5, sigmap)
    # find all candidate cosmic rays
    # this selection includes sharp features such as stars and HII regions

    arraysigmap, headersigmap = ntt.cosmics.fromfits(sigmap, verbose=False)
    arrayf = np.where(arraysigmap < sigclip, 0, arraysigmap)
    arrayf = np.where(arrayf > 0.1, 1, arrayf)
    ntt.cosmics.tofits(firstsel, np.float32(
        arrayf), headersigmap, verbose=False)

    # compare candidate CRs to median filtered image
    # this step rejects bright, compact sources from the initial CR list
    # subtract background and smooth component of objects
    delete(med3)
    iraf.median(oldoutput, med3, 3, 3, zlo='INDEF', zhi='INDEF', verbose='no')
    delete(med7)
    delete('_' + med3)
    iraf.median(med3, med7, 7, 7, zlo='INDEF', zhi='INDEF', verbose='no')
    iraf.imutil.imexpr(expr='(a-b)/c', a=med3, b=med7,
                       c=noise, output='_' + med3, verbose='no')
    iraf.imreplace('_' + med3, 0.01, upper=0.01, lower='INDEF')
    # compare CR flux to object flux
    delete(starreject)
    iraf.imutil.imexpr(expr='a+b+c', a=firstsel, b=sigmap,
                       c="_" + med3, output=starreject, verbose='no')
    # discard if CR flux <= objlim * object flux
    iraf.imreplace(starreject, 0, upper=objlim, lower='INDEF')
    iraf.imreplace(starreject, 1, lower=0.5, upper='INDEF')
    iraf.imarith(firstsel, "*", starreject, firstsel)

    # grow CRs by one pixel and check in original sigma map
    arrayfirst, headerfirst = ntt.cosmics.fromfits(firstsel, verbose=False)
    arraygfirst = ntt.cosmics.my_convolve_with_FFT2(arrayfirst, gkernel)

    arraygfirst = np.where(arraygfirst > 0.5, 1, arraygfirst)
    arraygfirst = arraygfirst * arraysigmap
    arraygfirst = np.where(arraygfirst < sigclip, 0, arraygfirst)
    arraygfirst = np.where(arraygfirst > 0.1, 1, arraygfirst)

    # grow CRs by one pixel and lower detection limit
    sigcliplow = sigfrac * sigclip
    # Finding neighbouring pixels affected by cosmic rays
    arrayfinal = ntt.cosmics.my_convolve_with_FFT2(arraygfirst, gkernel)
    arrayfinal = np.where(arrayfinal > 0.5, 1, arrayfinal)
    arrayfinal = arrayfinal * arraysigmap
    arrayfinal = np.where(arrayfinal < sigcliplow, 0, arrayfinal)
    arrayfinal = np.where(arrayfinal > 0.1, 1, arrayfinal)

    # determine number of CRs found in this iteration
    arraygfirst = (1 - (arrayfinal - arrayfinal)) * arrayfinal
    npix = [str(int(np.size(np.where(arraygfirst > 0.5)) / 2.))]
    # create cleaned output image; use 3x3 median with CRs excluded
    arrayoutmask = np.where(arrayfinal > 1, 1, arrayfinal)
    ntt.cosmics.tofits(outmask, np.float32(
        arrayoutmask), headerfirst, verbose=False)
    delete(inputmask)
    arrayinputmask = (1 - (10000 * arrayoutmask)) * arrayoldoutput
    ntt.cosmics.tofits(inputmask, np.float32(
        arrayinputmask), headerfirst, verbose=False)
    delete(med5)
    iraf.median(inputmask, med5, 5, 5, zloreject=-
                9999, zhi='INDEF', verbose='no')
    iraf.imarith(outmask, "*", med5, med5)
    delete('_yyy.fits')
    iraf.imutil.imexpr(expr='(1-a)*b+c', a=outmask, b=oldoutput,
                       c=med5, output='_yyy.fits', verbose='no')
    # add sky and object spectra back in
    iraf.imarith('_yyy.fits', "+", skymod, '_yyy.fits')
    # cleanup and get ready for next iteration
    if npix == 0:
        stop = yes
      # delete temp files
    iraf.imcopy('_yyy.fits', output + '[350:550,*]', verbose='no')
    delete(blk + "," + lapla + "," + deriv2 + "," + med5)
    delete(med3 + "," + med7 + "," + noise + "," + sigmap)
    delete(firstsel + "," + starreject)
    delete(finalsel + "," + inputmask)
    delete(oldoutput + "," + skymod + "," + galaxy)
    delete("_" + med3 + ",_" + sigmap)
    delete('_kernel' + "," + '_gkernel')
    delete(outmask)
    delete('_xxx.fits,_yyy.fits')
#################################################


def clean_image(img, cleanimg):
    # print "LOGX:: Entering `clean_image` method/function in %(__file__)s" %
    # globals()
    import ntt
    from ntt.util import readkey3, readhdr, delete
    array, header = ntt.cosmics.fromfits(img, verbose=False)
    import warnings

    def fxn():
        # print "LOGX:: Entering `fxn` method/function in %(__file__)s" %
        # globals()
        warnings.warn(" ", DeprecationWarning)

    original_filters = warnings.filters[:]
    # Ignore warnings.
    warnings.simplefilter("ignore")
    try:
        c = ntt.cosmics.cosmicsimage(array, gain=readkey3(header, 'gain'), readnoise=readkey3(
            header, 'ron'), sigclip=5.0, sigfrac=0.3, objlim=5.0, verbose=False)
        c.run(maxiter=4, verbose=False)
        fxn()
    finally:
        warnings.filters = original_filters

    if not cleanimg:
        delete(img)
        cleanimg = img
    ntt.cosmics.tofits(cleanimg, c.cleanarray, header, verbose=False)
    return cleanimg

##################################################


def my_convolve_with_FFT2(image1, kernel):
    # print "LOGX:: Entering `my_convolve_with_FFT2` method/function in
    # %(__file__)s" % globals()
    import numpy as np
    r1, c1 = image1.shape
    r2, c2 = kernel.shape
    r = r1 + r2 - 1
    c = c1 + c2 - 1
    padr = (r - r1) // 2  # add to each side
    padc = (c - c1) // 2  # add to each of top and bottom

    # pad the edges, with the same values found on the edges
    lside = np.empty((image1.shape[0], padc))
    rside = np.empty((image1.shape[0], padc))
    for i in range(padc):
        lside[:, i] = np.copy(image1[:, 0]).T
        rside[:, i] = np.copy(image1[:, image1.shape[1] - 1]).T
    # end for loop

    image1 = np.hstack((lside, image1, rside))
    top = np.empty((padr, image1.shape[1]))
    bot = np.empty((padr, image1.shape[1]))
    # pad the top and bottom
    for i in range(padr):
        top[i] = np.copy(image1[0, :]).T
        bot[i] = np.copy(image1[image1.shape[0] - 1, :]).T
    # end for loop
    image1 = np.vstack((top, image1, bot))
    rOrig = r
    cOrig = c
    pr2 = int(np.log(r) / np.log(2.0) + 1.0)
    pc2 = int(np.log(c) / np.log(2.0) + 1.0)
    r = 2**pr2
    c = 2**pc2
    fftimage = np.fft.fft2(image1, s=(r, c)) * \
        np.fft.fft2(kernel[::-1, ::-1], s=(r, c))
    ret = np.fft.ifft2(fftimage).real
    return ret[(rOrig - r1):rOrig, (cOrig - c1):cOrig]

###################################################


def lacos_im(_input, _output='clean.fits', outmask='mask.fits', gain=1.3, readn=9, xorder=9, yorder=9, sigclip=4.5, sigfrac=0.5, objlim=1, skyval=0, niter=2, verbose=True, interactive=False):
    # print "LOGX:: Entering `lacos_im` method/function in %(__file__)s" %
    # globals()
    import ntt
    from ntt.util import delete
    import sys
    import re
    import os
    import string
    from pyraf import iraf
    import numpy as np
    iraf.convolve.bilinear = 'no'
    iraf.convolve.radsym = 'no'
    # make temporary files
    oldoutput, galaxy, skymod, med5 = 'oldoutput.fits', 'galaxy.fits', 'skymod.fits', 'med5.fits'
    blk, lapla, med3, med7, sub5, sigima, finalsel = 'blk.fits', 'lapla.fits', 'med3.fits', 'med7.fits', 'sub5.fits', 'sigima.fits', 'finalsel.fits'
    deriv2, noise, sigmap, firstsel, starreject = 'deriv2.fits', 'noise.fits', 'sigmap.fits', 'firstsel.fits', 'starreject.fits'
    inputmask, gfirstsel = 'inputmask.fits', 'gfirstsel.fits'
    f = open('_kernel', 'w')
    f.write('0 -1 0;\n-1 4 -1;\n0 -1 0')
    f.close()
    # create growth kernel
    f = open('_gkernel', 'w')
    f.write('1 1 1;\n1 1 1;\n1 1 1')
    f.close()
    gkernel = np.array([[1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]])
    # initialize loop
    usegain = gain
    i = 1
    stop = 'no'
    previous = 0
    if not _output:
        _output = _input

    arrayinput, headerinput = ntt.cosmics.fromfits(_input, verbose=False)
    ntt.cosmics.tofits(outmask, np.float32(
        arrayinput - arrayinput), headerinput, verbose=False)

    delete(oldoutput)
    if skyval > 0:
        arrayoldoutput = arrayinput + skyval
    else:
        arrayoldoutput = arrayinput
    ntt.cosmics.tofits(oldoutput, np.float32(
        arrayoldoutput), headerinput, verbose=False)
    # start iterations
    while stop == 'no':
        # take second-order derivative (Laplacian) of input image
        # kernel is convolved with subsampled image, in order to remove negative
        # pattern around high pixels
        delete(blk)
        delete(lapla)
        delete(deriv2)
        iraf.blkrep(oldoutput, blk, 2, 2)
        iraf.convolve(blk, lapla, '_kernel')
        iraf.imreplace(lapla, 0, upper=0, lower='INDEF', radius=0)
        iraf.blkavg(lapla, deriv2, 2, 2, option="average")
        delete(med5)
        # create model of background flux - 5x5 box should exclude all CRs
        iraf.median(oldoutput, med5, 5, 5, zlo='INDEF',
                    zhi='INDEF', verbose='no')
        iraf.imreplace(med5, 0.0001, upper=0, lower='INDEF', radius=0)
        # create noise model
        delete(noise)
        iraf.imutil.imexpr(expr='sqrt(a*' + str(usegain) + '+' + str(readn) +
                           '**2)/' + str(usegain), a=med5, output=noise, verbose='no')
        # divide Laplacian by noise model
        delete(sigmap)
        iraf.imarith(deriv2, "/", noise, sigmap)
        # Laplacian of blkreplicated image counts edges twice:
        iraf.imarith(sigmap, "/", 2., sigmap)
        # removal of large structure (bright, extended objects)
        delete(med5)
        iraf.median(sigmap, med5, 5, 5, zlo='INDEF', zhi='INDEF', verbose='no')
        arraysigmap, headersigmap = ntt.cosmics.fromfits(sigmap, verbose=False)
        arraymed5, headermed5 = ntt.cosmics.fromfits(med5, verbose=False)
        arraysigmap = arraysigmap - arraymed5
        iraf.imarith(sigmap, "-", med5, sigmap)
        # find all candidate cosmic rays
        # this selection includes sharp features such as stars and HII regions

        delete(firstsel)
        iraf.imcopy(sigmap, firstsel, verbose='no')
        iraf.imreplace(firstsel, 0, upper=sigclip, lower='INDEF', radius=0)
        iraf.imreplace(firstsel, 1, lower=0.1, upper='INDEF', radius=0)
#		arraygfirst=arraysigmap
#		arraygfirst = np.where(arraygfirst<sigclip,0,arraygfirst)
#		arraygfirst = np.where(arraygfirst>0.1,1,arraygfirst)

        # compare candidate CRs to median filtered image
        # this step rejects bright, compact sources from the initial CR list
        # subtract background and smooth component of objects
        delete(med3)
        delete(med7)
        iraf.median(oldoutput, med3, 3, 3, zlo='INDEF',
                    zhi='INDEF', verbose='no')
        iraf.median(med3, med7, 7, 7, zlo='INDEF', zhi='INDEF', verbose='no')
        iraf.imarith(med3, "-", med7, med3)
        iraf.imarith(med3, "/", noise, med3)
        iraf.imreplace(med3, 0.01, upper=0.01, lower='INDEF', radius=0)
        # compare CR flux to object flux
        delete(starreject)
        iraf.imutil.imexpr(expr="(a*b)/c", a=firstsel, b=sigmap,
                           c=med3, output=starreject, verbose='no')
        # discard if CR flux <= objlim * object flux
        iraf.imreplace(starreject, 0, upper=objlim, lower='INDEF', radius=0)
        iraf.imreplace(starreject, 1, lower=0.5, upper='INDEF', radius=0)
        iraf.imarith(firstsel, "*", starreject, firstsel)
        # grow CRs by one pixel and check in original sigma map
        delete(gfirstsel)
        iraf.convolve(firstsel, gfirstsel, '_gkernel')
        iraf.imreplace(gfirstsel, 1, lower=0.5, upper='INDEF', radius=0)
        iraf.imarith(gfirstsel, "*", sigmap, gfirstsel)
        iraf.imreplace(gfirstsel, 0, upper=sigclip, lower='INDEF', radius=0)
        iraf.imreplace(gfirstsel, 1, lower=0.1, upper='INDEF', radius=0)
        # grow CRs by one pixel and lower detection limit
        sigcliplow = sigfrac * sigclip
        delete(finalsel)
        iraf.convolve(gfirstsel, finalsel, '_gkernel')
        iraf.imreplace(finalsel, 1, lower=0.5, upper='INDEF', radius=0)
        iraf.imarith(finalsel, "*", sigmap, finalsel)
        iraf.imreplace(finalsel, 0, upper=sigcliplow, lower='INDEF', radius=0)
        iraf.imreplace(finalsel, 1, lower=0.1, upper='INDEF', radius=0)
        # determine number of CRs found in this iteration
        delete(gfirstsel)
        iraf.imutil.imexpr(expr="(1-b)*a", a=outmask,
                           b=finalsel, output=gfirstsel, verbose='no')

        npix = iraf.imstat(gfirstsel, fields="npix",
                           lower=0.5, upper='INDEF', Stdout=1)
        # create cleaned output image; use 3x3 median with CRs excluded
        delete(med5)
        iraf.imarith(outmask, "+", finalsel, outmask)
        iraf.imreplace(outmask, 1, lower=1, upper='INDEF', radius=0)

        delete(inputmask)
        iraf.imutil.imexpr(expr="(1-10000*a)", a=outmask,
                           output=inputmask, verbose='no')
        iraf.imarith(oldoutput, "*", inputmask, inputmask)
        delete(med5)
        iraf.median(inputmask, med5, 5, 5, zloreject=-
                    9999, zhi='INDEF', verbose='no')
        iraf.imarith(outmask, "*", med5, med5)
        if i > 1:
            delete(_output)

        delete(_output)
        iraf.imutil.imexpr(expr="(1.-b)*a+c", a=oldoutput,
                           b=outmask, c=med5, output=_output, verbose='no')

        # cleanup and get ready for next iteration
        delete(oldoutput)
        iraf.imcopy(_output, oldoutput, verbose='no')

        if npix == 0:
            stop = 'yes'
        i = i + 1
        if i > niter:
            stop = 'yes'
        # delete temp files
        delete(blk + "," + lapla + "," + deriv2 + "," + med5)
        delete(med3 + "," + med7 + "," + noise + "," + sigmap)
        delete(firstsel + "," + starreject + "," + gfirstsel)
        delete(finalsel + "," + inputmask)

    if skyval > 0:
        iraf.imarith(_output, "-", skyval, _output)
    delete('_kernel' + "," + '_gkernel')
    delete(oldoutput)
#	delete(kernel+","+gkernel)
#################################################################
