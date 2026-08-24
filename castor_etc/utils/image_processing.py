

import numpy as np

def addflux2pix(px, py, pixels, fmod):
    """Usage: pixels=addflux2pix(px,py,pixels,fmod)

    Drizel Flux onto Pixels using a square PSF of pixel size unity
    px,py are the pixel position (integers)
    fmod is the flux calculated for (px,py) pixel
        and it has the same length as px and py
    pixels is the image.
    """
    # As far as I can tell, this file was built as a Python-ic version of this Fortran function
    # https://github.com/jasonfrowe/JWSTNIRISS/blob/master/ftools/specgen/utils/addflux2pix.f90

    xmax = pixels.shape[0]  # Size of pixel array
    ymax = pixels.shape[1]

    pxmh = px - 0.5  # location of reference corner of PSF square
    pymh = py - 0.5

    dx = np.floor(px + 0.5) - pxmh
    dy = np.floor(py + 0.5) - pymh

    # Supposing right-left as x axis and up-down as y axis:
    # Lower left pixel
    npx = int(pxmh)  # Numpy arrays start at zero
    npy = int(pymh)

    # print('n',npx,npy)

    # if (npx >= 0) & (npx < xmax) & (npy >= 0) & (npy < ymax) :
    #    pixels[npx,npy]=pixels[npx,npy]+fmod

    if (npx >= 0) & (npx < xmax) & (npy >= 0) & (npy < ymax):
        pixels[npx, npy] = pixels[npx, npy] + fmod * dx * dy

    # Same operations are done for the 3 pixels other neighbouring pixels

    # Lower right pixel
    npx = int(pxmh) + 1  # Numpy arrays start at zero
    npy = int(pymh)
    if (npx >= 0) & (npx < xmax) & (npy >= 0) & (npy < ymax):
        pixels[npx, npy] = pixels[npx, npy] + fmod * (1.0 - dx) * dy

    # Upper left pixel
    npx = int(pxmh)  # Numpy arrays start at zero
    npy = int(pymh) + 1
    if (npx >= 0) & (npx < xmax) & (npy >= 0) & (npy < ymax):
        pixels[npx, npy] = pixels[npx, npy] + fmod * dx * (1.0 - dy)

    # Upper right pixel
    npx = int(pxmh) + 1  # Numpy arrays start at zero
    npy = int(pymh) + 1
    if (npx >= 0) & (npx < xmax) & (npy >= 0) & (npy < ymax):
        pixels[npx, npy] = pixels[npx, npy] + fmod * (1.0 - dx) * (1.0 - dy)

    return pixels


def gen_unconv_image(pars, starmodel_flux, xcoo, ycoo):
    # This generates un-convolved images
    xpad = pars.xpad * pars.noversample
    ypad = pars.ypad * pars.noversample
    # array to hold synthetic image
    xmax = pars.xout * pars.noversample + xpad * 2
    ymax = pars.yout * pars.noversample + ypad * 2

    pixels = np.zeros((xmax, ymax))

    i = (xcoo + (pars.xout - pars.ccd_dim[0]) / 2) * pars.noversample
    j = (ycoo + (pars.yout - pars.ccd_dim[1]) / 2) * pars.noversample

    pixels = addflux2pix(i, j, pixels, starmodel_flux)

    return pixels