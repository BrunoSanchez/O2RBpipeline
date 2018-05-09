#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  phot_calibration.py
#
#  Copyright 2018 Bruno S <bruno@oac.unc.edu.ar>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#

import os
import numpy as np

from astropy.io import fits

from astropy import coordinates
from astropy import units as u
from astropy.wcs import WCS

from astropy.stats import sigma_clipped_stats
from astropy.stats import SigmaClip

from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize

from photutils import make_source_mask
from photutils import Background2D
from photutils import MedianBackground
from photutils import SExtractorBackground
from photutils import DAOStarFinder
from photutils import CircularAperture
from photutils import SkyCircularAperture
from photutils import CircularAnnulus
from photutils import SkyCircularAnnulus
from photutils import aperture_photometry

from photutils.psf import BasicPSFPhotometry

from ps_query import panstarrs_query

import matplotlib.pyplot as plt

basepath = '/home/bruno/Data/LIGO_O2/Jan04/newstacks'
image_path = os.path.join(basepath, 'PGC073926/pgc073926_170104.fits')

# load img with WCS
img = image_path.strip('.fits') + '.new'

hdr = fits.getheader(img)
data = fits.getdata(img)
wcs = WCS(hdr)

# Obtain image center
crval_ra= hdr['CRVAL1']
crval_dec = hdr['CRVAL2']
crval_x = hdr['CRPIX1']
crval_y = hdr['CRPIY1']

# query PanSTARRS DR1
ps1_tab = panstarrs_query(crval_ra, crval_dec, 0.5, maxsources=300)

# Get sources in image
mask = make_source_mask(data, snr=2, npixels=5, dilate_size=11)
mask_zeros = data==0
mask = mask | mask_zeros

mean, median, std = sigma_clipped_stats(data, sigma=3.0, mask=mask)
sigma_clip = SigmaClip(sigma=3., iters=10)
bkg_estimator = MedianBackground()
bkg_estimator = SExtractorBackground()

bkg = Background2D(data, (50, 50), filter_size=(3, 3),
                   sigma_clip=sigma_clip,
                   bkg_estimator=bkg_estimator,
                   mask=mask)

back = bkg.background * ~mask

daofind = DAOStarFinder(fwhm=5.0, threshold=8.*std)
sources = daofind(data - median)
print(sources)

positions = (sources['xcentroid'], sources['ycentroid'])
xypos = np.array([sources['xcentroid'], sources['ycentroid']]).T
sky_pos = wcs.all_pix2world(xypos, 0)
skycoords = coordinates.SkyCoord(sky_pos*u.deg, frame='icrs')

radii = [2.*u.arcsec, 3.*u.arcsec, 5.*u.arcsec]
annulus = [(5.*u.arcsec, 7.*u.arcsec),
           (6.*u.arcsec, 8.*u.arcsec),
           (8.*u.arcsec, 10.*u.arcsec)]

apertures = [SkyCircularAperture(skycoords, r=r).to_pixel(wcs)
             for r in radii]
annulii = [SkyCircularAnnulus(skycoords, r_in=ri, r_out=ro).to_pixel(wcs)
           for ri, ro in annulus]
apers = [apertures, annulii]
phot_table = aperture_photometry(data-back, apertures) # did not work annulus



norm = ImageNormalize(stretch=SqrtStretch())
plt.imshow(data - back, origin='lower', cmap='Greys_r', norm=norm)
plt.show()




apertures = CircularAperture(positions, r=5.)

norm = ImageNormalize(stretch=SqrtStretch())
plt.imshow(data, cmap='Greys', origin='lower', norm=norm)
apertures.plot(color='blue', lw=1.5, alpha=0.5)










def main(args):
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
