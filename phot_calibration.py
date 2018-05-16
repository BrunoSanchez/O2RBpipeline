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
from astropy.io import ascii

from astropy import coordinates as coord
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
#image_path = os.path.join(basepath, 'PGC073926/pgc073926_170104.fits')
def cal_image(img_path):
    image_path = os.path.join(basepath, img_path)

    # load img with WCS
    if not '.new' in image_path:
        img = image_path.strip('.fits') + '.new'
    else:
        img = image_path

    fits_file = fits.open(img, mode='update')
    hdr = fits_file[0].header #fits.getheader(img)
    data = fits_file[0].data  # fits.getdata(img)
    wcs = WCS(hdr)

    # Obtain image center
    crval_ra= hdr['CRVAL1']
    crval_dec = hdr['CRVAL2']
    crval_x = hdr['CRPIX1']
    crval_y = hdr['CRPIX2']

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

    daofind = DAOStarFinder(fwhm=5.0, threshold=5.*std)
    sources = daofind(data - median)
    #print(sources)

    positions = (sources['xcentroid'], sources['ycentroid'])
    xypos = np.array([sources['xcentroid'], sources['ycentroid']]).T
    sky_pos = wcs.all_pix2world(xypos, 0)
    skycoords = coord.SkyCoord(sky_pos*u.deg, frame='icrs')

    radii = [5.*u.arcsec, 7.*u.arcsec, 10.*u.arcsec]
    annulus = [(6.*u.arcsec, 8.*u.arcsec),
               (9.*u.arcsec, 11.*u.arcsec),
               (11.*u.arcsec, 12.*u.arcsec)]

    apertures = [SkyCircularAperture(skycoords, r=r).to_pixel(wcs)
                 for r in radii]
    annulii = [SkyCircularAnnulus(skycoords, r_in=ri, r_out=ro).to_pixel(wcs)
               for ri, ro in annulus]
    apers = apertures + annulii

    #~ tables = [aperture_photometry(data-back, aper) for aper in apers]
    phot_table = aperture_photometry(data-back, apers)

    bkg_sum_0 = apertures[0].area()*phot_table['aperture_sum_3']/annulii[0].area()
    bkg_sum_1 = apertures[0].area()*phot_table['aperture_sum_4']/annulii[1].area()
    bkg_sum_2 = apertures[0].area()*phot_table['aperture_sum_5']/annulii[2].area()

    final_sum_0 = phot_table['aperture_sum_0'] - bkg_sum_0
    final_sum_1 = phot_table['aperture_sum_1'] - bkg_sum_1
    final_sum_2 = phot_table['aperture_sum_2'] - bkg_sum_2

    phot_table['resid_aper_sum_0'] = final_sum_0
    phot_table['resid_aper_sum_1'] = final_sum_1
    phot_table['resid_aper_sum_2'] = final_sum_2

    phot_table['cmag_0'] = -2.5 * np.log10(final_sum_0)
    phot_table['cmag_1'] = -2.5 * np.log10(final_sum_1)
    phot_table['cmag_2'] = -2.5 * np.log10(final_sum_2)

    phot_table['mag_0'] = -2.5 * np.log10(phot_table['aperture_sum_0'])
    phot_table['mag_1'] = -2.5 * np.log10(phot_table['aperture_sum_1'])
    phot_table['mag_2'] = -2.5 * np.log10(phot_table['aperture_sum_2'])

    phot_table['RA'] = skycoords.ra
    phot_table['Dec'] = skycoords.dec

    # match the catalogs
    # attempt a xmatch
    phot_table.write('/tmp/sourcecat.csv', overwrite=True)
    try:
        from astroquery.xmatch import XMatch
        XMatch.TIMEOUT = 600
        try:
            table = XMatch.query(cat1=open('/tmp/sourcecat.csv'),
                             #cat2='vizier:II/336/apass9', #apass
                             cat2='vizier:I/322A/out',   # UCAC4
                             #cat2='vizier:II/246/out',   # 2MASS
                             max_distance=5 * u.arcsec,
                             colRA1='RA',
                             colDec1='Dec',
                             colRA2='RAJ2000',
                             colDec2='DEJ2000')
            print 'Using UCAC4'
            print table.colnames
            if 'Rmag' not in table.colnames and 'rmag' in table.colnames:
                table['Rmag'] = table['rmag']
        except:
            try:
                table = XMatch.query(cat1=open('/tmp/sourcecat.csv'),
                                 cat2='vizier:II/336/apass9', #apass
                                 #cat2='vizier:I/322A/out',   # UCAC4
                                 #cat2='vizier:II/246/out',   # 2MASS
                                 max_distance=5 * u.arcsec,
                                 colRA1='RA',
                                 colDec1='Dec',
                                 colRA2='RAJ2000',
                                 colDec2='DEJ2000')
                print 'Using APASS'
                print table.colnames
                if 'Rmag' not in table.colnames and 'rmag' in table.colnames:
                    table['Rmag'] = table['rmag']
            except:

                try:
                    table = XMatch.query(cat1=open('/tmp/sourcecat.csv'),
                                     #cat2='vizier:II/336/apass9', #apass
                                     #cat2='vizier:I/322A/out',   # UCAC4
                                     cat2='vizier:II/246/out',   # 2MASS
                                     max_distance=5 * u.arcsec,
                                     colRA1='RA',
                                     colDec1='Dec',
                                     colRA2='RAJ2000',
                                     colDec2='DEJ2000')
                    print 'Using 2MASS'
                    print table.colnames
                    if 'Rmag' not in table.colnames and 'rmag' in table.colnames:
                        table['Rmag'] = table['rmag']
                except:
                    print('\ntable 1 avaliable?:')
                    print XMatch.is_table_available('vizier:II/336/apass9')
                    print('\ntable 2 avaliable?:')
                    print XMatch.is_table_available('vizier:I/322A/out')
                    print('\ntable 3 avaliable?:')
                    print XMatch.is_table_available('vizier:II/246/out')
                    raise
    except:
        # query PanSTARRS DR1
        #~ catalog = panstarrs_query(crval_ra, crval_dec, 0.5, maxsources=300)
        #~ if len(ps1_tab) < 20:
            #~ from astroquery.vizier import Vizier
            #~ catalog = Vizier.query_region(coord.SkyCoord(crval_ra*u.deg,
                                                         #~ crval_dec*u.deg,
                                                         #~ frame='icrs'),
                                          #~ radius=coord.Angle(0.3, "deg"),
                                          #~ catalog='UCAC4')[0]
        #~ ra = catalog['RAJ2000']
        #~ dec = catalog['DEJ2000']
        #~ cat = coord.SkyCoord(ra, dec, frame='icrs')
        #~ # cat to sources
        #~ idx_c2s, d2d_c2s, _ = skycoords.match_to_catalog_sky(cat)
        #~ # sources to cat
        #~ idx_s2c, d2d_s2c, _ = cat.match_to_catalog_sky(skycoords)
        #~ matches = cat[idx_c2s] # same length than skycoords
        #~ idx_c2s[idx_s2c]
        raise

    from sklearn import linear_model
    if 'UCAC4' in table.colnames:
        table['color'] = table['Bmag'] - table['Vmag']
    if '2MASS' in table.colnames:
        table['color'] = table['Jmag'] - table['Kmag']
        table['Bmag'] = table['Jmag']


    # V Mag
    plt.subplot(221)
    reg = linear_model.LinearRegression()
    sub = table[['mag_1', 'Vmag']].to_pandas()
    sub = sub.dropna(axis=0, how='any')
    X = sub[['mag_1']].as_matrix()
    Y = sub['Vmag'].as_matrix()
    reg.fit(X,Y)
    print reg.coef_, reg.intercept_
    regV = reg
    table['pred_Vmag'] = reg.coef_*table['mag_1']+reg.intercept_
    plt.plot(table['mag_1'], table['Vmag'], 'r.')
    plt.plot(table['mag_1'], table['pred_Vmag'], 'k-', label='linear fit')
    plt.xlabel('mag_1')
    plt.ylabel('Vmag')
    plt.grid()

    # B mag
    plt.subplot(222)
    reg = linear_model.LinearRegression()
    sub = table[['mag_1', 'Bmag']].to_pandas()
    sub = sub.dropna(axis=0, how='any')
    X = sub[['mag_1']].as_matrix()
    Y = sub['Bmag'].as_matrix()
    reg.fit(X,Y)
    print reg.coef_, reg.intercept_
    regB = reg
    table['pred_Bmag'] = reg.coef_*table['mag_1']+reg.intercept_
    plt.plot(table['mag_1'], table['Bmag'], 'r.')
    plt.plot(table['mag_1'], table['pred_Bmag'], 'k-', label='linear fit')
    plt.xlabel('mag_1')
    plt.ylabel('Bmag')
    plt.grid()

    # R mag
    plt.subplot(223)
    reg = linear_model.LinearRegression()
    sub = table[['mag_1', 'Rmag']].to_pandas()
    sub = sub.dropna(axis=0, how='any')
    X = sub[['mag_1']].as_matrix()
    Y = sub['Rmag'].as_matrix()
    reg.fit(X,Y)
    print reg.coef_, reg.intercept_
    regR = reg
    table['pred_Rmag'] = reg.coef_*table['mag_1']+reg.intercept_
    plt.plot(table['mag_1'], table['Rmag'], 'r.')
    plt.plot(table['mag_1'], table['pred_Rmag'], 'k-', label='linear fit')
    plt.xlabel('mag_1')
    plt.ylabel('Rmag')
    plt.grid()

    # Residuals vs B-V
    plt.subplot(224)
    plt.plot(table['color'], table['pred_Vmag']-table['Vmag'], 'g.', label='Vmag')
    plt.plot(table['color'], table['pred_Bmag']-table['Bmag'], 'b.', label='Bmag')
    plt.plot(table['color'], table['pred_Rmag']-table['Rmag'], 'r.', label='Rmag')
    plt.xlabel('color')
    plt.ylabel('Residual')
    plt.legend(loc='best')
    plt.grid()

    plt.tight_layout()
    plt.savefig(image_path.strip('.new')+'_calibration.png')
    plt.clf()

    hdr['ZP_R'] = regR.intercept_
    hdr['ZP_B'] = regB.intercept_
    hdr['ZP_V'] = regV.intercept_

    hdr['SLOPE_R'] = regR.coef_[0]
    hdr['SLOPE_B'] = regB.coef_[0]
    hdr['SLOPE_V'] = regV.coef_[0]

    chosen = {'band': None, 'res': 10}
    for band in ['R', 'B', 'V']:
        res = table['pred_'+band+'mag']-table[band+'mag']
        residual = np.sqrt(np.sum(np.square(res)))
        if residual < chosen['res']:
            chosen['res'] = residual
            chosen['band'] = band

    hdr['CAL_BAND'] = chosen['band']
    hdr['CAL_RES'] = chosen['res']

    fits_file.flush()
    fits_file.close()

    band = chosen['band']
    plt.hist(table['pred_'+band+'mag'])
    plt.xlabel('instrumental '+band+' mag')
    plt.savefig(image_path.strip('.new')+'_magnitudes.png')

    table.write(image_path.strip('.new')+'_sources.csv', overwrite=True)
    with open('calibrations.txt', 'a') as cals:
        cals.write(image_path)
        cals.write(' ')
        cals.write('{}'.format(np.max(table['pred_'+band+'mag'])))
        cals.write(' ')
        cals.write('{}'.format(np.std(table['pred_'+band+'mag'])))
        cals.write(' ')
        cals.write('{}'.format(hdr['ZP_R']))
        cals.write(' ')
        cals.write('{}'.format(hdr['ZP_B']))
        cals.write(' ')
        cals.write('{}'.format(hdr['ZP_V']))
        cals.write(' ')
        cals.write('{}'.format(hdr['SLOPE_R']))
        cals.write(' ')
        cals.write('{}'.format(hdr['SLOPE_B']))
        cals.write(' ')
        cals.write('{}'.format(hdr['SLOPE_V']))
        cals.write(' ')
        cals.write('{}'.format(hdr['CAL_BAND']))
        cals.write(' ')
        cals.write('{}'.format(hdr['CAL_RES']))
        cals.write('\n')

    print 'finished'


if __name__ == '__main__':
    import sys
    sys.exit(cal_image(sys.argv[-1]))
