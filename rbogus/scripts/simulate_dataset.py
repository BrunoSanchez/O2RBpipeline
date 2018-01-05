#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  simulate_dataset.py
#
#  Copyright 2016 Bruno S <bruno@oac.unc.edu.ar>
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
import shutil
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
import sep

import ois
from properimage import single_image as si
from properimage import propersubtract as ps
from properimage import utils

import stuffskywrapper as w

from corral.conf import settings

def main(imgs_dir, ref_path, new_path, details):

    if not os.path.isdir(imgs_dir):
        os.makedirs(imgs_dir)

    ref_dest = os.path.join(imgs_dir, 'ref.fits')
    new_dest = os.path.join(imgs_dir, 'new.fits')

    os.rename(ref_path, ref_dest)
    os.rename(new_path, new_dest)

    ref = si.SingleImage(ref_dest, borders=True) #, crop=((150,150), (150, 150)))
    new = si.SingleImage(new_dest, borders=True) #, crop=((150,150), (150, 150)))

    #~ ##  Adding stars
    foo = new.cov_matrix
    srcs = new.best_sources

    rows = []
    for i in range(12):
        j = np.random.choice(len(srcs), 1, replace=False)
        flux = srcs['flux'][j][0]
        print flux
        star = flux*new.db.load(j)[0]
        sx, sy = star.shape
        x = np.random.choice(new.pixeldata.shape[0]-3*sx, 1)[0] + sx#* np.random.random())
        y = np.random.choice(new.pixeldata.shape[1]-3*sy, 1)[0] + sy
        # np.int((new.pixeldata.shape[1]-star.shape[1]) * np.random.random())

        #~ print x, x + sx
        #~ print y, y + sy
        #~ print new.pixeldata.shape
        #~ print new.pixeldata.data[x:x+sx, y:y+sy].shape
        #~ print star.shape, '\n'
        #~ if new.pixeldata.data[x:x+sx, y:y+sy].shape != star.shape:
            #~ import ipdb; ipdb.set_trace()
        new.pixeldata.data[x:x+sx, y:y+sy] = star
            #~ except:
                #~ continue
        xc = x+sx/2.
        yc = y+sy/2.
        app_mag = -2.5*np.log10(flux)

        rows.append([xc, yc, app_mag, flux])

    newcat = Table(rows=rows, names=['x', 'y', 'app_mag', 'flux'])
    fits.writeto(filename=new_dest, header=fits.getheader(new_dest),
                 data=new.pixeldata.data, overwrite=True)
    new._clean()
    new = si.SingleImage(new_dest, borders=False) #, crop=((150,150), (150, 150)))
    newcat.write(os.path.join(imgs_dir, 'transient.list'),
                           format='ascii.fast_no_header',
                           overwrite=True)

    try:
        print 'Images to be subtracted: {} {}'.format(ref_dest, new_dest)
        import time
        t0 = time.time()
        D, P, S, mask = ps.diff(new, ref, align=False,
                               iterative=False, shift=False, beta=True)
        dt_z = time.time() - t0
        new._clean()
        ref._clean()
        mea, med, std = sigma_clipped_stats(D.real)
        D = np.ma.MaskedArray(D.real, mask).filled(mea)

        fits.writeto(os.path.join(imgs_dir,'diff.fits'), D, overwrite=True)
        #~ utils.encapsule_R(D, path=os.path.join(imgs_dir, 'diff.fits'))
        utils.encapsule_R(P, path=os.path.join(imgs_dir, 'psf_d.fits'))
        utils.encapsule_R(S, path=os.path.join(imgs_dir, 's_diff.fits'))

        scorrdetected = utils.find_S_local_maxima(S, threshold=3.5)
        print 'S_corr found thath {} transients were above 3.5 sigmas'.format(len(scorrdetected))
        ascii.write(table=np.asarray(scorrdetected),
                output=os.path.join(imgs_dir, 's_corr_detected.csv'),
                names=['X_IMAGE', 'Y_IMAGE', 'SIGNIFICANCE'],
                format='csv')

        S = np.ascontiguousarray(S)
        #~ s_bkg = sep.Background(S)
        from astropy.stats import sigma_clipped_stats
        mean, median, std = sigma_clipped_stats(S)
        sdetected = sep.extract(S-median, 3.5*std,
                                filter_kernel=None)
        print 'S_corr with sep found thath {} transients were above 3.5 sigmas'.format(len(sdetected))
        ascii.write(table=sdetected,
                    output=os.path.join(imgs_dir, 'sdetected.csv'),
                    format='csv')

    ##  With OIS
        t0 = time.time()
        ois_d = ois.optimal_system(fits.getdata(new_dest), fits.getdata(ref_dest))[0]
        dt_o = time.time() - t0
        utils.encapsule_R(ois_d, path=os.path.join(imgs_dir, 'diff_ois.fits'))

    ##  With HOTPANTS
        t0 = time.time()
        os.system('hotpants -v 0 -inim {} -tmplim {} -outim {}'.format(new_dest, ref_dest,
            os.path.join(imgs_dir, 'diff_hot.fits')))
        dt_h = time.time() - t0

        return [newcat.to_pandas(), [dt_z, dt_o, dt_h]]
    except:
        raise


if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
