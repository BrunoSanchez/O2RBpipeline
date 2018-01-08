#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  new_pair.py
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
import json
import astroalign as aa
import numpy as np
from astropy.io import fits
from astropy.io import ascii
from properimage import single_image as si

STACK_PATH = '/home/bruno/Data/LIGO_O2/Jan04/newstacks'
REFERENCE_IMAGE = '/home/bruno/Data/O2RBpipeline/input_images/ref.fits'
NEW_IMAGE = '/home/bruno/Data/O2RBpipeline/input_images/new.fits'
DETAILS_FILE = '/home/bruno/Data/O2RBpipeline/input_images/details.json'

aa.PIXEL_TOL=1.7
aa.NUM_NEAREST_NEIGHBORS = 5

def main(ref_path, new_path, objname):
    # alineo las imagenes
    refdata = fits.getdata(os.path.join(STACK_PATH, ref_path))[300:-300, 300:-300]
    newdata = fits.getdata(os.path.join(STACK_PATH, new_path))[300:-300, 300:-300]

    try:
        new_aligned = aa.register(newdata.astype('<f8'),
                                  refdata.astype('<f8'))
        print new_aligned.shape
        print refdata.shape
        if new_aligned.shape != refdata.shape:
            import ipdb; ipdb.set_trace()
            new_aligned = new_aligned[:refdata.shape[0],
                                      :refdata.shape[1]]
    except:
        try:
            aa.MIN_MATCHES_FRACTION = 0.01
            ref = si.SingleImage(refdata.astype('<f8'))
            new = si.SingleImage(newdata.astype('<f8'))

            rs = np.empty((len(ref.best_sources), 2))
            j=0
            for x, y in ref.best_sources[['x', 'y']]:
                rs[j] = x, y
                j += 1

            ns = np.empty((len(new.best_sources), 2))
            j=0
            for x, y in new.best_sources[['x', 'y']]:
                ns[j] = x, y
                j += 1
            trf, _ = aa.find_transform(ns, rs)
            new_aligned = aa.apply_transform(trf, newdata.astype('<f8'),
                                                  refdata.astype('<f8'))
        except:
            raise
    # las copio al lugar designado en la pipeline
    ref_h = fits.getheader(os.path.join(STACK_PATH, ref_path))
    fits.writeto(data=refdata, header=ref_h, filename=REFERENCE_IMAGE, overwrite=True)

    new_h = fits.getheader(os.path.join(STACK_PATH, new_path))
    fits.writeto(data=new_aligned, header=new_h, filename=NEW_IMAGE, overwrite=True)

    # creo un file con los detalles
    meta = {}
    meta['object'] = objname
    meta['orig_ref_path'] = ref_path
    meta['orig_new_path'] = new_path

    with open(DETAILS_FILE, 'w') as fp:
        json.dump(meta, fp)

    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(*sys.argv[1:]))
