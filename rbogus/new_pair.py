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
import json
import astroalign as aa
from astropy.io import fits
from astropy.io import ascii

REFERENCE_IMAGE = '/home/bruno/Data/O2pipeline/input_images/ref.fits'
NEW_IMAGE = '/home/bruno/Data/O2pipeline/input_images/new.fits'


def main(ref_path, new_path, objname):
    # alineo las imagenes
    refdata = fits.getdata(ref_path)
    newdata = fits.getdata(newdata)

    new_aligned = aa.register(newdata, refdata)

    # las copio al lugar designado en la pipeline
    ref_h = fits.getheader(ref_path)
    fits.writeto(data=refdata, header=ref_h, filename=REFERENCE_IMAGE, overwrite=True)

    new_h = fits.getheader(new_path)
    fits.writeto(data=new_aligned, header=new_h, filename=NEW_IMAGE, overwrite=True)

    # creo un file con los detalles
    meta = {}
    meta['object'] = objname
    meta['orig_ref_path'] = ref_path
    meta['orig_new_path'] = new_path

    with open('details.json', 'w') as fp:
        json.dump(meta, fp)

    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
