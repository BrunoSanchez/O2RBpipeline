#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  stuffskywrapper.py
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
"""
Program that wraps SkyMaker and Stuff from Bertin, to
simulate astronomical images.

"""

import os
import shutil
import shlex
import subprocess

import jinja2
import numpy as np
from astropy.io import ascii

def write_stuffconf(dest_file, stuffconf_dict):
    """Writes a stuff configuration file using a template given in
    'stuffconf.j2', using jinja2 template rendering.

    Everything is parsed by a dictionary with the desired values settled.
    """
    loader = jinja2.FileSystemLoader(os.path.abspath('.'))
    jenv = jinja2.Environment(loader=loader,
                              trim_blocks=True,
                              lstrip_blocks=True)

    template = jenv.get_template('rbogus/scripts/templates/stuffconf.j2')

    #~ if os.path.exists(stuffconf_dict['cat_name']):
        #~ os.remove(stuffconf_dict['cat_name'])

    with open(dest_file, 'w') as f1:
        f1.write(template.render(stuffconf_dict))
    return


def write_skyconf(dest_file, skyconf_dict):
    """Writes a sky configuration file using a template given in
    'skyconf.j2', using jinja2 template rendering.

    Everything is parsed by a dictionary with the desired values settled.
    """
    loader = jinja2.FileSystemLoader(os.path.abspath('.'))
    jenv = jinja2.Environment(loader=loader,
                              trim_blocks=True,
                              lstrip_blocks=True)

    template = jenv.get_template('rbogus/scripts/templates/skyconf.j2')

    with open(dest_file, 'w') as f1:
        f1.write(template.render(skyconf_dict))
    return


def run_stuff(stuffconf):
    cmd = "stuff -c {conf} ".format(conf=stuffconf)
    cmd = shlex.split(cmd)

    subprocess.call(cmd)
    return

def run_sky(skyconf, stuff_cat='', img_path=None, t_exp=None):
    if img_path is None:
    #    img_path = skyconf['image_name']
        cmd = "sky "+stuff_cat+" -c {conf} ".format(conf=skyconf)
    else:
        cmd = "sky "+stuff_cat+" -c {conf} -IMAGE_NAME {img_path}".format(
                                                                conf=skyconf,
                                                            img_path=img_path)
    if t_exp is not None:
        cmd = cmd + " -EXPOSURE_TIME {}".format(t_exp)

    cmd = shlex.split(cmd)

    subprocess.call(cmd)
    return img_path


def run_sex(sexconf, img_path=None, cat_output=None):
    print os.path.abspath('.')
    if img_path is None:
        print 'Provide an image path'
        return

    else:
        cmd = "sextractor {img_path} -c {conf} ".format(conf=sexconf,
                                                 img_path=img_path)
    if cat_output is not None:
        cmd = cmd + " -CATALOG_NAME {}".format(cat_output)

    cmd = shlex.split(cmd)

    subprocess.call(cmd)
    return cat_output



def main(args):
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
