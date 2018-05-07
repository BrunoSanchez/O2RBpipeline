#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  ps_query.py
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

import requests
from astropy.io.votable import parse_single_table

def panstarrs_query(ra_deg, dec_deg, rad_deg, mindet=1,
                    maxsources=10000,
                    server='http://archive.stsci.edu/panstarrs/search.php',
                    path='panstarrs.xml'):
    """
    Query Pan-STARRS DR1 @ MAST
    parameters: ra_deg, dec_deg, rad_deg: RA, Dec, field
                                          radius in degrees
                mindet: minimum number of detection (optional)
                maxsources: maximum number of sources
                server: servername
    returns: astropy.table object
    """
    r = requests.get(server,
    params= {'RA': ra_deg, 'DEC': dec_deg,
             'SR': rad_deg, 'max_records': maxsources,
             'outputformat': 'VOTable',
             'ndetections': ('>%d' % mindet)})

    # write query data into local file
    outf = open(path, 'w')
    outf.write(r.text)
    outf.close()

    # parse local file into astropy.table object
    data = parse_single_table(path)
    return data.to_table(use_names_over_ids=True)



def main(ra_deg, dec_deg, rad_deg, mindet=1, maxsources=100, path=None):
    return panstarrs_query(ra_deg, dec_deg, rad_deg,
                           mindet=mindet,
                           maxsources=maxsources,
                           path=path)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("ra_deg", help="RA in degrees", type=float)
    parser.add_argument("dec_deg", help="Dec in degrees", type=float)
    parser.add_argument("rad_deg", help="Rdius of search in degrees", type=float)
    parser.add_argument("--mindet", help="minimum number of detections",
                        type=int, default=1)
    parser.add_argument("--maxsources", help="Maximum number of sources",
                        type=int, default=100)
    parser.add_argument("--path", help="Path for query storage",
                        type=str, default='panstarrs.xml')
    args = parser.parse_args()

    main(args.ra_deg, args.dec_deg, args.rad_deg, args.mindet, args.maxsources,
         args.path)
    print 'Catalog stored in {}'.format(args.path)


# Example query
#print(panstarrs_query(12.345, 67.89, 0.1))
