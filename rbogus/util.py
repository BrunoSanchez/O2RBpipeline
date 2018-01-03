#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""Module with simple tools for specific
rbogus data handling routines.
"""
import numpy as np
from astroML import crossmatch as cx


def matching(detected, simulated, radius=1., masked=False, sep=False):
    """Function to match stars between frames.
    """

    masterXY = np.empty((len(detected), 2), dtype=np.float64)
    master_idx = np.zeros(len(detected))
    for i in range(len(detected)):
        if sep:
            masterXY[i, 0] = detected[i].x
            masterXY[i, 1] = detected[i].y
        else:
            masterXY[i, 0] = detected[i].X_IMAGE
            masterXY[i, 1] = detected[i].Y_IMAGE
        master_idx[i] = detected[i].id


    imXY = np.empty((len(simulated), 2), dtype=np.float64)
    for i in range(len(simulated)):
        imXY[i, 0] = simulated[i].x
        imXY[i, 1] = simulated[i].y

    try:
        dist, ind = cx.crossmatch(masterXY, imXY, max_distance=radius)
        dist_, ind_ = cx.crossmatch(imXY, masterXY, max_distance=radius)
        match = ~np.isinf(dist)
        match_ = ~np.isinf(dist_)
    except ValueError:
        print "masterXY.shape"
        print masterXY.shape
        print "\n\nimXY.shape"
        print imXY.shape
        print "\n"
        try:
            dist, ind = cx.crossmatch(masterXY, imXY, max_distance=radius*3)
            dist_, ind_ = cx.crossmatch(imXY, masterXY, max_distance=radius*3)
            match = ~np.isinf(dist)
            match_ = ~np.isinf(dist_)
        except:
            # import ipdb; ipdb.set_trace()
            print 'There were cross match exceptions'
            return np.zeros_like(imXY) - 13

    IDs = np.zeros_like(ind_) - 13133
    # IDs has length = len(imXY) = len(simulated)

    for i in xrange(len(ind_)):
        if dist_[i] != np.inf:
            dist_o = dist_[i]
            ind_o = ind_[i]

            if dist[ind_o] != np.inf:
                dist_s = dist[ind_o]
                ind_s = ind[ind_o]

                if ind_s == i:
                    IDs[i] = master_idx[ind_o]
    #                print master_idx[ind_o], ind_o

    print len(IDs), len(ind_), len(ind)
    if masked:
        mask = IDs > 0
        return(IDs, mask)
    return(IDs)

