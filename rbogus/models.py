#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This module contains the table model definitions for toritos database.
# It uses the sqlalchemy orm to define the column structures

import os
import shutil

from corral import db
from corral.conf import settings


class Object(db.Model):
    """Model for galaxies observed"""

    __tablename__ = "Object"

    id = db.Column(db.Integer, primary_key=True)

    name = db.Column(db.String(12), nullable=False)

    def __repr__(self):
        return self.name



class RefImages(db.Model):

    __tablename__ = "RefImages"

    id = db.Column(db.Integer, primary_key=True)

    path = db.Column(db.String(100), nullable=False)

    object_id = db.Column(db.Integer, db.ForeignKey('Object.id'))
    object = db.relationship('Object',
                            backref=db.backref('refimage', order_by=id))

    def __repr__(self):
        return self.path


class NewImages(db.Model):

    __tablename__ = "NewImages"

    id = db.Column(db.Integer, primary_key=True)

    path = db.Column(db.String(100), nullable=False)

    object_id = db.Column(db.Integer, db.ForeignKey('Object.id'))
    object = db.relationship('Object',
                            backref=db.backref('newimage', order_by=id))

    def __repr__(self):
        return self.path


class Simulated(db.Model):
    """Model for transients that are simulated"""

    __tablename__ = "Simulated"

    id = db.Column(db.Integer, primary_key=True)

    x = db.Column(db.Float, nullable=False)
    y = db.Column(db.Float, nullable=False)
    app_mag = db.Column(db.Float, nullable=False)
    flux = db.Column(db.Float, nullable=False)

    image_id = db.Column(db.Integer, db.ForeignKey('Images.id'))
    image = db.relationship('Images',
                            backref=db.backref('simulateds', order_by=id))

    simage_id = db.Column(db.Integer, db.ForeignKey('SImages.id'))
    simage = db.relationship('SImages',
                             backref=db.backref('simulateds', order_by=id))

    scorrimage_id = db.Column(db.Integer, db.ForeignKey('SCorrImages.id'))
    scorrimage = db.relationship('SCorrImages',
                             backref=db.backref('simulateds', order_by=id))

    image_id_ois = db.Column(db.Integer, db.ForeignKey('ImagesOIS.id'))
    image_ois = db.relationship('ImagesOIS',
                            backref=db.backref('simulateds', order_by=id))

    image_id_hot = db.Column(db.Integer, db.ForeignKey('ImagesHOT.id'))
    image_hot = db.relationship('ImagesHOT',
                            backref=db.backref('simulateds', order_by=id))

    def __repr__(self):
        return str(self.id)


# =============================================================================
# Propersubtract tables
# =============================================================================
class Images(db.Model):

    __tablename__ = "Images"

    id = db.Column(db.Integer, primary_key=True)

    path = db.Column(db.String(100), nullable=False)

    crossmatched = db.Column(db.Boolean, nullable=False)

    refimage_id = db.Column(db.Integer, db.ForeignKey('RefImages.id'))
    refimage = db.relationship('RefImages',
                             backref=db.backref('Images', order_by=id))

    newimage_id = db.Column(db.Integer, db.ForeignKey('NewImages.id'))
    newimage = db.relationship('NewImages',
                             backref=db.backref('Images', order_by=id))

    exec_time = db.Column(db.Float, nullable=False)

    def __repr__(self):
        return self.path


class Reals(db.Model):

    __tablename__ = "Reals"

    id = db.Column(db.Integer, primary_key=True)

    detected_id = db.Column(db.Integer, db.ForeignKey('Detected.id'))
    detected = db.relationship('Detected',
                               backref=db.backref('true_pos', order_by=id))

    simulated_id = db.Column(db.Integer, db.ForeignKey('Simulated.id'))
    simulated = db.relationship('Simulated',
                                backref=db.backref('true_pos'), order_by=id)

    def __repr__(self):
        return str(self.id)


class Bogus(db.Model):

    __tablename__ = "Bogus"

    id = db.Column(db.Integer, primary_key=True)

    detected_id = db.Column(db.Integer, db.ForeignKey('Detected.id'))
    detected = db.relationship('Detected',
                               backref=db.backref('true_neg', order_by=id))

    def __repr__(self):
        return str(self.id)


class Undetected(db.Model):

    __tablename__ = "Undetected"

    id = db.Column(db.Integer, primary_key=True)

    simulated_id = db.Column(db.Integer, db.ForeignKey('Simulated.id'))
    simulated = db.relationship('Simulated',
                                backref=db.backref('false_neg', order_by=id))

    def __repr__(self):
        return str(self.id)


class Detected(db.Model):

    __tablename__ = "Detected"

    id = db.Column(db.Integer, primary_key=True)

    NUMBER = db.Column(db.Integer, nullable=False)
    FLUX_ISO = db.Column(db.Float, nullable=False)
    FLUXERR_ISO = db.Column(db.Float, nullable=False)
    MAG_ISO = db.Column(db.Float, nullable=False)
    MAGERR_ISO = db.Column(db.Float, nullable=False)
    FLUX_APER = db.Column(db.Float, nullable=False)
    FLUXERR_APER = db.Column(db.Float, nullable=False)
    MAG_APER = db.Column(db.Float, nullable=False)
    MAGERR_APER = db.Column(db.Float, nullable=False)
    FLUX_AUTO = db.Column(db.Float, nullable=False)
    FLUXERR_AUTO = db.Column(db.Float, nullable=False)
    MAG_AUTO = db.Column(db.Float, nullable=False)
    MAGERR_AUTO = db.Column(db.Float, nullable=False)
    BACKGROUND = db.Column(db.Float, nullable=False)
    THRESHOLD = db.Column(db.Float, nullable=False)
    FLUX_MAX = db.Column(db.Float, nullable=False)
    XMIN_IMAGE = db.Column(db.Float, nullable=False)
    YMIN_IMAGE = db.Column(db.Float, nullable=False)
    XMAX_IMAGE = db.Column(db.Float, nullable=False)
    YMAX_IMAGE = db.Column(db.Float, nullable=False)
    XPEAK_IMAGE = db.Column(db.Float, nullable=False)
    YPEAK_IMAGE = db.Column(db.Float, nullable=False)
    X_IMAGE = db.Column(db.Float, nullable=False)
    Y_IMAGE = db.Column(db.Float, nullable=False)
    X2_IMAGE = db.Column(db.Float, nullable=False)
    Y2_IMAGE = db.Column(db.Float, nullable=False)
    XY_IMAGE = db.Column(db.Float, nullable=False)
    CXX_IMAGE = db.Column(db.Float, nullable=False)
    CYY_IMAGE = db.Column(db.Float, nullable=False)
    CXY_IMAGE = db.Column(db.Float, nullable=False)
    A_IMAGE = db.Column(db.Float, nullable=False)
    B_IMAGE = db.Column(db.Float, nullable=False)
    THETA_IMAGE = db.Column(db.Float, nullable=False)
    MU_MAX = db.Column(db.Float, nullable=False)
    FLAGS = db.Column(db.Float, nullable=False)
    FWHM_IMAGE = db.Column(db.Float, nullable=False)
    ELONGATION = db.Column(db.Float, nullable=False)
    ELLIPTICITY = db.Column(db.Float, nullable=False)
    CLASS_STAR = db.Column(db.Float, nullable=False)
    MU_THRESHOLD = db.Column(db.Float, nullable=False)
    SNR_WIN = db.Column(db.Float, nullable=False)

    DELTAX = db.Column(db.Float, nullable=False)
    DELTAY = db.Column(db.Float, nullable=False)
    RATIO = db.Column(db.Float, nullable=False)
    ROUNDNESS = db.Column(db.Float, nullable=False)
    PEAK_CENTROID = db.Column(db.Float, nullable=False)
    IS_REAL = db.Column(db.Boolean, nullable=True)

    image_id = db.Column(db.Integer, db.ForeignKey('Images.id'))
    image = db.relationship('Images',
                            backref=db.backref('detected_srcs', order_by=id))

    def __repr__(self):
        return '{}::{}'.format(self.image, self.NUMBER)

# =============================================================================
#  SDetect Propersubtract tables
# =============================================================================
class SImages(db.Model):

    __tablename__ = "SImages"

    id = db.Column(db.Integer, primary_key=True)

    path = db.Column(db.String(100), nullable=False)

    crossmatched = db.Column(db.Boolean, nullable=False)


    refimage_id = db.Column(db.Integer, db.ForeignKey('RefImages.id'))
    refimage = db.relationship('RefImages',
                             backref=db.backref('SImages', order_by=id))

    newimage_id = db.Column(db.Integer, db.ForeignKey('NewImages.id'))
    newimage = db.relationship('NewImages',
                             backref=db.backref('SImages', order_by=id))

    exec_time = db.Column(db.Float, nullable=False)

    def __repr__(self):
        return self.path


class SReals(db.Model):

    __tablename__ = "SReals"

    id = db.Column(db.Integer, primary_key=True)

    detected_id = db.Column(db.Integer, db.ForeignKey('SDetected.id'))
    detected = db.relationship('SDetected',
                               backref=db.backref('true_pos_s', order_by=id))

    simulated_id = db.Column(db.Integer, db.ForeignKey('Simulated.id'))
    simulated = db.relationship('Simulated',
                                backref=db.backref('true_pos_s'), order_by=id)

    def __repr__(self):
        return str(self.id)


class SBogus(db.Model):

    __tablename__ = "SBogus"

    id = db.Column(db.Integer, primary_key=True)

    detected_id = db.Column(db.Integer, db.ForeignKey('SDetected.id'))
    detected = db.relationship('SDetected',
                               backref=db.backref('true_neg_s', order_by=id))

    def __repr__(self):
        return str(self.id)


class SUndetected(db.Model):

    __tablename__ = "SUndetected"

    id = db.Column(db.Integer, primary_key=True)

    simulated_id = db.Column(db.Integer, db.ForeignKey('Simulated.id'))
    simulated = db.relationship('Simulated',
                                backref=db.backref('false_neg_s', order_by=id))

    def __repr__(self):
        return str(self.id)


class SDetected(db.Model):

    __tablename__ = "SDetected"

    id = db.Column(db.Integer, primary_key=True)

    thresh = db.Column(db.Float, nullable=False)
    npix = db.Column(db.Float, nullable=False)
    tnpix = db.Column(db.Float, nullable=False)
    xmin = db.Column(db.Float, nullable=False)
    xmax = db.Column(db.Float, nullable=False)
    ymin = db.Column(db.Float, nullable=False)
    ymax = db.Column(db.Float, nullable=False)
    x = db.Column(db.Float, nullable=False)
    y = db.Column(db.Float, nullable=False)
    x2 = db.Column(db.Float, nullable=False)
    y2 = db.Column(db.Float, nullable=False)
    xy = db.Column(db.Float, nullable=False)
    errx2 = db.Column(db.Float, nullable=False)
    erry2 = db.Column(db.Float, nullable=False)
    errxy = db.Column(db.Float, nullable=False)
    a = db.Column(db.Float, nullable=False)
    b = db.Column(db.Float, nullable=False)
    theta = db.Column(db.Float, nullable=False)
    cxx = db.Column(db.Float, nullable=False)
    cyy = db.Column(db.Float, nullable=False)
    cxy = db.Column(db.Float, nullable=False)
    cflux = db.Column(db.Float, nullable=False)
    flux = db.Column(db.Float, nullable=False)
    cpeak = db.Column(db.Float, nullable=False)
    peak = db.Column(db.Float, nullable=False)
    xcpeak = db.Column(db.Float, nullable=False)
    ycpeak = db.Column(db.Float, nullable=False)
    xpeak = db.Column(db.Float, nullable=False)
    ypeak = db.Column(db.Float, nullable=False)
    flag = db.Column(db.Float, nullable=False)
    MU_THRESHOLD = db.Column(db.Float, nullable=False)
    SNR_WIN = db.Column(db.Float, nullable=False)
    DELTAX = db.Column(db.Float, nullable=False)
    DELTAY = db.Column(db.Float, nullable=False)
    RATIO = db.Column(db.Float, nullable=False)
    ROUNDNESS = db.Column(db.Float, nullable=False)
    PEAK_CENTROID = db.Column(db.Float, nullable=False)
    IS_REAL = db.Column(db.Boolean, nullable=True)

    image_id = db.Column(db.Integer, db.ForeignKey('SImages.id'))
    image = db.relationship('SImages',
                            backref=db.backref('Sdetected_srcs', order_by=id))

    def __repr__(self):
        return '{}::{}'.format(self.image, self.id)

# =============================================================================
#  SCOrrDetect Propersubtract tables
# =============================================================================
class SCorrImages(db.Model):

    __tablename__ = "SCorrImages"

    id = db.Column(db.Integer, primary_key=True)

    path = db.Column(db.String(100), nullable=False)

    crossmatched = db.Column(db.Boolean, nullable=False)

    refimage_id = db.Column(db.Integer, db.ForeignKey('RefImages.id'))
    refimage = db.relationship('RefImages',
                             backref=db.backref('SCorrImages', order_by=id))

    newimage_id = db.Column(db.Integer, db.ForeignKey('NewImages.id'))
    newimage = db.relationship('NewImages',
                             backref=db.backref('SCorrImages', order_by=id))
    exec_time = db.Column(db.Float, nullable=False)

    def __repr__(self):
        return self.path


class SCorrReals(db.Model):

    __tablename__ = "SCorrReals"

    id = db.Column(db.Integer, primary_key=True)

    detected_id = db.Column(db.Integer, db.ForeignKey('SCorrDetected.id'))
    detected = db.relationship('SCorrDetected',
                               backref=db.backref('true_pos_scorr', order_by=id))

    simulated_id = db.Column(db.Integer, db.ForeignKey('Simulated.id'))
    simulated = db.relationship('Simulated',
                                backref=db.backref('true_pos_scorr'), order_by=id)

    def __repr__(self):
        return str(self.id)


class SCorrBogus(db.Model):

    __tablename__ = "SCorrBogus"

    id = db.Column(db.Integer, primary_key=True)

    detected_id = db.Column(db.Integer, db.ForeignKey('SCorrDetected.id'))
    detected = db.relationship('SCorrDetected',
                               backref=db.backref('true_neg_scorr', order_by=id))

    def __repr__(self):
        return str(self.id)


class SCorrUndetected(db.Model):

    __tablename__ = "SCorrUndetected"

    id = db.Column(db.Integer, primary_key=True)

    simulated_id = db.Column(db.Integer, db.ForeignKey('Simulated.id'))
    simulated = db.relationship('Simulated',
                                backref=db.backref('false_neg_scorr', order_by=id))

    def __repr__(self):
        return str(self.id)


class SCorrDetected(db.Model):

    __tablename__ = "SCorrDetected"

    id = db.Column(db.Integer, primary_key=True)

    X_IMAGE = db.Column(db.Float, nullable=False)
    Y_IMAGE = db.Column(db.Float, nullable=False)
    SIGNIFICANCE = db.Column(db.Float, nullable=False)

    IS_REAL = db.Column(db.Boolean, nullable=True)

    image_id = db.Column(db.Integer, db.ForeignKey('SCorrImages.id'))
    image = db.relationship('SCorrImages',
                            backref=db.backref('SCorrdetected_srcs', order_by=id))

    def __repr__(self):
        return '{}::{}'.format(self.image, self.id)


# =============================================================================
# OIS tables
# =============================================================================
class ImagesOIS(db.Model):

    __tablename__ = "ImagesOIS"

    id = db.Column(db.Integer, primary_key=True)

    path = db.Column(db.String(100), nullable=False)

    crossmatched = db.Column(db.Boolean, nullable=False)

    refimage_id = db.Column(db.Integer, db.ForeignKey('RefImages.id'))
    refimage = db.relationship('RefImages',
                             backref=db.backref('OISImages', order_by=id))

    newimage_id = db.Column(db.Integer, db.ForeignKey('NewImages.id'))
    newimage = db.relationship('NewImages',
                             backref=db.backref('OISImages', order_by=id))
    exec_time = db.Column(db.Float, nullable=False)

    def __repr__(self):
        return self.path


class DetectedOIS(db.Model):

    __tablename__ = "DetectedOIS"

    id = db.Column(db.Integer, primary_key=True)

    NUMBER = db.Column(db.Integer, nullable=False)
    FLUX_ISO = db.Column(db.Float, nullable=False)
    FLUXERR_ISO = db.Column(db.Float, nullable=False)
    MAG_ISO = db.Column(db.Float, nullable=False)
    MAGERR_ISO = db.Column(db.Float, nullable=False)
    FLUX_APER = db.Column(db.Float, nullable=False)
    FLUXERR_APER = db.Column(db.Float, nullable=False)
    MAG_APER = db.Column(db.Float, nullable=False)
    MAGERR_APER = db.Column(db.Float, nullable=False)
    FLUX_AUTO = db.Column(db.Float, nullable=False)
    FLUXERR_AUTO = db.Column(db.Float, nullable=False)
    MAG_AUTO = db.Column(db.Float, nullable=False)
    MAGERR_AUTO = db.Column(db.Float, nullable=False)
    BACKGROUND = db.Column(db.Float, nullable=False)
    THRESHOLD = db.Column(db.Float, nullable=False)
    FLUX_MAX = db.Column(db.Float, nullable=False)
    XMIN_IMAGE = db.Column(db.Float, nullable=False)
    YMIN_IMAGE = db.Column(db.Float, nullable=False)
    XMAX_IMAGE = db.Column(db.Float, nullable=False)
    YMAX_IMAGE = db.Column(db.Float, nullable=False)
    XPEAK_IMAGE = db.Column(db.Float, nullable=False)
    YPEAK_IMAGE = db.Column(db.Float, nullable=False)
    X_IMAGE = db.Column(db.Float, nullable=False)
    Y_IMAGE = db.Column(db.Float, nullable=False)
    X2_IMAGE = db.Column(db.Float, nullable=False)
    Y2_IMAGE = db.Column(db.Float, nullable=False)
    XY_IMAGE = db.Column(db.Float, nullable=False)
    CXX_IMAGE = db.Column(db.Float, nullable=False)
    CYY_IMAGE = db.Column(db.Float, nullable=False)
    CXY_IMAGE = db.Column(db.Float, nullable=False)
    A_IMAGE = db.Column(db.Float, nullable=False)
    B_IMAGE = db.Column(db.Float, nullable=False)
    THETA_IMAGE = db.Column(db.Float, nullable=False)
    MU_MAX = db.Column(db.Float, nullable=False)
    FLAGS = db.Column(db.Float, nullable=False)
    FWHM_IMAGE = db.Column(db.Float, nullable=False)
    ELONGATION = db.Column(db.Float, nullable=False)
    ELLIPTICITY = db.Column(db.Float, nullable=False)
    CLASS_STAR = db.Column(db.Float, nullable=False)
    MU_THRESHOLD = db.Column(db.Float, nullable=False)
    SNR_WIN = db.Column(db.Float, nullable=False)

    DELTAX = db.Column(db.Float, nullable=False)
    DELTAY = db.Column(db.Float, nullable=False)
    RATIO = db.Column(db.Float, nullable=False)
    ROUNDNESS = db.Column(db.Float, nullable=False)
    PEAK_CENTROID = db.Column(db.Float, nullable=False)
    IS_REAL = db.Column(db.Boolean, nullable=True)

    image_id = db.Column(db.Integer, db.ForeignKey('ImagesOIS.id'))
    image = db.relationship('ImagesOIS',
                            backref=db.backref('detected_srcs', order_by=id))

    def __repr__(self):
        return '{}::{}'.format(self.image, self.NUMBER)


class RealsOIS(db.Model):

    __tablename__ = "RealsOIS"

    id = db.Column(db.Integer, primary_key=True)

    detected_id = db.Column(db.Integer, db.ForeignKey('DetectedOIS.id'))
    detected = db.relationship('DetectedOIS',
                               backref=db.backref('true_pos_ois', order_by=id))

    simulated_id = db.Column(db.Integer, db.ForeignKey('Simulated.id'))
    simulated = db.relationship('Simulated',
                                backref=db.backref('true_pos_ois'), order_by=id)

    def __repr__(self):
        return str(self.id)


class UndetectedOIS(db.Model):

    __tablename__ = "UndetectedOIS"

    id = db.Column(db.Integer, primary_key=True)

    simulated_id = db.Column(db.Integer, db.ForeignKey('Simulated.id'))
    simulated = db.relationship('Simulated',
                                backref=db.backref('false_neg_ois', order_by=id))

    def __repr__(self):
        return str(self.id)


class BogusOIS(db.Model):

    __tablename__ = "BogusOIS"

    id = db.Column(db.Integer, primary_key=True)

    detected_id = db.Column(db.Integer, db.ForeignKey('DetectedOIS.id'))
    detected = db.relationship('DetectedOIS',
                               backref=db.backref('true_neg_ois', order_by=id))

    def __repr__(self):
        return str(self.id)


# =============================================================================
# HOT tables
# =============================================================================
class ImagesHOT(db.Model):

    __tablename__ = "ImagesHOT"

    id = db.Column(db.Integer, primary_key=True)

    path = db.Column(db.String(100), nullable=False)

    crossmatched = db.Column(db.Boolean, nullable=False)

    refimage_id = db.Column(db.Integer, db.ForeignKey('RefImages.id'))
    refimage = db.relationship('RefImages',
                             backref=db.backref('ImagesHOT', order_by=id))

    newimage_id = db.Column(db.Integer, db.ForeignKey('NewImages.id'))
    newimage = db.relationship('NewImages',
                             backref=db.backref('ImagesHOT', order_by=id))

    exec_time = db.Column(db.Float, nullable=False)

    def __repr__(self):
        return self.path


class DetectedHOT(db.Model):

    __tablename__ = "DetectedHOT"

    id = db.Column(db.Integer, primary_key=True)

    NUMBER = db.Column(db.Integer, nullable=False)
    FLUX_ISO = db.Column(db.Float, nullable=False)
    FLUXERR_ISO = db.Column(db.Float, nullable=False)
    MAG_ISO = db.Column(db.Float, nullable=False)
    MAGERR_ISO = db.Column(db.Float, nullable=False)
    FLUX_APER = db.Column(db.Float, nullable=False)
    FLUXERR_APER = db.Column(db.Float, nullable=False)
    MAG_APER = db.Column(db.Float, nullable=False)
    MAGERR_APER = db.Column(db.Float, nullable=False)
    FLUX_AUTO = db.Column(db.Float, nullable=False)
    FLUXERR_AUTO = db.Column(db.Float, nullable=False)
    MAG_AUTO = db.Column(db.Float, nullable=False)
    MAGERR_AUTO = db.Column(db.Float, nullable=False)
    BACKGROUND = db.Column(db.Float, nullable=False)
    THRESHOLD = db.Column(db.Float, nullable=False)
    FLUX_MAX = db.Column(db.Float, nullable=False)
    XMIN_IMAGE = db.Column(db.Float, nullable=False)
    YMIN_IMAGE = db.Column(db.Float, nullable=False)
    XMAX_IMAGE = db.Column(db.Float, nullable=False)
    YMAX_IMAGE = db.Column(db.Float, nullable=False)
    XPEAK_IMAGE = db.Column(db.Float, nullable=False)
    YPEAK_IMAGE = db.Column(db.Float, nullable=False)
    X_IMAGE = db.Column(db.Float, nullable=False)
    Y_IMAGE = db.Column(db.Float, nullable=False)
    X2_IMAGE = db.Column(db.Float, nullable=False)
    Y2_IMAGE = db.Column(db.Float, nullable=False)
    XY_IMAGE = db.Column(db.Float, nullable=False)
    CXX_IMAGE = db.Column(db.Float, nullable=False)
    CYY_IMAGE = db.Column(db.Float, nullable=False)
    CXY_IMAGE = db.Column(db.Float, nullable=False)
    A_IMAGE = db.Column(db.Float, nullable=False)
    B_IMAGE = db.Column(db.Float, nullable=False)
    THETA_IMAGE = db.Column(db.Float, nullable=False)
    MU_MAX = db.Column(db.Float, nullable=False)
    FLAGS = db.Column(db.Float, nullable=False)
    FWHM_IMAGE = db.Column(db.Float, nullable=False)
    ELONGATION = db.Column(db.Float, nullable=False)
    ELLIPTICITY = db.Column(db.Float, nullable=False)
    CLASS_STAR = db.Column(db.Float, nullable=False)
    MU_THRESHOLD = db.Column(db.Float, nullable=False)
    SNR_WIN = db.Column(db.Float, nullable=False)

    DELTAX = db.Column(db.Float, nullable=False)
    DELTAY = db.Column(db.Float, nullable=False)
    RATIO = db.Column(db.Float, nullable=False)
    ROUNDNESS = db.Column(db.Float, nullable=False)
    PEAK_CENTROID = db.Column(db.Float, nullable=False)
    IS_REAL = db.Column(db.Boolean, nullable=True)

    image_id = db.Column(db.Integer, db.ForeignKey('ImagesHOT.id'))
    image = db.relationship('ImagesHOT',
                            backref=db.backref('detected_srcs', order_by=id))

    def __repr__(self):
        return '{}::{}'.format(self.image, self.NUMBER)


class RealsHOT(db.Model):

    __tablename__ = "RealsHOT"

    id = db.Column(db.Integer, primary_key=True)

    detected_id = db.Column(db.Integer, db.ForeignKey('DetectedHOT.id'))
    detected = db.relationship('DetectedHOT',
                               backref=db.backref('true_pos_hot', order_by=id))

    simulated_id = db.Column(db.Integer, db.ForeignKey('Simulated.id'))
    simulated = db.relationship('Simulated',
                                backref=db.backref('true_pos_hot'), order_by=id)

    def __repr__(self):
        return str(self.id)


class UndetectedHOT(db.Model):

    __tablename__ = "UndetectedHOT"

    id = db.Column(db.Integer, primary_key=True)

    simulated_id = db.Column(db.Integer, db.ForeignKey('Simulated.id'))
    simulated = db.relationship('Simulated',
                                backref=db.backref('false_neg_hot', order_by=id))

    def __repr__(self):
        return str(self.id)


class BogusHOT(db.Model):

    __tablename__ = "BogusHOT"

    id = db.Column(db.Integer, primary_key=True)

    detected_id = db.Column(db.Integer, db.ForeignKey('DetectedHOT.id'))
    detected = db.relationship('DetectedHOT',
                               backref=db.backref('true_neg_hot', order_by=id))

    def __repr__(self):

        return str(self.id)
