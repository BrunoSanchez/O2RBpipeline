#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Created at 2016-12-14T10:42:52.224316 by corral 0.0.1


# =============================================================================
# DOCS
# =============================================================================

"""rbogus main loader

"""


# =============================================================================
# IMPORTS
# =============================================================================

from corral import run
from . import models
import json
from scripts import gen_diff
from corral.conf import settings as stgs

# =============================================================================
# LOADER
# =============================================================================

class Load(run.Loader):

    def setup(self):
        self.session.autocommit = False
        self.current_index = self.session.query(models.Images).order_by(
                                models.Images.id.desc()).first()
        #~ self.session.buff = []

    def generate(self):
        ref_path = stgs.REFERENCE_IMAGE
        new_path = stgs.NEW_IMAGE
        with open(stgs.DETAILS_FILE) as fp:
            details = json.load(fp)

        obj = self.session.query(models.Object).filter(
                models.Object.name==details['object']).first()
        if obj is None:
            obj = models.Object()
            obj.name = details['object']
            self.session.add(obj)
            self.session.commit()

        ref = models.RefImages()
        ref.object = obj
        ref.path = details['orig_ref_path']

        new = models.NewImages()
        new.object = obj
        new.path = details['orig_new_path']

        self.session.add(ref)
        self.session.add(new)
        self.session.commit()
        if self.current_index is not None:
            index = self.current_index.id
        else:
            index = None
        results = gen_diff.main(ref_path, new_path, details, index)
        #~ results = gen_diff.main(self.current_index, **self.current_params)

        diff_path      = results[0]
        detections     = results[1]
        diff_ois_path  = results[2]
        detections_ois = results[3]
        diff_hot_path  = results[4]
        detections_hot = results[5]
        transients     = results[6]
        sdetections    = results[7]
        scorrdetections= results[8]
        times          = results[9]

# =============================================================================
# properimage
# =============================================================================
        image = models.Images()
        image.path = diff_path
        image.refimage = ref
        image.newimage = new
        image.crossmatched = False
        image.exec_time = times[0]

        self.session.add(image)
        self.session.commit()

        detections['image_id'] = gen_diff.np.repeat(image.id, len(detections))
        detections.to_sql('Detected', self.session.get_bind(),
                           if_exists='append', index=False)

# =============================================================================
# sdetect
# =============================================================================
        simage = models.SImages()
        simage.path = diff_path
        simage.ref = ref
        simage.new = new
        simage.crossmatched = False
        simage.exec_time = times[0]

        self.session.add(simage)
        self.session.commit()

        sdetections['image_id'] = gen_diff.np.repeat(simage.id, len(sdetections))
        sdetections.to_sql('SDetected', self.session.get_bind(),
                            if_exists='append', index=False)
# =============================================================================
# scorrdetect
# =============================================================================
        scorrimage = models.SCorrImages()
        scorrimage.path = diff_path
        scorrimage.ref = ref
        scorrimage.new = new
        scorrimage.crossmatched = False
        scorrimage.exec_time = times[0]

        self.session.add(scorrimage)
        self.session.commit()

        scorrdetections['image_id'] = gen_diff.np.repeat(scorrimage.id, len(scorrdetections))
        scorrdetections.to_sql('SCorrDetected', self.session.get_bind(),
                            if_exists='append', index=False)
#------------------------------------------------------------------------------
# =============================================================================
# OIS
# =============================================================================
        image_ois = models.ImagesOIS()
        image_ois.path = diff_ois_path
        image_ois.ref = ref
        image_ois.new = new
        image_ois.crossmatched = False
        image_ois.exec_time = times[1]

        self.session.add(image_ois)
        self.session.commit()

        detections_ois['image_id'] = gen_diff.np.repeat(image_ois.id,
                                                        len(detections_ois))
        detections_ois.to_sql('DetectedOIS', self.session.get_bind(),
                           if_exists='append', index=False)
#------------------------------------------------------------------------------
# =============================================================================
# HOTPANTS
# =============================================================================
        image_hot = models.ImagesHOT()
        image_hot.path = diff_hot_path
        image_hot.ref = ref
        image_hot.new = new
        image_hot.crossmatched = False
        image_hot.exec_time = times[2]

        self.session.add(image_hot)
        self.session.commit()

        detections_hot['image_id'] = gen_diff.np.repeat(image_hot.id,
                                                        len(detections_hot))
        detections_hot.to_sql('DetectedHOT', self.session.get_bind(),
                           if_exists='append', index=False)
#------------------------------------------------------------------------------

        transients['image_id'] = gen_diff.np.repeat(image.id, len(transients))

        transients['simage_id'] = gen_diff.np.repeat(simage.id, len(transients))

        transients['image_id_ois'] = gen_diff.np.repeat(image_ois.id,
                                                        len(transients))
        transients['image_id_hot'] = gen_diff.np.repeat(image_hot.id,
                                                        len(transients))
        # print transients
        transients.to_sql('Simulated', self.session.get_bind(),
                          if_exists='append', index=False)

    def teardown(self, type, value, traceback):
        if not type:
            self.session.commit()



