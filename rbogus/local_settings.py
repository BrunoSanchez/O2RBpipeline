"""Local settings for rbpipeline"""

IMGS_PATH = '/home/bruno/Data/NewRBpipeline/images'
CONFIG_PATH = '/home/bruno/Data/NewRBpipeline/config'
CATS_PATH = '/home/bruno/Data/NewRBpipeline/cats'


SIM_CUBE = []

for zp in [4e3, 8e3, 16e3, 32e3, 64e3]:
    for slope in [0.1, 0.3, 0.5, 0.7, 0.9]:
        for fwhm in [0.8, 1., 1.2, 1.4]:
            SIM_CUBE.append({'zp': zp, 'slope': slope, 'fwhm': fwhm})
