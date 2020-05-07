import sys
import os

dirname = os.path.split(__file__)[0]

sys.path.append(dirname + '/../lrg')
from eboss_dr14_lrg import Eboss14LRGLikelihood


class Eboss14QSOLikelihood(Eboss14LRGLikelihood):
    data_type = 'QSO'
    like_name = "eboss14_qso"
    mean = 3843. 
    sigma = 17.0
    fiducial_z = 1.52


setup, execute, cleanup = Eboss14QSOLikelihood.build_module()