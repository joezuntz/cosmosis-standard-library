import numpy as np
import copy
import pickle
import tqdm
import hashlib
from ._version import __version__

def _md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def _transform_space(x, space_rotation=False, rotation=None, bounds=None):
    """Normalize coordinates to [0,1] intervals and if necessary apply a rotation

    :param x: coordinates in parameter space
    :type x: ndarray
    :param space_rotation: whether to apply the rotation matrix defined through
                           the rotation keyword, defaults to False
    :type space_rotation: bool, optional
    :param rotation: rotation matrix, defaults to None
    :type rotation: ndarray, optional
    :param bounds: ranges within which the emulator hypervolume is defined,
                   defaults to None
    :type bounds: ndarray, optional
    :return: normalized and (if required) rotated coordinates
    :rtype: ndarray
    """
    if space_rotation:
        #Get x into the eigenbasis
        R = rotation['rotation_matrix'].T
        xR = copy.deepcopy(np.array([np.dot(R, xi)
                                     for xi in x]))
        xR = xR - rotation['rot_points_means']
        xR = xR/rotation['rot_points_stddevs']
        return xR
    else:
        return (x - bounds[:, 0])/(bounds[:, 1] - bounds[:, 0])

def accuracy_exp_002(y_true, y_pred):
    dataset = K.abs(K.exp(y_pred)/K.exp(y_true)-1)
    tot = dataset >= 0
    sel = dataset <= 0.002
    return K.shape(dataset[sel])[0] /K.shape(dataset[tot])[0]

def accuracy_exp_005(y_true, y_pred):
    dataset = K.abs(K.exp(y_pred)/K.exp(y_true)-1)
    tot = dataset >= 0
    sel = dataset <= 0.005
    return K.shape(dataset[sel])[0] /K.shape(dataset[tot])[0]

def accuracy_exp_01(y_true, y_pred):
        dataset = K.abs(K.exp(y_pred)/K.exp(y_true)-1)
        tot = dataset >= 0
        sel = dataset <= 0.01
        return K.shape(dataset[sel])[0] /K.shape(dataset[tot])[0]

def mean_absolute_exp_percentage_error(y_true, y_pred):
    diff = K.abs((K.exp(y_true) - K.exp(y_pred)) / K.clip(K.exp(y_true),
                                            K.epsilon(),None))
    return K.mean(diff, axis=-1)

class MyProgressBar():
    def __init__(self):
        self.pbar = None

    def __call__(self, block_num, block_size, total_size):
        if not self.pbar:
            self.pbar=tqdm.trange(total_size)
            self.pbar.set_description("Downloading")
        downloaded = block_num * block_size
        if downloaded < total_size:
            self.pbar.update(downloaded)
        else:
            self.pbar.close()
