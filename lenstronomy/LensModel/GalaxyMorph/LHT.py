__author__ = 'not sibirrer'


import numpy as np
from lenstronomy.LensModel.Profiles.base_profile import LensProfileBase

__all__ = ['LHT']


class LHT(LensProfileBase):
    """
    class to compute the Log-Hyperbolic-Tanh coordinate rotation.
    """
    param_names = ['r_in', 'r_scale', 'theta_out', 'r_wind', 'center_x', 'center_y']
    lower_limit_default = {'r_in':0, 'r_scale':0, 'theta_out': -4*np.pi, 'r_wind':0.01, 'center_x': -100, 'center_y': -100}
    lower_limit_default = {'r_in':100, 'r_scale':100, 'theta_out': -4*np.pi, 'r_wind':100, 'center_x': -100, 'center_y': -100}

    def __init__(self):
        self.r_min = 10**(-25)
        super(LHT, self).__init__()

    def mapping(self, x, y, r_in, r_scale, theta_out, r_wind, center_x=0, center_y=0):
        """

        :param x: x-coord of the image
        :param y: y-coord of the image
        :param r_in: the inner radius of the spiral arm
        :param r_scale: the scale of spiral arm in radius; r_in + r_scale = r_out
        :param theta_out: the total rotation angle of spiral arm
        :param r_wind: wind radius of the spiral arm        
        :return: x, y-coord of the source
        """

        x_ = x - center_x
        y_ = y - center_y
        
        # See the GALFIT II paper for more detail
        CDEF = 0.23
        A = (2 * CDEF) / (np.abs(theta_out) + CDEF) - 1.00001
        B = (2 - np.arctanh(A)) * (r_in + r_scale) / r_scale
        
        r_ = ((x_ ** 2) + (y ** 2)) ** (1/2)
        theta_rot = theta_out * 0.5 * (np.tanh(B * (r_ / (r_in + r_scale) - 1) + 2) + 1) * (np.log(1 + r_ / r_wind) / np.log(1 + (r_in + r_scale) / r_wind))
        
        x_rot = x_ * np.cos(theta_rot) + y_ * np.sin(theta_rot)
        y_rot =-x_ * np.sin(theta_rot) + y_ * np.cos(theta_rot)
        return x_rot + center_x, y_rot + center_y
        
    def reverse_mapping(self, x, y, r_in, r_scale, theta_out, r_wind, center_x=0, center_y=0):
        """

        :param x: x-coord of the source
        :param y: y-coord of the source
        :param r_in: the inner radius of the spiral arm
        :param r_scale: the scale of spiral arm in radius; r_in + r_scale = r_out
        :param theta_out: the total rotation angle of spiral arm
        :param r_wind: wind radius of the spiral arm        
        :return: x, y-coord of the image
        """

        x_ = x - center_x
        y_ = y - center_y
        
        # See the GALFIT II paper for more detail
        CDEF = 0.23
        A = (2 * CDEF) / (np.abs(theta_out) + CDEF) - 1.00001
        B = (2 - np.arctanh(A)) * (r_in + r_scale) / r_scale

        r_ = ((x_ ** 2) + (y ** 2)) ** (1/2)
        theta_rot = theta_out * 0.5 * (np.tanh(B * (r_ / (r_in + r_scale) - 1) + 2) + 1) * (np.log(1 + r_ / r_wind) / np.log(1 + (r_in + r_scale) / r_wind))
        
        x_rot = x_ * np.cos(theta_rot) - y_ * np.sin(theta_rot)
        y_rot = x_ * np.sin(theta_rot) + y_ * np.cos(theta_rot)
        return x_rot + center_x, y_rot + center_y
        
    def param_dict(self, index):
        p_dict = {}
        for param in self.param_names:
            if param not in ['center_x', 'center_y']:
                p_dict[param] = param  + '_' + str(index)
        return p_dict