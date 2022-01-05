__author__ = 'not sibirrer'


import numpy as np
from lenstronomy.LensModel.Profiles.base_profile import LensProfileBase

__all__ = ['Bending3']


class Bending3(LensProfileBase):
    """
    class to compute the Bending Correction with power-law (alpha = 3).
    """
    param_names = ['Strength', 'center_x', 'center_y', 'phi']
    lower_limit_default = {'Strength': -100, 'center_x': -100, 'center_y': -100}
    upper_limit_default = {'Strength':  100, 'center_x': 100, 'center_y': 100}

    def __init__(self):
        self.r_min = 10**(-25)
        super(Bending3, self).__init__()

    def mapping(self, x, y, Strength, phi = 0, center_x=0, center_y=0):
        """

        :param x: x-coord of the image
        :param y: y-coord of the image
        :param Strength: strength of bending (y* = y + strength * x**3)
        :return: x, y-coord of the source
        """
        x_ = x - center_x
        y_ = y - center_y
        
        x_rot = x_ * np.cos(phi) + y_ * np.sin(phi)
        y_rot = -x_ * np.sin(phi) + y_ * np.cos(phi)
        
        x_ = x_rot
        y_ = y_rot + Strength * x_rot ** 3
        return x_ + center_x, y_ + center_y
        
    def reverse_mapping(self, x, y, Strength, center_x=0, center_y=0):
        """

        :param x: x-coord of source
        :param y: y-coord of source
        :param Strength: strength of bending (y* = y + strength * x**3)
        :return: x, y-coord of the image
        """
        x_ = x - center_x
        y_ = y - center_y
        
        x_ = x_
        y_ = y_ - Strength * x_ ** 3
        
        x_rot = x_ * np.cos(phi) - y_ * np.sin(phi)
        y_rot = x_ * np.sin(phi) + y_ * np.cos(phi)
        return x_rot + center_x, y_rot + center_y
        
    def param_dict(self, index):
        p_dict = {}
        for param in self.param_names:
            if param not in ['center_x', 'center_y']:
                p_dict[param] = param  + '_' + str(index)
        return p_dict
        
