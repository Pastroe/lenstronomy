__author__ = 'not sibirrer'


import numpy as np
from lenstronomy.LensModel.Profiles.base_profile import LensProfileBase

__all__ = ['Bending3']


class Bending3(LensProfileBase):
    """
    class to compute the Bending Correction with power-law (alpha = 3).
    """
    param_names = ['Strength', 'center_x', 'center_y']
    lower_limit_default = {'Strength': -100, 'center_x': -100, 'center_y': -100}
    upper_limit_default = {'Strength':  100, 'center_x': 100, 'center_y': 100}

    def __init__(self):
        self.r_min = 10**(-25)
        super(Bending3, self).__init__()
        # alpha = 4*const.G * (mass*const.M_sun)/const.c**2/(r*const.Mpc)

    def mapping(self, x, y, Strength, center_x=0, center_y=0):
        """

        :param x: x-coord (in angles)
        :param y: y-coord (in angles)
        :param Strength: Einstein radius (in angles)
        :return: deflection angle (in angles)
        """
        x_ = x - center_x
        y_ = y - center_y
        return x_ + center_x, y_ + Strength * x_ ** 3 + center_y
