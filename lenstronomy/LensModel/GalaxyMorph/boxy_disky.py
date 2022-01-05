__author__ = 'not sibirrer'


import numpy as np
from lenstronomy.LensModel.Profiles.base_profile import LensProfileBase

__all__ = ['BoxyDisky']


class BoxyDisky(LensProfileBase):
    """
    class to compute the Boxy-Disky Correction.
    """
    param_names = ['C_0', 'center_x', 'center_y']
    lower_limit_default = {'C_0': -5, 'center_x': -100, 'center_y': -100}
    upper_limit_default = {'C_0':  5, 'center_x': 100, 'center_y': 100}

    def __init__(self):
        self.r_min = 10**(-25)
        super(BoxyDisky, self).__init__()

    def mapping(self, x, y, C_0, center_x=0, center_y=0):
        """

        :param x: x-coord of the image
        :param y: y-coord of the image
        :param C_0: index account for the boxy-disky shape
        :return: x, y-coord of the source
        """

        x_ = np.abs(x - center_x)
        y_ = np.abs(y - center_y)
        r_ = (x_**2 + y_**2) ** (1/2)
        r_BD = (x_**(C_0+2) + y_**(C_0+2)) ** (1/(C_0+2))
        
        return (x - center_x)*(r_BD/r_) + center_x, (y - center_y)*(r_BD/r_) + center_y
        
    def reverse_mapping(self, x, y, C_0, center_x=0, center_y=0):
        """

        :param x: x-coord of the source
        :param y: y-coord of the source
        :param C_0: index account for the boxy-disky shape
        :return: x, y-coord of the image
        """

        x_ = np.abs(x - center_x)
        y_ = np.abs(y - center_y)
        r_ = (x_**2 + y_**2) ** (1/2)
        r_BD = (x_**(C_0+2) + y_**(C_0+2)) ** (1/(C_0+2))
        
        return (x - center_x)*(r_/r_BD) + center_x, (y - center_y)*(r_/r_BD) + center_y
        
    def param_dict(self, index):
        p_dict = {}
        for param in self.param_names:
            if param not in ['center_x', 'center_y']:
                p_dict[param] = param  + '_' + str(index)
        return p_dict