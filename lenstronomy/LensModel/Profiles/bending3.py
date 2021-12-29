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

    def function(self, x, y, Strength, center_x=0, center_y=0):
        """

        :param x: x-coord (in angles)
        :param y: y-coord (in angles)
        :param Strength: Einstein radius (in angles)
        :return: lensing potential
        """
        print('Referred Func')
        x_ = x - center_x
        y_ = y - center_y
        a = np.sqrt(x_**2 + y_**2)
        if isinstance(a, int) or isinstance(a, float):
            r = max(self.r_min, a)
        else:
            r = np.empty_like(a)
            r[a > self.r_min] = a[a > self.r_min]  #in the SIS regime
            r[a <= self.r_min] = self.r_min
        phi = Strength**2*np.log(r)
        return phi


    def derivatives(self, x, y, Strength, center_x=0, center_y=0):
        """

        :param x: x-coord (in angles)
        :param y: y-coord (in angles)
        :param Strength: Einstein radius (in angles)
        :return: deflection angle (in angles)
        """
        #print('Referred deri')
        x_ = x - center_x
        y_ = y - center_y
        return 0, - Strength * x_ ** 3


    def hessian(self, x, y, Strength, center_x=0, center_y=0):
        """

        :param x: x-coord (in angles)
        :param y: y-coord (in angles)
        :param Strength: Einstein radius (in angles)
        :return: hessian matrix (in angles)
        """
        print('Referred Hess')

        f_xx = 0
        f_xy = 0
        f_yy = 0
        
        
        return f_xx, f_xy, f_xy, f_yy