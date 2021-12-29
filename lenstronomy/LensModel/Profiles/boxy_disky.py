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
        # alpha = 4*const.G * (mass*const.M_sun)/const.c**2/(r*const.Mpc)

    def function(self, x, y, C_0, center_x=0, center_y=0):
        """

        :param x: x-coord (in angles)
        :param y: y-coord (in angles)
        :param C_0: Einstein radius (in angles)
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
        phi = C_0**2*np.log(r)
        return phi


    def derivatives(self, x, y, C_0, center_x=0, center_y=0):
        """

        :param x: x-coord (in angles)
        :param y: y-coord (in angles)
        :param C_0: Einstein radius (in angles)
        :return: deflection angle (in angles)
        """
        print('Referred deri')

        x_ = np.abs(x - center_x)
        y_ = np.abs(y - center_y)
        r_ = (x_**2 + y_**2) ** (1/2)
        r_BD = (x_**(C_0+2) + y_**(C_0+2)) ** (1/(C_0+2))
        
        return -(x - center_x)*(r_BD/r_-1), -(y - center_y)*(r_BD/r_-1)


    def hessian(self, x, y, C_0, center_x=0, center_y=0):
        """

        :param x: x-coord (in angles)
        :param y: y-coord (in angles)
        :param C_0: Einstein radius (in angles)
        :return: hessian matrix (in angles)
        """
        print('Referred Hess')

        x_ = x - center_x
        y_ = y - center_y
        C = C_0**2
        a = x_**2 + y_**2
        if isinstance(a, int) or isinstance(a, float):
            r2 = max(self.r_min**2, a)
        else:
            r2 = np.empty_like(a)
            r2[a > self.r_min**2] = a[a > self.r_min**2]  #in the SIS regime
            r2[a <= self.r_min**2] = self.r_min**2
        f_xx = C * (y_**2-x_**2)/r2**2
        f_yy = C * (x_**2-y_**2)/r2**2
        f_xy = -C * 2*x_*y_/r2**2
        return f_xx, f_xy, f_xy, f_yy