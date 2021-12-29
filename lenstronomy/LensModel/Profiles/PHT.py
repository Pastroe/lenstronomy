__author__ = 'not sibirrer'


import numpy as np
from lenstronomy.LensModel.Profiles.base_profile import LensProfileBase

__all__ = ['PHT']


class PHT(LensProfileBase):
    """
    class to compute the Boxy-Disky Correction.
    """
    param_names = ['r_in', 'r_scale', 'theta_out', 'alpha', 'center_x', 'center_y']
    lower_limit_default = {'r_in':0, 'r_scale':0, 'theta_out': -4*np.pi, 'alpha':-10, 'center_x': -100, 'center_y': -100}
    lower_limit_default = {'r_in':100, 'r_scale':100, 'theta_out': -4*np.pi, 'alpha':10, 'center_x': -100, 'center_y': -100}

    def __init__(self):
        self.r_min = 10**(-25)
        super(PHT, self).__init__()
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


    def derivatives(self, x, y, r_in, r_scale, theta_out, alpha, center_x=0, center_y=0):
        """

        :param x: x-coord (in angles)
        :param y: y-coord (in angles)
        :param C_0: Einstein radius (in angles)
        :return: deflection angle (in angles)
        """
        #print('Referred deri')

        x_ = x - center_x
        y_ = y - center_y
        r_ = ((x_ ** 2) + (y ** 2)) ** (1/2)
        theta_rot = theta_out * np.tanh((r_ - r_in) / r_scale) * (0.5 * (1 + r_ / (r_in + r_scale)))**alpha
        x_rot = x_ * np.cos(theta_rot) + y_ * np.sin(theta_rot)
        y_rot =-x_ * np.sin(theta_rot) + y_ * np.cos(theta_rot)
        return x_ - x_rot, y_ - y_rot


    def hessian(self, x, y, r_in, r_scale, theta_out, alpha, center_x=0, center_y=0):
        """

        :param x: x-coord (in angles)
        :param y: y-coord (in angles)
        :param Strength: Einstein radius (in angles)
        :return: hessian matrix (in angles)
        """
        
        print('Referred Hess')
        # Have no choice but use numerical solution
        
        dx = 1E-4
        dy = 1E-4
        df_x1, df_y1 = self.derivatives(x+dx, y, r_in, r_scale, theta_out, alpha, center_x, center_y)
        df_x2, df_y2 = self.derivatives(x-dx, y, r_in, r_scale, theta_out, alpha, center_x, center_y)
        df_x3, df_y3 = self.derivatives(x, y+dy, r_in, r_scale, theta_out, alpha, center_x, center_y)
        df_x4, df_y4 = self.derivatives(x, y-dy, r_in, r_scale, theta_out, alpha, center_x, center_y)

        f_xx = (df_x1 - df_x2) / 2 / dx
        f_yy = (df_y3 - df_y4) / 2 / dy
        f_xy = (df_y1 - df_y2) / 4 / dx + (df_x3 - df_x4) / 4 / dy
        return f_xx, f_xy, f_xy, f_yy