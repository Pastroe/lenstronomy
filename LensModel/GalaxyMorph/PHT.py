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

    def mapping(self, x, y, r_in, r_scale, theta_out, alpha, center_x=0, center_y=0):
        """

        :param x: x-coord (in angles)
        :param y: y-coord (in angles)
        :param C_0: Einstein radius (in angles)
        :return: deflection angle (in angles)
        """
        #print('Referred deri')
        CDEF = 0.23
        A = (2 * CDEF) / (np.abs(theta_out) + CDEF) - 1.00001
        B = (2 - np.arctanh(A)) * (r_in + r_scale) / r_scale
        x_ = x - center_x
        y_ = y - center_y
        r_ = ((x_ ** 2) + (y ** 2)) ** (1/2)
        theta_rot = theta_out * 0.5 * (np.tanh(B * (r_ / (r_in + r_scale) - 1) + 2) + 1) * (0.5 * (1 + r_ / (r_in + r_scale)))**alpha
        x_rot = x_ * np.cos(theta_rot) + y_ * np.sin(theta_rot)
        y_rot =-x_ * np.sin(theta_rot) + y_ * np.cos(theta_rot)
        return x_rot + center_x, y_rot + center_y
        
    def param_dict(self, index):
        p_dict = {}
        for param in self.param_names:
            if param not in ['center_x', 'center_y']:
                p_dict[param] = param  + '_' + str(index)
        return p_dict