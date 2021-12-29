__author__ = 'not sibirrer'


import numpy as np
from lenstronomy.LensModel.Profiles.base_profile import LensProfileBase

__all__ = ['Fourier']


class Fourier(LensProfileBase):
    """
    class to compute the Boxy-Disky Correction.
    """
    param_names = ['m_1', 'm_3', 'm_4', 'm_5', 'm_6', 'phi_1', 'phi_3', 'phi_4', 'phi_5', 'phi_6', 'center_x', 'center_y']
    lower_limit_default = {'m_1':-1, 'm_3':-1, 'm_4':-1, 'm_5':-1, 'm_6':-1,'phi_1':0, 'phi_3':0, 'phi_4':0, 'phi_5':0, 'phi_6':0, 'center_x': -100, 'center_y': -100}
    upper_limit_default = {'m_1':1, 'm_3':1, 'm_4':1, 'm_5':1, 'm_6':1,'phi_1':np.pi, 'phi_3':np.pi, 'phi_4':np.pi, 'phi_5':np.pi, 'phi_6':np.pi, 'center_x': -100, 'center_y': -100}

    def __init__(self):
        self.r_min = 10**(-25)
        super(Fourier, self).__init__()
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


    def derivatives(self, x, y, m_1, m_3, m_4, m_5, m_6, phi_1, phi_3, phi_4, phi_5, phi_6, center_x=0, center_y=0):
        """

        :param x: x-coord (in angles)
        :param y: y-coord (in angles)
        :param C_0: Einstein radius (in angles)
        :return: deflection angle (in angles)
        """
        print('Referred deri')

        x_ = x - center_x
        y_ = y - center_y
        theta_ = np.arctan2(y_, x_)
        f_mode_list = [1,3,4,5,6]
        f_coef_list = [m_1, m_3, m_4, m_5, m_6]
        f_angl_list = [phi_1, phi_3, phi_4, phi_5, phi_6]
        F_scale = np.zeros_like(x)
        for i, m in enumerate(f_mode_list):
            F_scale = F_scale + f_coef_list[i] * np.cos(m * (theta_ + f_angl_list[i]))
        return -(x - center_x) * F_scale, -(y - center_y) * F_scale

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