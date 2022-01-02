import numpy as np
import scipy.special as ss

__all__ = ['EdgeOnDisk']

class EdgeOnDisk(object):
    """
    class of the edge-on disk (van der Kruit & Searle 1981)

    """
    param_names = ['amp', 'r_s', 'h_s', 'PA', 'center_x', 'center_y']
    lower_limit_default = {'amp': 0, 'r_s': 0, 'h_s': 0, 'PA':0, 'center_x': -100, 'center_y': -100}
    upper_limit_default = {'amp': 100, 'r_s': 100, 'h_s': 100, 'PA':np.pi, 'center_x': 100, 'center_y': 100}

    def __init__(self):
        pass

    def function(self, x, y, amp, r_s, h_s, PA, center_x=0, center_y=0):
        """

        :param x: ra-coordinate
        :param y: dec-coordinate
        :param w_c:
        :param w_t:
        :param amp: amplitude of first power-law flux
        :param e1: eccentricity parameter
        :param e2: eccentricity parameter
        :param center_x: center
        :param center_y: center
        :return: flux of chameleon profile
        """
        x_, y_ = x - center_x, y - center_y
        r =  np.abs(np.cos(PA) * x_ + np.sin(PA) * y_)
        h =  np.abs(-np.sin(PA) * x_ + np.cos(PA) * y_)

        flux = amp * (r / r_s) * ss.kv(1, r / r_s) / np.cosh(h / h_s) ** 2
        return flux