__author__ = 'not_sibirrer'

import numpy as np
from lenstronomy.LensModel.profile_list_base import ProfileListBase

__all__ = ['GalaxyMorph']


class GalaxyMorph(ProfileListBase):
    """
    class to handle an arbitrary list of lens models in a single lensing plane
    """

    def ray_shooting(self, x, y, kwargs, k=None):
        """
        maps image to source position (inverse deflection)
        :param x: x-position (preferentially arcsec)
        :type x: numpy array
        :param y: y-position (preferentially arcsec)
        :type y: numpy array
        :param kwargs: list of keyword arguments of lens model parameters matching the lens model classes
        :param k: only evaluate the k-th lens model
        :return: source plane positions corresponding to (x, y) in the image plane
        """
        x = np.array(x, dtype=float)
        y = np.array(y, dtype=float)
        if isinstance(k, int):
            return self.func_list[k].mapping(x, y, **kwargs[k])
        bool_list = self._bool_list(k)
        x_, y_ = x, y
        for i, func in enumerate(self.func_list):
            if bool_list[i] is True:
                x_, y_ = func.mapping(x_, y_, **kwargs[i])

        return x_, y_
        
    def reverse_ray_shooting(self, x, y, kwargs, k=None):
        """
        maps source to image position (deflection)
        :param x: x-position (preferentially arcsec)
        :type x: numpy array
        :param y: y-position (preferentially arcsec)
        :type y: numpy array
        :param kwargs: list of keyword arguments of lens model parameters matching the lens model classes
        :param k: only evaluate the k-th lens model
        :return: source plane positions corresponding to (x, y) in the image plane
        """
        x = np.array(x, dtype=float)
        y = np.array(y, dtype=float)
        if isinstance(k, int):
            return self.func_list[k].reverse_mapping(x, y, **kwargs[k])
        bool_list = self._bool_list(k)
        x_, y_ = x, y
        for i, func in enumerate(self.func_list[::-1]):
            if bool_list[-i+1] is True:
                x_, y_ = func.reverse_mapping(x_, y_, **kwargs[-i+1])

        return x_, y_

        
    def polar_chessboard(self, x, y, spatial_freq_r=1, spatial_freq_phi=np.pi/6):
        r = (x**2 + y**2)**0.5
        phi = np.arctan2(y, x)
        return (r % spatial_freq_r > spatial_freq_r / 2) ^ (phi % spatial_freq_phi > spatial_freq_phi / 2)
    
    def cart_chessboard(self, x, y, spatial_freq_x = 1, spatial_freq_y = 1):
        return (x % spatial_freq_x > spatial_freq_x / 2) ^ (y % spatial_freq_y > spatial_freq_y / 2)