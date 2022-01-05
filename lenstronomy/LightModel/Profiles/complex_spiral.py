import numpy as np
import copy
from lenstronomy.LensModel.galaxy_morph import GalaxyMorph
from lenstronomy.LightModel.Profiles.sersic import SersicElliptic


__all__ = ['ComplexSpiral']

class ComplexSpiral(SersicElliptic):
    """
    class of the edge-on disk (van der Kruit & Searle 1981)

    """
    param_names_source = ['amp', 'R_sersic', 'n_sersic', 'e1', 'e2', 'theta_PA', 'q', 'center_x', 'center_y']
    _s = 0.00001


    def __init__(self, model_list, smoothing=_s, sersic_major_axis=False, is_print_eg_param=True):

        self._smoothing = smoothing
        self._sersic_major_axis = sersic_major_axis

        self.lens_model = GalaxyMorph(model_list)
        self.param_dict_maps = [func.param_dict(i) for i, func in enumerate(self.lens_model.func_list)]
        self.param_example_dict = {param_name:0 for param_name in self.param_names_source}
        for param_dict_map in self.param_dict_maps:
            self.param_example_dict.update({value:0 for value in param_dict_map.values()})
        if is_print_eg_param:
            print("Example of parameter list:\n", self.param_example_dict)

    def ray_shooting(self, x, y, theta_PA, q, center_x = 0, center_y = 0, **kwargs):
        """

        :param x: ra-coordinate
        :param y: dec-coordinate
        :param kwargs: other parameters, depends on the model
        :return: flux of complex spiral profile
        """
        # x, y are coord on the image
        # remove the offset
        x_, y_ = x - center_x, y - center_y
        # rotate to the standrd position (x-major axis, y-minor axis)
        x_rot = x_ * np.cos(theta_PA) + y_ * np.sin(theta_PA)
        y_rot =-x_ * np.sin(theta_PA) + y_ * np.cos(theta_PA)
        # deprojection
        x_ = x_rot
        y_ = y_rot * q

        lens_kwargs_list = [{param:kwargs[dict_map[param]] for j, param in enumerate(dict_map)} for i, dict_map in enumerate(self.param_dict_maps)]
        x_, y_ = self.lens_model.ray_shooting(x_, y_, lens_kwargs_list)
        return x_ + center_x, y_ + center_y
        
    def reverse_ray_shooting(self, x, y, theta_PA, q, center_x = 0, center_y = 0, **kwargs):
        """

        :param x: ra-coordinate
        :param y: dec-coordinate
        :param kwargs: other parameters, depends on the model
        :return: flux of complex spiral profile
        """
        
        # x, y are coord on the source
        # remove the offset
        x_, y_ = x - center_x, y - center_y
        
        lens_kwargs_list = [{param:kwargs[dict_map[param]] for j, param in enumerate(dict_map)} for i, dict_map in enumerate(self.param_dict_maps)]
        x_, y_ = self.lens_model.reverse_ray_shooting(x_, y_, lens_kwargs_list)
        
        # deprojection
        x_ = x_
        y_ = y_ / q
        
        # rotate to the standrd position (x-major axis, y-minor axis)
        x_rot = x_ * np.cos(theta_PA) - y_ * np.sin(theta_PA)
        y_rot = x_ * np.sin(theta_PA) + y_ * np.cos(theta_PA)
        

        return x_rot + center_x, y_rot + center_y

    def function(self, x, y, amp, R_sersic, n_sersic, e1, e2, theta_PA = 0, q = 1, center_x=0, center_y=0, max_R_frac=100.0, **kwargs):

        R_sersic = np.maximum(0, R_sersic)
        x_, y_ = self.ray_shooting(x, y, theta_PA, q, center_x, center_y, **kwargs)
        R = self.get_distance_from_center(x_, y_, e1, e2, center_x, center_y)
        result = self._r_sersic(R, R_sersic, n_sersic, max_R_frac)
        return amp * result
        
    def mapping_image_to_source(self, theta_PA, q, coord = 'cart', numPix = 500, deltaPix = 0.05, center_x = 0, center_y = 0, **kwargs):
        '''
        
        Mapping a grid of the source model into the image model.
        :param coord: can be 'cart' or 'polar'. This parameter defines the pattern of the grid.
        '''
        Scale = numPix * deltaPix
        x, y = np.meshgrid(np.arange(center_x-Scale/2, center_x+Scale/2, deltaPix), 
                           np.arange(center_y-Scale/2, center_y+Scale/2, deltaPix))
        x_, y_ = self.ray_shooting(x, y, theta_PA, q, **kwargs)
        
        assert coord in ['cart', 'polar'], 'param coord must be one of "cart" or "polar".'
        if coord is 'cart':
            origin = self.lens_model.cart_chessboard(x,y)
            mapping = self.lens_model.cart_chessboard(x_,y_)
        elif coord is 'polar':
            origin = self.lens_model.polar_chessboard(x,y)
            mapping = self.lens_model.polar_chessboard(x_,y_)

        return origin, mapping
    
        
    def mapping_source_to_image(self, theta_PA, q, coord = 'cart', numPix = 500, deltaPix = 0.05, center_x = 0, center_y = 0, **kwargs):
        '''
        
        Mapping a grid of the source model into the image model.
        :param coord: can be 'cart' or 'polar'. This parameter defines the pattern of the grid.
        '''
        Scale = numPix * deltaPix
        x, y = np.meshgrid(np.arange(center_x-Scale/2, center_x+Scale/2, deltaPix), 
                           np.arange(center_y-Scale/2, center_y+Scale/2, deltaPix))
                           
        x_, y_ = self.reverse_ray_shooting(x, y, theta_PA, q, **kwargs)
        
        assert coord in ['cart', 'polar'], 'param coord must be one of "cart" or "polar".'
        if coord is 'cart':
            origin = self.lens_model.cart_chessboard(x,y)
            mapping = self.lens_model.cart_chessboard(x_,y_)
        elif coord is 'polar':
            origin = self.lens_model.polar_chessboard(x,y)
            mapping = self.lens_model.polar_chessboard(x_,y_)

        return origin, mapping
    







