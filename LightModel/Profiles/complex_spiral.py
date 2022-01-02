import numpy as np
import copy
from lenstronomy.LensModel.galaxy_morph import GalaxyMorph
from lenstronomy.LightModel.Profiles.sersic import SersicElliptic


__all__ = ['ComplexSpiral']

class ComplexSpiral(SersicElliptic):
    """
    class of the edge-on disk (van der Kruit & Searle 1981)

    """
    param_names_source = ['amp', 'R_sersic', 'n_sersic', 'e1', 'e2', 'center_x', 'center_y']
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
            print(self.param_example_dict)

    def ray_shooting(self, x, y, center_x = 0, center_y = 0, **kwargs):
        """

        :param x: ra-coordinate
        :param y: dec-coordinate
        :param kwargs: other parameters, depends on the model
        :return: flux of complex spiral profile
        """

        lens_kwargs_list = [{param:kwargs[dict_map[param]] for j, param in enumerate(dict_map)} for i, dict_map in enumerate(self.param_dict_maps)]
        x_, y_ = self.lens_model.ray_shooting(x - center_x, y - center_y, lens_kwargs_list)
        return x_ + center_x, y_ + center_y
        
    def function(self, x, y, amp, R_sersic, n_sersic, e1, e2, center_x=0, center_y=0, max_R_frac=100.0, **kwargs):
        
        R_sersic = np.maximum(0, R_sersic)
        x_, y_ = self.ray_shooting(x, y, center_x, center_y, **kwargs)
        R = self.get_distance_from_center(x_, y_, e1, e2, center_x, center_y)
        result = self._r_sersic(R, R_sersic, n_sersic, max_R_frac)
        return amp * result
        

        
    
        
        
        
