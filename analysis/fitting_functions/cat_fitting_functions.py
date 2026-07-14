import numpy as np

def _construct_mass_spec_phi(   delta, alpha, t, t0, delta_0, phi):
    delta_shift = 2 * np.pi * (delta - delta_0)
    x = delta_shift * t / 2
    mass_spec_phi = 2 * alpha * np.sinc(x / np.pi) * np.sin(x + delta_shift * t0 - phi)
    return mass_spec_phi

def get_single_ion_cat_lineshape():
    return [_fit_func_single_ion_d, _fit_func_single_ion_b]

def _fit_func_single_ion_b(delta, d, alpha, t, t0, delta_0, phi):
    mass_spec_phi =  _construct_mass_spec_phi(delta, alpha, t, t0, delta_0, phi)
    return (1 - 2 * d) * np.cos(mass_spec_phi) ** 2 + d

def _fit_func_single_ion_d(delta, d, alpha, t, t0, delta_0, phi):
    mass_spec_phi =  _construct_mass_spec_phi(delta, alpha, t, t0, delta_0, phi)
    return (1 - 2 * d) * np.sin(mass_spec_phi) ** 2 + d

def get_unentangled_two_ion_cat_lineshape():
    return [_fit_func_unentangled_two_ion_dd, _fit_func_unentangled_two_ion_bd,
             _fit_func_unentangled_two_ion_bb]

def _fit_func_unentangled_two_ion_bb(delta, d, alpha, t, t0, delta_0, phi):
    mass_spec_phi =  _construct_mass_spec_phi(delta, alpha, t, t0, delta_0, phi)
    return (1 - 2 * d) * np.cos(mass_spec_phi) ** 4 + d

def _fit_func_unentangled_two_ion_dd(delta, d, alpha, t, t0, delta_0, phi):
    mass_spec_phi =  _construct_mass_spec_phi(delta, alpha, t, t0, delta_0, phi)
    return (1 - 2 * d) * np.sin(mass_spec_phi) ** 4 + d

def _fit_func_unentangled_two_ion_bd(delta, d, alpha, t, t0, delta_0, phi):
    mass_spec_phi =  _construct_mass_spec_phi(delta, alpha, t, t0, delta_0, phi)
    fit_bb = (1 - 2 * d) * np.cos(mass_spec_phi) ** 4 + d
    fit_dd = (1 - 2 * d) * np.sin(mass_spec_phi) ** 4 + d
    return 1 - fit_bb - fit_dd

def get_entangled_two_ion_cat_lineshape():
    return [ _fit_func_entangled_two_ion_dd,  _fit_func_entangled_two_ion_bd,
             _fit_func_entangled_two_ion_bb]

def _fit_func_entangled_two_ion_bb(delta, d, alpha, t, t0, delta_0, phi):
    mass_spec_phi =  _construct_mass_spec_phi(delta, alpha, t, t0, delta_0, phi)
    return (1 - 2 * d) * np.cos(mass_spec_phi) ** 2 + d

def _fit_func_entangled_two_ion_dd(delta, d, alpha, t, t0, delta_0, phi):
    mass_spec_phi =  _construct_mass_spec_phi(delta, alpha, t, t0, delta_0, phi)
    return (1 - 2 * d) * np.sin(mass_spec_phi) ** 2 + d

def _fit_func_entangled_two_ion_bd(delta, d, alpha, t, t0, delta_0, phi):
    mass_spec_phi =  _construct_mass_spec_phi(delta, alpha, t, t0, delta_0, phi)
    fit_bb = (1 - 2 * d) * np.cos(mass_spec_phi) ** 2 + d
    fit_dd = (1 - 2 * d) * np.sin(mass_spec_phi) ** 2 + d
    return 1 - fit_bb - fit_dd
