"""
LAX.analysis

Contains all standard analysis modules for the EGGS Experiment at UCLA.
"""
__all__ = []

# processing
from LAX_exp.analysis import processing
from LAX_exp.analysis.processing import *
__all__.extend(processing.__all__)

# fitting
from LAX_exp.analysis import fitting
from LAX_exp.analysis.fitting import *
__all__.extend(fitting.__all__)

# optimization
from LAX_exp.analysis import optimize
from LAX_exp.analysis.optimize import *
__all__.extend(optimize.__all__)
