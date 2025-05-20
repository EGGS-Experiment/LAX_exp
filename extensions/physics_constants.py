"""
LAX.extensions.physics_constants

Contains physics-related constants.
"""
__all__ = []


# fundamental constants
hbar =      1.0545718e-34   # reduced Planck Constant
kb =        1.38064852e-23  # boltzmann's constant
c =         2.99792458e8    # speed of light (in vacuum)
qe =        1.60217662e-19  # charge of electron
me =        9.10938356e-31  # electron mass
amu =       1.66053904e-27  # atomic mass unit
epsilon_0 = 8.8541878128-12 # vacuum permittivity
__all__.extend(['hbar', 'kb', 'c', 'qe', 'me', 'amu', 'epsilon_0'])

# derived constants
debye =     1e-21/c         # debye
alpha_fs =  2.30707751e-28 / (hbar * c) # fine structure constant
bohr_magneton = hbar ** 2. / (me * 2.30707751e-28) # bohr magneton
__all__.extend(['debye', 'alpha_fs', 'bohr_magneton'])

# experiment-specific constants
r0 =        550e-6      # trap radius
mCa =       40 * amu    # mass of Calcium
__all__.extend(['r0', 'mCa'])

