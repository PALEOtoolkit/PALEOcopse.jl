# Setup PALEOcopse/examples environment to use local downloaded COPSE package code from PALEOcopse/scalar
# This only needs to be run once, after cloning the github repository

import Pkg

Pkg.activate(".") # use the PALEOcopse/examples environment
Pkg.develop(path="../")   # use the local version of PALEOcopse packages to allow local modifications
Pkg.instantiate()