# ========================================================================================
# File   : test/constants.jl
# Brief  : Tests for constants
# Author : Rolfe Power <rpower@purdue.edu>
# ========================================================================================

using Symple: GRAVITATIONAL_CONSTANT

@test GRAVITATIONAL_CONSTANT > 0.0
