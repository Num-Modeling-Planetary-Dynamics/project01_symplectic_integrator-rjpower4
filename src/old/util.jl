# ========================================================================================
# File: util.jl
# Brief: utility functions
# Author: Rolfe Power <rpower@purdue.edu>
# ========================================================================================

"""
    angular_difference(a, b)

Return the correctly signed counter-clockwise angular difference from a to b.
"""
function angular_difference(a, b)
    ay, ax = sincos(a)
    by, bx = sincos(b)
    theta = atan(ax * by - ay * bx, ax * bx + ay * by)
    return theta > 0 ? theta : 2π + theta
end
