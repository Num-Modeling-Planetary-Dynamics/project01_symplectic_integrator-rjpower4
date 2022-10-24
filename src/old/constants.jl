# ========================================================================================
# File: dimensions.jl
# Brief: definition of useful physical constants
# Author: Rolfe Power <rpower@purdue.edu>
#
# NOTE: All dimensions defined in the Kilometer-Kilogram-Second system
# ========================================================================================

"""
GRAVITATIONAL_CONSTANT

Universal Constant of Gravitation

This is also known as the Newtonian Constant of Gravitation and the Cavendish
Gravitational Constant. In the classical (Newtonian) formulation of gravity,
this value serves as the constant of proportionality between the force
exterted between two bodies and the sum of their products divided by the
square of the distance between them.

Note that the value of this constant is only really known to about four
significant figures. Therefore, if its use can be avoided it should be.

* Dimensions : [M^(-1) L^3 T^(-2)]
* Units      : [km^3 / (kg . s^2)]
* Sources    : [Nist](https://physics.nist.gov/cuu/Constants/)
               [ArXiv](https://arxiv.org/abs/1507.07956)
"""
const GRAVITATIONAL_CONSTANT = 6.6743e-20

"""
    ASTRONOMICAL_UNIT

Astronomical unit

Unit of length that is about equal to the distance from the Earth to the Sun.
According to WikiPedia, the distance from the Earth to the Sun varies by about
3% over the course of a year so this value could be (in theory) about 3% off
the distance for any given epoch. However, the use case of this unit is that
it is a convenient scaling factor for distances of solar system magnitudes as
it is officially defined.

* Dimensions : [L]
* Units      : [km]
* Sources    : [IAU](https://www.iau.org/static/resolutions/IAU2012_English.pdf)
"""
const ASTRONOMICAL_UNIT = 149597870.700

"""
    SECONDS_PER_DAY

Seconds in a Julian Day

Number of seconds in a Julian day. It's what it says on the tin...

* Dimensions : [T]
* Units      : [s]
* Sources    : SPICE
"""
const SECONDS_PER_DAY = 86400.0

"""
    SECONDS_PER_YEAR

Seconds in a Julian Year

Number of seconds in a Julian year. Did you also get confused about the
SECONDS_PER_DAY one...?

* Dimensions : [T]
* Units      : [s]
* Sources    : SPICE
"""
const SECONDS_PER_YEAR = 31557600.0
