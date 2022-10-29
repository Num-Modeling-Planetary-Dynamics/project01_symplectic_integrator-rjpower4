# ========================================================================================
# File   : src/path.jl
# Brief  : Path utilities
# Author : Rolfe Power <rpower@purdue.edu>
# ========================================================================================
"""
    project_dir()
    project_dir(subs...)

Return the path to the Symple project optionally appending sub paths.
"""
project_dir() = dirname(@__DIR__)
project_dir(subs...) = joinpath(project_dir(), subs...)

"""
    data_dir()

Return the path to the Symple data directory optionally appending sub paths.
"""
data_dir(subs...) = project_dir("data", subs...)
