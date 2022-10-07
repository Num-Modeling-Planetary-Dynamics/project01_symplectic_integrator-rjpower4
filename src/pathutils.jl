# ========================================================================================
# File: pathutils.jl
# Author: Rolfe Power <rpower@purdue.edu>
# Brief: Functions and types for dealing with paths on disk
# ========================================================================================

"""
    project_directory(subpaths...)

Return path to the root directory of the EAPS project optionally appending path segments.
"""
project_directory(args...) = joinpath(dirname(@__DIR__), args...)

"""
    data_directory(subpaths...)

Return path to the data directory of the EAPS project optionally appending path segments.
"""
data_directory(args...) = project_directory("data", args...)

"""
    plot_directory(subpaths...)

Return path to the plot directory of the EAPS project optionally appending path segments.
"""
plot_directory(args...) = project_directory("plots", args...)