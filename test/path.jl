# ========================================================================================
# File   : test/path_tests.jl
# Brief  : Tests for path stuff
# Author : Rolfe Power <rpower@purdue.edu>
# ========================================================================================
using Symple: project_dir, data_dir

@test isdir(project_dir())
@test isdir(project_dir("src"))
@test isdir(project_dir("test"))
@test isfile(project_dir("src", "Symple.jl"))
@test isfile(project_dir("test", "runtests.jl"))

@test isdir(data_dir())
