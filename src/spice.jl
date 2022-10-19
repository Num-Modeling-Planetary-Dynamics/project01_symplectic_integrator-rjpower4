# ========================================================================================
# File: spice.jl
# Brief: wrappers around spice calls
# Author: Rolfe Power <rpower@purdue.edu>
# ========================================================================================

"""
    SPICE_KERNEL_EXTENSIONS

Tuple of extensions assoicated with spice kernels
"""
const SPICE_KERNEL_EXTENSIONS = (".bsp", ".tpc", ".tls")

"""
    is_kernel_path(path; require_exists=true)

Determine if the path given points to a valid spice kernel based on path only.

If `require_exists` is `false`, then only extension of file is checked. If it is true,
then also return `false` if the file does not exist.
"""
function is_kernel_path(path::String; require_exists = true)
    _, ext = splitext(path)
    return (isfile(path) || !require_exists) && (ext in SPICE_KERNEL_EXTENSIONS)
end

"""
    list_kernels(path::String)

Return a vector of paths to all kernels in directory.
"""
function list_kernels(path::String)
    isdir(path) || ArgumentError("path not a directory: " + path)
    (r, _, files) = walkdir(normpath(path)) |> first
    root = abspath(r)
    kernels = String[]
    for p in files
        pth = joinpath(root, p)
        if is_kernel_path(pth)
            push!(kernels, pth)
        end
    end
    return kernels
end

"""
    load_kernel(pth::String)

Load the kernel at the given path throwing error if not found.
"""
function load_kernel(pth::String)
    isfile(path) || ArgumentError("kernel file not found: " + pth)
    furnsh(pth)
    return nothing
end

"""
    load_kernels(paths::AbstractVector)

Load all kernels with path in `paths`. Return error and don't load any if any not found.
"""
function load_kernels(pths::AbstractVector)
    if any(x -> !isfile(x), pths)
        error("kernel(s) not found")
    end
    foreach(furnsh, pths)
    return nothing
end

"""
    ephemeris_time(s::String)

Compute the ephemeris time from the time string.
"""
ephemeris_time(s::String) = str2et(s)

"""
    spice_state(id::Integer, et; frame="J2000", observer=0)

Retrieve the state of the body with id `id` at the given ephemeris time.
"""
function spice_state(id::Integer, et; frame = "J2000", observer = 0)
    out, _ = spkez(id, et, frame, "NONE", observer)
    return out
end

"""
    spice_state(name::String, et; frame="J2000", observer="SOLAR_SYSTEM_BARYCENTER")

Retrieve the state of the body with name `name` at the given ephemeris time.
"""
function spice_state(name::String, et; frame = "J2000",
                     observer = "SOLAR_SYSTEM_BARYCENTER")
    out, _ = spkezr(name, et, frame, "NONE", observer)
    return out
end
