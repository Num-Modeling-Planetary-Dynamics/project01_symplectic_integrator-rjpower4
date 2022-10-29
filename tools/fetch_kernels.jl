# ========================================================================================
# File   : tools/fetch_kernels.jl
# Brief  : Utility script to download the kernels
# Author : Rolfe Power <rpower@purdue.edu>
# ========================================================================================
using Symple: data_dir
using Logging
using Downloads

# ----------------------------------------------------------------------------------------
# Data Set Type Definition
# ----------------------------------------------------------------------------------------
"""
    DataSet

Representation of a fetchable dataset.
"""
struct DataSet
    filename::String
    url::String
    md5_hash::String
end

"""
    filename

Retrieve the filename of the data set.
"""
filename(ds::DataSet) = ds.filename

"""
    url

Retrieve the url of the data set.
"""
url(ds::DataSet) = ds.url

"""
    md5_hash

Retrieve the md5 hash of the data set.
"""
md5_hash(ds::DataSet) = ds.md5_hash

# ----------------------------------------------------------------------------------------
# Define Kernels
# ----------------------------------------------------------------------------------------
KERNELS = Dict{String,DataSet}()

KERNELS["naif0012"] = DataSet(
    "naif0012.tls",
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls",
    "25a2fff30b0dedb4d76c06727b1895b1",
)

KERNELS["de440"] = DataSet(
    "de440.bsp",
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp",
    "c9d581bfd84209dbeee8b1583939b148",
)

KERNELS["gm_de431"] = DataSet(
    "gm_de431.tpc",
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de431.tpc",
    "6445f12003d1effcb432ea2671a51f18",
)

KERNELS["nep097"] = DataSet(
    "nep097.bsp",
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/nep097.bsp",
    "937e05462aac5e34cb5cfddae6c99199",
)

KERNELS["plu058"] = DataSet(
    "plu058.bsp",
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/plu058.bsp",
    "f1173cd80bd13feb91d0e29487292bbf",
)

# ----------------------------------------------------------------------------------------
# Fetching
# ----------------------------------------------------------------------------------------
function fetch(ds::DataSet, dir::String; force=false)
    output_path = joinpath(dir, filename(ds))

    if isfile(output_path) && !force
        @warn("file exists at '$output_path' and not forcing, skipping")
        return false
    end

    isdir(dir) || mkpath(dir)

    Downloads.download(url(ds), output_path)

    return true
end

# ----------------------------------------------------------------------------------------
# Entry
# ----------------------------------------------------------------------------------------
for (key, dset) in KERNELS
    @info "Fetching $key"
    fetch(dset, data_dir("kernels"))
end
