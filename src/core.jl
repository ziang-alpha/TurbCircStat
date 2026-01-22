using FourierFlows, HDF5, JLD2, Measurements, CUDA
include("circulation.jl")
include("diagnostics.jl")
include("statistics.jl")


"""
    Estimate the mean value and standard deviation of a series of data.
    Each datum frame could be an Array or CuArray.
    The results is given by 'mean .± sqrt.(var)'.
    Here, '±' is provided by the Measurements.jl.
"""
function measure(data_series)
    ndata = length(data_series)
    mean = sum(data_series) / ndata
    var = sum(datum -> (datum - mean) .^ 2, data_series) / ndata^2
    return mean .± sqrt.(var)
end
