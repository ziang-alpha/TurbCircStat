using CUDA, Measurements
# Compute the probability distribution funtion for samples given by an Array
function pdf(sp::Array, bin)
    pdf = map(bin[1:(end-1)], bin[2:end]) do l, r
        count(x -> l ≤ x < r, sp) / ((r - l) * length(sp))
    end

    return pdf
end

# Compute the probability distribution funtion for samples given by a CuArray
function pdf(sp::CuArray, bin)
    min_val, bin_width, nbins = bin[1], bin[2] - bin[1], length(bin) - 1

    hist = CUDA.zeros(nbins)

    threads_per_block = 256

    blocks_per_grid = cld(length(sp), threads_per_block)

    @cuda threads = threads_per_block blocks = blocks_per_grid histogram_kernel!(
        hist, sp, min_val, bin_width, nbins,
    )

    pdf = hist / (bin_width * length(sp))

    return pdf
end

function histogram_kernel!(hist, data, min_val, bin_width, nbins)
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if idx <= length(data)
        value = data[idx]
        if value >= min_val && value <= min_val + bin_width * nbins
            bin_idx = Int(floor((value - min_val) / bin_width)) + 1
            bin_idx = max(1, min(bin_idx, nbins))
            CUDA.@atomic hist[bin_idx] += 1
        end
    end
    return nothing
end

moment(sp, k) = sum(x -> x^k, sp) / length(sp)

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