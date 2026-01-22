include((@__DIR__) * "/../../core/core.jl")

filenames = readdir((@__DIR__) * "/../.output/")
dev = CPU()

getζhs(fid) = begin
    nframe = maximum(parse.(Int, keys(fid)))
    Iterators.map(1:nframe) do n
        ζh = read(fid, "$(n)")
        device_array(dev)(ζh)
    end
end

raw_data = map(filenames) do name
    fid = h5open((@__DIR__) * "/../.output/" * name, "r")
    ζhs = getζhs(fid)
    regex = r"Re(\d+\.\d+)_N(\d+)\.h5"
    m = match(regex, name)
    Re, N = parse(Float64, m[1]), parse(Int, m[2])
    grid = TwoDGrid(dev; nx=N, Lx=2π)
    fh = float(grid.Krsq .< (N / 3sqrt(Re) - 5)^2)
    (ζhs, grid, fh, fid)
end
