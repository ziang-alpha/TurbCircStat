using HDF5, ArgParse, Measurements
include((@__DIR__) * "/../src/circulation.jl")
include((@__DIR__) * "/../src/diagnostics.jl")
include((@__DIR__) * "/../src/statistics.jl")

aps = ArgParseSettings()
@add_arg_table! aps begin
	"--path", "-p"
	help = "Path of the raw data."
	arg_type = String
	default = (@__DIR__)
	"--mode", "-m"
	help = "Mode of the output .h5 file."
	arg_type = String
	default = "w"
    "--gpu"
	help = "Using GPU acceleration."
	action = :store_true
end
args = parse_args(aps)

dev = args["gpu"] ? GPU() : CPU()
filenames = readdir(args["path"] * "/.output/")

rfiles = map(filenames) do name
	h5open(args["path"] * "/.output/" * name, "r")
end

rawdata = map(rfiles, filenames) do fid, name
	ζhs = begin
		nframe = maximum(parse.(Int, keys(fid)))
		Iterators.map(1:nframe) do n
			ζh = read(fid, "$(n)")
			device_array(dev)(ζh)
		end
	end
	regex = r"Re(\d+\.\d+)_N(\d+)\.h5"
	m = match(regex, name)
	Re, N = parse(Float64, m[1]), parse(Int, m[2])
	grid = TwoDGrid(dev; nx = N, Lx = 2π)
	fh = device_array(dev)(float(grid.Krsq .< (N / 3sqrt(Re) - 5)^2))
	fζhs = Iterators.map(ζhs) do ζh
		ζh .* fh
	end

	Re => (ζhs = ζhs, fζhs = fζhs, grid = grid)
end
rawdata = Dict(rawdata)

pfile = h5open(args["path"] * "/postdata.h5", args["mode"])

include("datalist.jl")

close(pfile)
close.(rfiles)


