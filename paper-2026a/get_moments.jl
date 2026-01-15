using FourierFlows, CUDA, ProgressMeter, JLD2, HDF5

include("../core/circulation.jl")
include("../core/statistics.jl")

grid = TwoDGrid(GPU(), nx = 1024, Lx = 2π)

fid = h5open((@__DIR__)*"/vorticity_hat.h5", "r")
ζhs = Iterators.map(1:1024) do n
	device_array(grid)(fid["data"][:, :, n])
end

Γ_moment(hsh, k) = begin
	sum(ζhs) do ζh
		Γ = getΓ(ζh, hsh, grid)
		moment(abs.(Γ), k)
	end / length(ζhs)
end

hshs = [L_shapedloop(grid, 0.1rand(), 0.1rand(), 0.1rand(), 0.1rand()) for _ in 1:200]

orders = -0.5:0.5:4

moments = map(orders) do k
	res = zeros(length(hshs))
	@showprogress for (n, hsh) in enumerate(hshs)
		res[n] = Γ_moment(hsh, k)
	end
	res
end
