using HDF5, FourierFlows

grid = TwoDGrid(GPU(), nx = 1024, Lx = 2π)

fid = h5open((@__DIR__)*"/vorticity_hat.h5", "w")
create_dataset(fid, "data", ComplexF64, (513, 1024, 1024); chunk = (513, 1024, 1))

for n in 1:1024

	(u, v) = h5open((@__DIR__)*"/.giverny_output/isotropic1024coarse_z$(n).h5", "r") do file
		vecU = read(file, "velocity_0001")
		(device_array(grid)(vecU[1, :, :, 1]), device_array(grid)(vecU[2, :, :, 1]))
	end

	ωh = begin
		uh, vh = grid.rfftplan * u, grid.rfftplan * v
		(grid.kr .* vh - grid.l .* uh)
	end

	ωh[grid.Krsq .> 8*512^2/9] .= 0

	fid["vorticity"][:, :, n] = Array(ωh)
end

close(fid)
