@time begin
	@info "wavenumber-energy_spectra"
	group = create_group(pfile, "wavenumber-energy_spectra")
	for (Re, datum) in rawdata
		kr = radialk(datum.grid)
		Ehrs = map(datum.ζhs) do ζh
			Eh = energy(ζh, datum.grid)
			radialspectrum(Eh, datum.grid)
		end
		mean, std = measure(Ehrs)
		group["Re$(Re)/x"] = Array(kr)
		group["Re$(Re)/y"] = Array(mean)
		group["Re$(Re)/Δy"] = Array(std)
	end
end

@time begin
	@info "wavenumber-normalized_energy_fluxes"
	group = create_group(pfile, "wavenumber-normalized_energy_fluxes")
	for (Re, datum) in rawdata
		kr = radialk(datum.grid)
		Fhrs = map(datum.ζhs) do ζh
			Nh = energyInjectionRate(ζh, datum.grid)
			cumradialspectrum(Nh, datum.grid)
		end
		mean, std = measure(Fhrs)
		group["Re$(Re)/x"] = Array(kr)
		group["Re$(Re)/y"] = Array(mean)
		group["Re$(Re)/Δy"] = Array(std)
	end
end

@time begin
	@info "wavenumber-enstrophy_fluxes"
	group = create_group(pfile, "wavenumber-enstrophy_fluxes")
	for (Re, datum) in rawdata
		kr = radialk(datum.grid)
		kf = datum.grid.nx / 3sqrt(Re)
		Fhrs = map(datum.ζhs) do ζh
			Nh = enstrophyInjectionRate(ζh, datum.grid)
			cumradialspectrum(Nh, datum.grid) ./ kf^2
		end
		mean, std = measure(Fhrs)
		group["Re$(Re)/x"] = Array(kr)
		group["Re$(Re)/y"] = Array(mean)
		group["Re$(Re)/Δy"] = Array(std)
	end
end

@time begin
	@info "velocity-pdf:primitive"
	group = create_group(pfile, "velocity-pdf:primitive")
	for (Re, datum) in rawdata
		grid = datum.grid

		umax = maximum(datum.ζhs) do ζh
			uh = im * ζh .* grid.invKrsq .* grid.l
			u = grid.rfftplan \ uh
			maximum(u)
		end

		bin = range(-umax, umax, 100)
		bc = (bin[1:(end-1)] .+ bin[2:end]) / 2

		pdfs = map(datum.ζhs) do ζh
			uh = im * ζh .* grid.invKrsq .* grid.l
			u = grid.rfftplan \ uh
			pdf(u, bin)
		end
		mean, std = measure(pdfs)
		group["Re$(Re)/x"] = Array(bc)
		group["Re$(Re)/y"] = Array(mean)
		group["Re$(Re)/Δy"] = Array(std)
	end
end

@time begin
	@info "velocity-pdf:filtered"
	group = create_group(pfile, "velocity-pdf:filtered")
	for (Re, datum) in rawdata
		grid = datum.grid

		umax = maximum(datum.fζhs) do ζh
			uh = im * ζh .* grid.invKrsq .* grid.l
			u = grid.rfftplan \ uh
			maximum(u)
		end

		bin = range(-umax, umax, 100)
		bc = (bin[1:(end-1)] .+ bin[2:end]) / 2

		pdfs = map(datum.fζhs) do ζh
			uh = im * ζh .* grid.invKrsq .* grid.l
			u = grid.rfftplan \ uh
			pdf(u, bin)
		end
		mean, std = measure(pdfs)
		group["Re$(Re)/x"] = Array(bc)
		group["Re$(Re)/y"] = Array(mean)
		group["Re$(Re)/Δy"] = Array(std)
	end
end

@time begin
	@info "loop_size-var_ratio:equal area"
	group = create_group(pfile, "loop_size-var_ratio:equal area")

	for (Re, datum) in rawdata
		grid = datum.grid
		dl = 2π / grid.nx
		loopsizes = (10:10:floor(Int, grid.nx/2)) .* dl
		squares = [rectsloop(grid, l, l, -π, -π) for l in loopsizes]
		rects = [rectsloop(grid, 2 * l, l / 2, -π, -π) for l in loopsizes]

		varSquares = map(squares) do hsh
			vars = map(datum.fζhs) do fζh
				Γ = getΓ(fζh, hsh, grid)
				moment(Γ, 2)
			end
			mean, std = measure(vars)
			mean ± std
		end

		varRects = map(rects) do hsh
			vars = map(datum.fζhs) do fζh
				Γ = getΓ(fζh, hsh, grid)
				moment(Γ, 2)
			end
			mean, std = measure(vars)
			mean ± std
		end

		varRatios = varRects ./ varSquares
		group["Re$(Re)/x"] = Array(loopsizes)
		group["Re$(Re)/y"] = Array(Measurements.value.(varRatios))
		group["Re$(Re)/Δy"] = Array(Measurements.uncertainty.(varRatios))
	end
end

@time begin
	@info "loop_size-var_ratio:equal perimeter"
	group = create_group(pfile, "loop_size-var_ratio:equal perimeter")

	for (Re, datum) in rawdata
		grid = datum.grid
		dl = 2π / grid.nx
		loopsizes = (10:10:floor(Int, grid.nx/2)) .* dl
		squares = [rectsloop(grid, l, l, -π, -π) for l in loopsizes]
		rects = [rectsloop(grid, 8l / 5, 2l / 5, -π, -π) for l in loopsizes]

		varSquares = map(squares) do hsh
			vars = map(datum.fζhs) do fζh
				Γ = getΓ(fζh, hsh, grid)
				moment(Γ, 2)
			end
			mean, std = measure(vars)
			mean ± std
		end

		varRects = map(rects) do hsh
			vars = map(datum.fζhs) do fζh
				Γ = getΓ(fζh, hsh, grid)
				moment(Γ, 2)
			end
			mean, std = measure(vars)
			mean ± std
		end

		varRatios = varRects ./ varSquares
		group["Re$(Re)/x"] = Array(loopsizes)
		group["Re$(Re)/y"] = Array(Measurements.value.(varRatios))
		group["Re$(Re)/Δy"] = Array(Measurements.uncertainty.(varRatios))
	end
end

@time begin
	Re = 4.0
	datum = rawdata[Re]
	@info "aspect_ratio-moments_ratio:fixed area"
	group = create_group(pfile, "aspect_ratio-moments_ratio:fixed area")

	grid = datum.grid
	dl = 2π / grid.nx
	fixarea = 8100dl^2
	loopsizes = [15; 20; 30; 40; 45; 60; 90] .* dl
	rects = [rectsloop(grid, l, fixarea / l, -π, -π) for l in loopsizes]
	orders = 2:2:10

	HDF5.attributes(group)["area"] = fixarea

	momentsRatio = map(orders) do order
		momRects = map(rects) do hsh
			moms = map(datum.fζhs) do fζh
				Γ = getΓ(fζh, hsh, grid)
				moment(Γ, order)
			end
			mean, std = measure(moms)
			mean ± std
		end
		y_mea = momRects ./ last(momRects)
		group["Re$(Re) order$(order)/x"] = Array(loopsizes .^ 2 ./ fixarea)
		group["Re$(Re) order$(order)/y"] = Array(Measurements.value.(y_mea))
		group["Re$(Re) order$(order)/Δy"] = Array(Measurements.uncertainty.(y_mea))
	end
end

@time begin
	Re = 3.5
	datum = rawdata[Re]
	@info "aspect_ratio-moments_ratio:fixed perimeter"
	group = create_group(pfile, "aspect_ratio-moments_ratio:fixed perimeter")

	grid = datum.grid
	dl = 2π / grid.nx
	fixperi = 180dl
	loopsizes = (10:10:90) .* dl
	rects = [rectsloop(grid, l, fixperi - l, -π, -π) for l in loopsizes]
	orders = 2:2:10

	HDF5.attributes(group)["perimeter"] = fixperi

	momentsRatio = map(orders) do order
		momRects = map(rects) do hsh
			moms = map(datum.fζhs) do fζh
				Γ = getΓ(fζh, hsh, grid)
				moment(Γ, order)
			end
			mean, std = measure(moms)
			mean ± std
		end
		y_mea = momRects ./ last(momRects)
		group["Re$(Re) order$(order)/x"] = Array(loopsizes .^ 2 ./ fixarea)
		group["Re$(Re) order$(order)/y"] = Array(Measurements.value.(y_mea))
		group["Re$(Re) order$(order)/Δy"] = Array(Measurements.uncertainty.(y_mea))
	end
end
