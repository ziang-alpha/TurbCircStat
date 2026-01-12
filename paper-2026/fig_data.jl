include("../core/core.jl")

# Open all datasets
DATASETS = [h5open("./output/Re$(Re)N64.h5", "r") for Re in [3.0; 3.5; 4.0]]

# Build the grid. Modify 'ngrid' for data with different resolution.
grid = TwoDGrid(CPU(); nx=64, Lx=2π)

# read each frame of the data as an iterator
getζhs(fid) = begin
    nframe = maximum(parse.(Int, keys(fid)))
    Iterators.map(1:nframe) do n
        ζh = read(fid, "$(n)")
        device_array(grid)(ζh)
    end
end



energy_spectra = map(DATASETS) do fid
    ζhs = getζhs(fid)
    kr = radialk(grid)
    # Compute the averaged 2D energy spectrum
    Ehrs = map(ζhs) do ζh
        Eh = energy(ζh, grid)
        radialspectrum(Eh, grid)
    end
    Ehr_mearsure = measure(Ehrs)
    (kr, Ehr_mearsure)
end

energy_fluxes = map(DATASETS) do fid
    ζhs = getζhs(fid)
    kr = radialk(grid)
    Fhrs = map(ζhs) do ζh
        Nh = energyInjectionRate(ζh, grid)
        cumradialspectrum(Nh, grid)
    end
    Fhr_measure = measure(Fhrs)
    (kr, Fhr_measure)
end

enstrophy_fluxes = map(DATASETS) do fid
    ζhs = getζhs(fid)
    kr = radialk(grid)
    Fhrs = map(ζhs) do ζh
        Nh = enstrophyInjectionRate(ζh, grid)
        cumradialspectrum(Nh, grid)
    end
    Fhr_measure = measure(Fhrs)
    (kr, Fhr_measure)
end

u_pdfs = map(DATASETS) do fid
    ζhs = getζhs(fid)
    bin = -50:1:50
    pdfs = map(ζhs) do ζh
        uh = im * ζh .* grid.invKrsq .* grid.l
        u = grid.rfftplan \ uh
        pdf(u, bin)
    end
    pdf_measure = measure(pdfs)
    bc = (bin[1:end-1] .+ bin[2:end]) / 2
    (bc, pdf_measure)
end

u_pdfs_filtered = map(DATASETS) do fid
    ζhs = getζhs(fid)
    bin = -50:1:50
    pdfs = map(ζhs) do ζh
        uh = im * ζh .* grid.invKrsq .* grid.l
        uh .*= (grid.Krsq .< (grid.nx / 6)^2)
        u = grid.rfftplan \ uh
        pdf(u, bin)
    end
    pdf_measure = measure(pdfs)
    bc = (bin[1:end-1] .+ bin[2:end]) / 2
    (bc, pdf_measure)
end

eqarea_variance_ratio = map(DATASETS) do fid
    ζhs = getζhs(fid)
    dl = 2π / grid.nx
    loopsizes = (10:10:floor(Int, grid.nx / 2)) .* dl
    loops_square = [rectsloop(grid, l, l, -π, -π) for l in loopsizes]
    loops_rect4 = [rectsloop(grid, 2 * l, l / 2, -π, -π) for l in loopsizes]

    var_square_measure = map(loops_square) do hsh
        vars = map(ζhs) do ζh
            Γ = getΓ(ζh, hsh, grid)
            moment(Γ, 2)
        end
        measure(vars)
    end

    var_rect4_measure = map(loops_rect4) do hsh
        vars = map(ζhs) do ζh
            Γ = getΓ(ζh, hsh, grid)
            moment(Γ, 2)
        end
        measure(vars)
    end

    var_ratio_measure = var_rect4_measure ./ var_square_measure
    (loopsizes, var_ratio_measure)
end

eqperi_variance_ratio = map(DATASETS) do fid
    ζhs = getζhs(fid)
    dl = 2π / grid.nx
    loopsizes = (10:10:floor(Int, grid.nx / 2)) .* dl
    loops_square = [rectsloop(grid, l, l, -π, -π) for l in loopsizes]
    loops_rect4 = [rectsloop(grid, 8l / 5, 2l / 5, -π, -π) for l in loopsizes]

    var_square_measure = map(loops_square) do hsh
        vars = map(ζhs) do ζh
            Γ = getΓ(ζh, hsh, grid)
            moment(Γ, 2)
        end
        measure(vars)
    end

    var_rect4_measure = map(loops_rect4) do hsh
        vars = map(ζhs) do ζh
            Γ = getΓ(ζh, hsh, grid)
            moment(Γ, 2)
        end
        measure(vars)
    end

    var_ratio_measure = var_rect4_measure ./ var_square_measure
    (loopsizes, var_ratio_measure)
end

order = 2:2:10
eqarea_moment_ratio = map(order) do n
    ζhs = getζhs(DATASETS[1])
    dl = 2π / grid.nx
    loopsizes = [15; 20; 30; 40; 45; 60; 90] .* dl
    loop_square = rectsloop(grid, 90dl, 90dl, -π, -π)
    loops_rect4 = [rectsloop(grid, l, 8100dl^2 / l, -π, -π) for l in loopsizes]

    var_square_measure = begin
        vars = map(ζhs) do ζh
            Γ = getΓ(ζh, loop_square, grid)
            moment(Γ, 2)
        end
        measure(vars)
    end

    var_rect4_measure = map(loops_rect4) do hsh
        vars = map(ζhs) do ζh
            Γ = getΓ(ζh, hsh, grid)
            moment(Γ, 2)
        end
        measure(vars)
    end

    var_ratio_measure = var_rect4_measure ./ var_square_measure
    (loopsizes .^ 2 / 8100dl^2, var_ratio_measure)
end

eqperi_moment_ratio = map(order) do n
    ζhs = getζhs(DATASETS[1])
    dl = 2π / grid.nx
    loopsizes = (10:10:90) .* dl
    loop_square = rectsloop(grid, 90dl, 90dl, -π, -π)
    loops_rect4 = [rectsloop(grid, l, 180dl - l, -π, -π) for l in loopsizes]

    var_square_measure = begin
        vars = map(ζhs) do ζh
            Γ = getΓ(ζh, loop_square, grid)
            moment(Γ, 2)
        end
        measure(vars)
    end

    var_rect4_measure = map(loops_rect4) do hsh
        vars = map(ζhs) do ζh
            Γ = getΓ(ζh, hsh, grid)
            moment(Γ, 2)
        end
        measure(vars)
    end

    var_ratio_measure = var_rect4_measure ./ var_square_measure
    (loopsizes ./ (180dl .- loopsizes), var_ratio_measure)
end

@save "./output/figure_data.jld2" energy_spectra energy_fluxes enstrophy_fluxes u_pdfs u_pdfs_filtered eqarea_variance_ratio eqperi_variance_ratio eqarea_moment_ratio eqperi_moment_ratio

close.(DATASETS)
