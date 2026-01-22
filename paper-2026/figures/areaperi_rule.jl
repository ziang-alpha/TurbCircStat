eqarea_variance_ratio = map(raw_data) do (ζhs, grid, fh)
    dl = 2π / grid.nx
    loopsizes = (10:10:floor(Int, grid.nx / 2)) .* dl
    loops_square = [rectsloop(grid, l, l, -π, -π) for l in loopsizes]
    loops_rect4 = [rectsloop(grid, 2 * l, l / 2, -π, -π) for l in loopsizes]

    var_square_measure = map(loops_square) do hsh
        vars = map(ζhs) do ζh
            Γ = getΓ(ζh .* fh, hsh, grid)
            moment(Γ, 2)
        end
        measure(vars)
    end

    var_rect4_measure = map(loops_rect4) do hsh
        vars = map(ζhs) do ζh
            Γ = getΓ(ζh .* fh, hsh, grid)
            moment(Γ, 2)
        end
        measure(vars)
    end

    var_ratio_measure = var_rect4_measure ./ var_square_measure
    (loopsizes, var_ratio_measure)
end

eqperi_variance_ratio = map(raw_data) do (ζhs, grid, fh)
    dl = 2π / grid.nx
    loopsizes = (10:10:floor(Int, grid.nx / 2)) .* dl
    loops_square = [rectsloop(grid, l, l, -π, -π) for l in loopsizes]
    loops_rect4 = [rectsloop(grid, 8l / 5, 2l / 5, -π, -π) for l in loopsizes]

    var_square_measure = map(loops_square) do hsh
        vars = map(ζhs) do ζh
            Γ = getΓ(ζh .* fh, hsh, grid)
            moment(Γ, 2)
        end
        measure(vars)
    end

    var_rect4_measure = map(loops_rect4) do hsh
        vars = map(ζhs) do ζh
            Γ = getΓ(ζh .* fh, hsh, grid)
            moment(Γ, 2)
        end
        measure(vars)
    end

    var_ratio_measure = var_rect4_measure ./ var_square_measure
    (loopsizes, var_ratio_measure)
end

order = 2:2:10

eqarea_moment_ratio = map(order) do n
    ζhs, grid, fh = raw_data[1]
    dl = 2π / grid.nx
    loopsizes = [15; 20; 30; 40; 45; 60; 90] .* dl
    loop_square = rectsloop(grid, 90dl, 90dl, -π, -π)
    loops_rect4 = [rectsloop(grid, l, 8100dl^2 / l, -π, -π) for l in loopsizes]
    mom_square_measure = begin
        vars = map(ζhs) do ζh
            Γ = getΓ(ζh .* fh, loop_square, grid)
            moment(Γ, n)
        end
        measure(vars)
    end
    mom_rect4_measure = map(loops_rect4) do hsh
        vars = map(ζhs) do ζh
            Γ = getΓ(ζh .* fh, hsh, grid)
            moment(Γ, n)
        end
        measure(vars)
    end
    mom_ratio_measure = mom_rect4_measure ./ mom_square_measure
    (loopsizes .^ 2 / 8100dl^2, mom_ratio_measure)
end

eqperi_moment_ratio = map(order) do n
    ζhs, grid, fh = raw_data[1]
    dl = 2π / grid.nx
    loopsizes = (10:10:90) .* dl
    loop_square = rectsloop(grid, 90dl, 90dl, -π, -π)
    loops_rect4 = [rectsloop(grid, l, 180dl - l, -π, -π) for l in loopsizes]

    mom_square_measure = begin
        vars = map(ζhs) do ζh
            Γ = getΓ(ζh .* fh, loop_square, grid)
            moment(Γ, 2)
        end
        measure(vars)
    end

    mom_rect4_measure = map(loops_rect4) do hsh
        vars = map(ζhs) do ζh
            Γ = getΓ(ζh .* fh, hsh, grid)
            moment(Γ, 2)
        end
        measure(vars)
    end

    mom_ratio_measure = mom_rect4_measure ./ mom_square_measure
    (loopsizes ./ (180dl .- loopsizes), mom_ratio_measure)
end

@save (@__DIR__) * "/areaperi_rule.jld2" eqarea_variance_ratio eqperi_variance_ratio eqarea_moment_ratio eqperi_moment_ratio order

