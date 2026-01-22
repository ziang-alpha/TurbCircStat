fid = h5open("diagnostics.jl", "w")
energy_spectra = map(raw_data) do (ζhs, grid, fh)
    kr = radialk(grid)
    # Compute the averaged 2D energy spectrum
    Ehrs = map(ζhs) do ζh
        Eh = energy(ζh, grid)
        radialspectrum(Eh, grid)
    end
    Ehr_mearsure = measure(Ehrs)

end

energy_fluxes = map(raw_data) do (ζhs, grid, fh)
    kr = radialk(grid)
    Fhrs = map(ζhs) do ζh
        Nh = energyInjectionRate(ζh, grid)
        cumradialspectrum(Nh, grid)
    end
    Fhr_measure = measure(Fhrs)
    (kr, Fhr_measure)
end

enstrophy_fluxes = map(raw_data) do (ζhs, grid, fh)
    kr = radialk(grid)
    Fhrs = map(ζhs) do ζh
        Nh = enstrophyInjectionRate(ζh, grid)
        cumradialspectrum(Nh, grid)
    end
    Fhr_measure = measure(Fhrs)
    (kr, Fhr_measure)
end

u_pdfs = map(raw_data) do (ζhs, grid, fh)
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

u_pdfs_filtered = map(raw_data) do (ζhs, grid, fh)
    bin = -50:1:50
    pdfs = map(ζhs) do ζh
        uh = im * ζh .* grid.invKrsq .* grid.l
        uh .*= fh
        u = grid.rfftplan \ uh
        pdf(u, bin)
    end
    pdf_measure = measure(pdfs)
    bc = (bin[1:end-1] .+ bin[2:end]) / 2
    (bc, pdf_measure)
end

@save (@__DIR__) * "/diagnostics.jld2" energy_spectra energy_fluxes enstrophy_fluxes u_pdfs u_pdfs_filtered
