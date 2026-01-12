using FourierFlows

# Compute the radial spectrum
radial(fh, grid) = FourierFlows.radialspectrum(fh, grid; refinement=1)
# Compute the cumulate radial spectrum (e.g. fluxes)
radialcum(fh, grid) = begin
    (kr, fhr) = radialspectrum(fh, grid)
    fhrc = cumsum(fhr) * kr.step
    (kr, fhr)
end

# Diagnostic Quantities for 2D Navier Stokes equations
enstrophy(ζh, grid) = abs2.(ζh) / (2 * grid.nx * grid.ny)

energy(ζh, grid) = abs2.(ζh) .* grid.invKrsq / (2 * grid.nx * grid.ny)

enstrophyInjectionRate(ζh, grid) = begin
    ζx = grid.rfftplan \ (@. im * grid.kr * ζh)
    ζy = grid.rfftplan \ (@. im * grid.l * ζh)
    u = grid.rfftplan \ (@. im * grid.l * grid.invKrsq * ζh)
    v = grid.rfftplan \ (@. -im * grid.kr * grid.invKrsq * ζh)
    Nh = grid.rfftplan * (@. u * ζx + v * ζy)
    real.(Nh .* conj(ζh))
end

energyInjectionRate(ζh, grid) = begin
    ζx = grid.rfftplan \ (@. im * grid.kr * ζh)
    ζy = grid.rfftplan \ (@. im * grid.l * ζh)
    u = grid.rfftplan \ (@. im * grid.l * grid.invKrsq * ζh)
    v = grid.rfftplan \ (@. -im * grid.kr * grid.invKrsq * ζh)
    Nh = grid.rfftplan * (@. u * ζx + v * ζy)
    real.(Nh .* conj(ζh) .* grid.invKrsq)
end


