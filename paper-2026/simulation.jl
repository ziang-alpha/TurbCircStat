using GeophysicalFlows, Random, HDF5, ProgressMeter

include("parameters.jl")
include("../core/diagnostics.jl")

# Get the forcing spectrum
forcing_spectrum = begin
    forcing_wavenumber = kf * 2π / L  # the forcing wavenumber, `k_f`, for a spectrum that is a ring in wavenumber space
    forcing_bandwidth = 1.5 * 2π / L  # the width of the forcing spectrum, `δ_f`

    grid = TwoDGrid(CPU(); nx=ngrid, Lx=L)

    K = @. sqrt(grid.Krsq)             # a 2D array with the total wavenumber

    forcing_spectrum = @. exp(-(K - forcing_wavenumber)^2 / (2 * forcing_bandwidth^2))

    forcing_spectrum[grid.Krsq.==0] .= 0 # ensure forcing has zero domain-average

    ε0 = FourierFlows.parsevalsum(forcing_spectrum .* grid.invKrsq / 2, grid) / (grid.Lx * grid.Ly)
    forcing_spectrum .*= ε / ε0        # normalize forcing to inject energy at rate ε	
    device_array(dev)(forcing_spectrum)
end

# Functions describing the random forcing
function calcF!(Fh, sol, t, clock, vars, params, grid)
    randn!(Fh)
    @. Fh *= sqrt(forcing_spectrum) / sqrt(clock.dt)
    return nothing
end

# Instantiate the TwoDNavierStokes problem
prob = TwoDNavierStokes.Problem(dev;
    nx=ngrid,
    Lx=L,
    ν,
    dt=ε^(-1 / 3) * kf^(-2 / 3) * CFL_forcing,
    stepper="ETDRK4",
    calcF=calcF!,
    stochastic=true,
)

# Initialize out put file & get the initial value.
if isfile(output_path)
    sol = h5open(output_path, "r") do f
        frames = [parse(Int, k) for k in keys(f)]
        read(f, "$(maximum(frames))")
    end
    ζ₀ = prob.grid.rfftplan \ device_array(dev)(sol)
    set_ζ!(prob, ζ₀)
else
    set_ζ!(prob, device_array(dev)(zeros(ngrid, ngrid)))
end

fid = h5open(output_path, "w")
fid_diag = h5open(diag_path, "w")
diag = create_dataset(fid_diag, "energy", Float64, (ndata))

p = Progress(ndata, dt=0.1, desc="Progress:", barlen=0, color=:white)
for nframe in 1:ndata
    stepforward!(prob, nstep)
    write_dataset(fid, "$(nframe)", Array(prob.sol))
    diag[nframe] = FourierFlows.parsevalsum(energy(prob.sol, prob.grid), prob.grid)
    next!(p)
end
finish!(p)

close(fid_diag)
close(fid)
