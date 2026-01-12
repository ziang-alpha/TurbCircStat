PATH = dirname(@__FILE__)
ENV["JULIA_STDOUT_LINE_BUFFERED"] = "1"
# Parameters for the simulation

ngrid       = length(ARGS) ≥ 1 ? parse(Int, ARGS[1])     : 64         # Grid number on each direction
Re          = length(ARGS) ≥ 2 ? parse(Float64, ARGS[2]) : 4.0    # Reynolds number, Re = kf^(-4/3) * ε^(1/3) * ν^(-1) 
ndata       = length(ARGS) ≥ 3 ? parse(Int, ARGS[3])     : 400        # total frame of data to be saved
nstep       = length(ARGS) ≥ 4 ? parse(Int, ARGS[4])     : 100        # number of steps between each frame saved
CFL_forcing = length(ARGS) ≥ 5 ? parse(Float64, ARGS[5]) : 1e-2 # CFL defined by the forcing parameters. dt = ε^(-1/3) * kf^(-2/3) * CFL_forcing

ε  = 1.0     # Energy injection rate
L  = 2π  # domain size
kf = ngrid / (3 * sqrt(Re))       # Forcing wavenumber 
ν  = ε^(1 / 3) * kf^(-4 / 3) / Re    # Viscosity

dev = CPU() # Device of computation: GPU() or CPU()

output_path = PATH * "/.output/Re$(Re)_N$(ngrid).h5"   # path of the output file
diag_path   = PATH * "/output/Re$(Re)_N$(ngrid)_diag.h5"

isdir(PATH * "/.output/") || mkdir(PATH * "/.output/")
isdir(PATH * "/output/")  || mkdir(PATH * "/output/")