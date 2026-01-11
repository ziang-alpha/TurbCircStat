# Parameters for the simulation

ngrid = 64  # Grid number on each direction
Re = 3.0    # Reynolds number, Re = kf^(-4/3) * ε^(1/3) * ν^(-1) 
ε = 1.0     # Energy injection rate

CFL_forcing = 1e-2 # CFL defined by the forcing parameters. dt = ε^(-1/3) * kf^(-2/3) * CFL_forcing

ndata = 400 # total frame of data to be saved
nstep = 100 # number of steps between each frame saved

output_path = "F:/Re$(Re).h5"   # path of the output file

L = 2π  # domain size

kf = ngrid / (3*sqrt(Re))       # Forcing wavenumber 
ν = ε^(1/3) * kf^(-4/3) / Re    # Viscosity

dev = GPU() # Device of computation: GPU() or CPU()






