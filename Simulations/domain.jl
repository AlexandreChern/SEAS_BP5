# This file specifies the computational domain
# we are looking at a 3D problem with 256 grid spaces in each direction
# N_x = N_y = N_z = 256
# Number of grid points is 257 for each direction, beware of the trivia
# Nx = Ny = Nz = 257


include("Assembling_3D_matrices.jl")

N_x = N_y = N_z = 32   # Assuming the same number of grids in each direction
SBPp = 2                # SBPp order
M, RHS, H_tilde, HI_tilde, analy_sol, source = Assembling_3D_matrices(N_x, N_y, N_z;p=SBPp);

nothing # avoid printing out results 