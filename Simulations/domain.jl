# This file specifies the computational domain
# we are looking at a 3D problem with 256 grid spaces in each direction
# N_x = N_y = N_z = 256
# Number of grid points is 257 for each direction, beware of the trivia
# Nx = Ny = Nz = 257


include("Assembling_3D_matrices.jl")

N_x = N_y = N_z = 32   # Assuming the same number of grids in each direction
Nx = N_x + 1
Ny = N_y + 1
Nz = N_z + 1

# RS friction region test purpose at this moment
# fN2 and fN3 represents number of points in the fault region
fN2 = 16
fN3 = 9

# fNy and fNz represents the indices of the fault region
fNy = (1, 16)                    # y direction 
fNz = (2, 10)                   # z direction

SBPp = 2                # SBPp order
M, RHS, H_tilde, HI_tilde, analy_sol, source = Assembling_3D_matrices(N_x, N_y, N_z;p=SBPp);


# setting up dψV, ψδ in odefun for ODEProblem()

# size of ψ, dψ: (rate-and-state portion of the fault): (fN2 + 1) * (fN3 + 1)
# fN2 + 1 ≈ 41      fN3 + 1 ≈ 101

# size of V: 2 * (N3 + 1) * (N2 + 1)
# size of δ: 2 * (N3 + 1) * (N2 + 1) 
# size of dψV = size(dψ) + size(V) = (fN2 + 1) * (fN3 + 1) + 2 * (N3 + 1) * (N2 + 1)

# dψV = zeros((fN2 + 1) * (fN3 + 1) + 2 * (Ny) * (Nz))
# ψδ = zeros((fN2 + 1) * (fN3 + 1) + 2 * (Ny) * (Nz))

function create_ψVδ()
    dψV = zeros((fN2 + 1) * (fN3 + 1) + 2 * (Ny) * (Nz))
    ψδ = zeros((fN2 + 1) * (fN3 + 1) + 2 * (Ny) * (Nz))
    return dψV, ψδ
end

function create_view(dψV, ψδ)
    index_1 = (fN2 + 1)*(fN3 + 1)

    dψ = @view dψV[1:index_1]
    V = @view dψV[index_1 + 1:end]
    ψ = @view ψδ[1:index_1]
    δ = @view ψδ[index_1 + 1:end]

    return dψ, V, ψ, δ
end

function get_ψ_indices(Nx, Ny, Nz, fNy, fNz)
    # fNy is a tuple of the start index and end index of fault in y(x2) direction 
    # fNz is a tuple of the start index and end index of fault in z(x3) direction 
    x_idx = spzeros(Nx)
    y_idx = spzeros(Ny)
    z_idx = spzeros(Nz)
    x_idx[1] = 1 # the fault is on the first face
    y_idx[fNy[1]:fNy[2]] .= 1
    z_idx[fNz[1]:fNz[2]] .= 1
    return kron(z_idx, y_idx, x_idx) # x->y->z is the index order, hence the kron order is z_idx<-y_idx<-x_idx
end

test_ψ_indices = get_ψ_indices(Nx, Ny, Nz, fNy, fNz);
reshape(test_ψ_indices.nzind,fN2,fN3)'

nothing # avoid printing out results 