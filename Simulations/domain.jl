# This file specifies the computational domain
# we are looking at a 3D problem with 256 grid spaces in each direction
# N_x = N_y = N_z = 256
# Number of grid points is 257 for each direction, beware of the trivia
# Nx = Ny = Nz = 257


include("Assembling_3D_matrices.jl")
include("coefficients.jl")

N_x = N_y = N_z =64   # Assuming the same number of grids in each direction
Nx = N_x + 1
Ny = N_y + 1
Nz = N_z + 1

# RS friction region test purpose at this moment
# fN2 and fN3 represents number of points in the fault region

# fN2 = 51 # l_f / (2 * Δz) + 1
# fN3 = 41 # W_f / (2 * Δz) + 1

fN2 = BP5_coeff.lf / (2 * BP5_coeff.Δz / 1000) + 1
fN2 = round(Int, fN2, RoundUp)
fN3 = BP5_coeff.Wf / (BP5_coeff.Δz / 1000) + 1
fN3 = round(Int, fN3, RoundUp)

# fNy and fNz represents the indices of the fault region
fNy = (1, fN2)                    # y direction 
fNz = (1, fN3)                   # z direction

# VW region of the RS

fN2_VW = BP5_coeff.l /(2 * BP5_coeff.Δz / 1000) + 1
fN2_VW = round(Int, fN2_VW, RoundUp)
fN3_VW = BP5_coeff.H / (BP5_coeff.Δz / 1000) + 1
fN3_VW = round(Int, fN3_VW, RoundUp)

# fNy_VW and fNz_VW represents the indices of the fault region
fNy_VW = (1, fN2_VW)
fNz_VW_start = round(Int, (BP5_coeff.hs + BP5_coeff.ht) / (BP5_coeff.Δz / 1000), RoundUp) + 1
fNz_VW = (fNz_VW_start, fNz_VW_start + fN3_VW - 1)


# VW-VS transition region
fN2_VW_VS = (BP5_coeff.l /2 + BP5_coeff.ht) / (BP5_coeff.Δz / 1000) 
fN2_VW_VS = round(Int, fN2_VW_VS, RoundUp)

fN3_VW_VS = (BP5_coeff.H + BP5_coeff.ht * 2) / (BP5_coeff.Δz / 1000) + 1
fN3_VW_VS = round(Int, fN3_VW_VS, RoundUp)

fNy_VW_VS = (1, fN2_VW_VS)
fNz_VW_VS_start = round(Int, BP5_coeff.hs / (BP5_coeff.Δz / 1000), RoundUp) + 1
fNz_VW_VS = (fNz_VW_VS_start, fNz_VW_VS_start + fN3_VW_VS - 1)

# fNy_VW_VS and fNz_VW represents the indices of the fault region

# Assembling matrices for 3D SBP-SAT
SBPp = 2                # SBPp order
M, RHS, H_tilde, HI_tilde, analy_sol, source = Assembling_3D_matrices(N_x, N_y, N_z;p=SBPp);


# setting up dψV, ψδ in odefun for ODEProblem()

# size of ψ, dψ: (rate-and-state portion of the fault): (fN2 + 1) * (fN3 + 1)
# fN2 + 1 ≈ 41      fN3 + 1 ≈ 101

# size of V: 2 * (N3 + 1) * (N2 + 1)
# size of δ: 2 * (N3 + 1) * (N2 + 1) 
# size of dψV = size(dψ) + size(V) = (fN2 + 1) * (fN3 + 1) + 2 * (N3 + 1) * (N2 + 1)
# dψV = zeros((fN2 + 1) * (fN3 + 1) + 2 * (Ny) * (Nz))

# size of ψδ = size(ψ) + size(δ) =  (fN2 + 1) * (fN3 + 1) + 2 * (N3 + 1) * (N2 + 1)
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

# getting rate-and-state 
function get_RS_indices(Nx, Ny, Nz, fNy, fNz)
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

function get_uniform_indices(Nx, Ny, Nz, fNy_VW, fNz_VW)
    x_idx = spzeros(Nx)
    y_idx = spzeros(Ny)
    z_idx = spzeros(Nz)
    x_idx[1] = 1 # the fault is on the first face
    y_idx[fNy_VW[1]:fNy_VW[2]] .= 1
    z_idx[fNz_VW[1]:fNz_VW[2]] .= 1
    return kron(z_idx, y_idx, x_idx)
end

function get_transition_indices(Nx, Ny, Nz, fNy_VW, fNz_VW, fNy_VW_VS, fNz_VW_VS)
    x_idx = spzeros(Nx)
    y_idx = spzeros(Ny)
    z_idx = spzeros(Nz)
    x_idx[1] = 1 # the fault is on the first face
    y_idx[fNy_VW_VS[1]:fNy_VW_VS[2]] .= 1
    z_idx[fNz_VW_VS[1]:fNz_VW_VS[2]] .= 1
    return kron(z_idx, y_idx, x_idx) - get_uniform_indices(Nx, Ny, Nz, fNy_VW, fNz_VW)
end

let 
    test_RS_indices = get_RS_indices(Nx, Ny, Nz, fNy, fNz);
    reshape(test_RS_indices.nzind,fN2,fN3)'

    test_uniform_indices = get_uniform_indices(Nx, Ny, Nz, fNy_VW, fNz_VW);
    reshape(test_uniform_indices.nzind, fN2_VW, fN3_VW)'

    test_transition_indices = get_transition_indices(Nx, Ny, Nz, fNy_VW, fNz_VW, fNy_VW_VS, fNz_VW_VS)
    test_uniform_indices.nzind

    reshape((test_transition_indices + test_uniform_indices).nzind, fN2_VW_VS, fN3_VW_VS)'
end

RS_filter = get_RS_indices(Nx, Ny, Nz, fNy, fNz)
VW_filter = get_uniform_indices(Nx, Ny, Nz, fNy_VW, fNz_VW)
VW_VS_transition_filter = get_transition_indices(Nx, Ny, Nz, fNy_VW, fNz_VW, fNy_VW_VS, fNz_VW_VS)
VS_filter = RS_filter - VW_filter - VW_VS_transition_filter

# assert these three regions have the total points of the entire rate and state
# friction region
@assert length(VW_filter.nzind) + length(VW_VS_transition_filter.nzind) + length(VS_filter.nzind) == fN2 * fN3

nothing # avoid printing out results 