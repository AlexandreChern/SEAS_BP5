# 3D array indexing in Julia
#
#    3
#    |
#    |
#    |___________  2
#   /
#  /
# 1
#
# Example: julia> test_3D
# 3×3×3 Array{Float64, 3}:
# [:, :, 1] =
# 1.99217   -0.419911   1.27037
# -0.616689  -0.615541  -0.452607
# -0.144164  -1.93718    0.418154

# [:, :, 2] =
# 3.03034   0.37151    0.964773
# 1.40377  -0.338357  -0.2258
# 1.21607   0.57358   -0.364077

# [:, :, 3] =
# -0.417334   -0.675461   -1.81201
# -0.0343112   0.0751002   0.562684
# -0.326663    0.949942    0.587156
#
#

#   julia> test_3D[:]
#   27-element Vector{Float64}:
#   1.9921749921668979
#  -0.6166891690372321
#  -0.1441644800574074
#  -0.4199112178717434
#  -0.6155411265727337
#  -1.9371793724089905
#   1.270365705663953
#  -0.45260724271602304
#   0.41815434826401493
#   3.0303354534340596
#   1.4037652552774593
#   1.2160687591726536
#   0.3715100914804111
#  -0.3383569289529389
#   0.5735799032430334
#   0.9647725600192497
#  -0.22579997325200207
#  -0.3640769675163103
#  -0.4173335304749411
#  -0.03431124808637245
#  -0.326663027771571
#  -0.6754607583657409
#   0.07510016190239746
#   0.9499416235735255
#  -1.8120098781175735
#   0.5626840107886724
#   0.5871558963442842

using LinearAlgebra
using SparseArrays

function e(i,n)
    A = Matrix{Float64}(I,n,n)
    return A[:,i]
end

function eyes(n)
    return Matrix{Float64}(I,n,n)
end


function Diag(A)
    # Self defined function that is similar to Matlab Diag
    return Diagonal(A[:])
end


function get_bottom_face(Nx,Ny,Nz)
    mat = kron(e(1,Nz)',sparse(eyes(Nx)),sparse(eyes(Ny)))
    return mat
end

function get_top_face(Nx,Ny,Nz)
    mat = kron(e(Nz,Nz)',sparse(eyes(Nx)),sparse(eyes(Ny)))
    return mat
end

function get_left_face(Nx,Ny,Nz)
    mat = kron(sparse(eyes(Nz)),e(1,Ny)',sparse(eyes(Nx)))
    return mat
end


function get_right_face(Nx,Ny,Nz)
    mat = kron(sparse(eyes(Nz)),e(Ny,Ny)',sparse(eyes(Nx)))
    return mat
end

function get_front_face(Nx,Ny,Nz)
    mat = kron(sparse(eyes(Ny)),sparse(eyes(Nz)),sparse(e(1,Nx)'))
    return mat
end


function get_back_face(Nx,Ny,Nz)
    mat = kron(sparse(eyes(Ny)),sparse(eyes(Nz)),sparse(e(Nx,Nx)'))
    return mat
end


test_matrix = randn(2,3,4)