# Solving a Poisson Equation in 3D
#
#
#   Governing equations: -(Dxx + Dyy + Dzz) u = f
#   Boundary conditions:
#   1. Free surface on top and bottom face (denoted T and D)
#   2. Dirichlet boundary conditions on left (fault) and right face (denoted L and R)
#   3. Free boundary conditions on front and back (denoted F and B)
#
#   Manufactured solutions



include("diagonal_sbp.jl")
include("3D_face.jl")
include("analy_sol.jl")

using LinearAlgebra
using IterativeSolvers
using BenchmarkTools
using Plots

p = 2
i = j = k = 2
h_list_x = [1/2^2, 1/2^3, 1/2^4, 1/2^5, 1/2^6, 1/2^7, 1/2^8,1/2^9,1/2^10]
h_list_y = [1/2^2, 1/2^3, 1/2^4, 1/2^5, 1/2^6, 1/2^7, 1/2^8,1/2^9,1/2^10]
h_list_z = [1/2^2, 1/2^3, 1/2^4, 1/2^5, 1/2^6, 1/2^7, 1/2^8,1/2^9,1/2^10]


hx = h_list_x[i];
hy = h_list_y[j];
hz = h_list_z[k];

x = range(0,step=hx,1);
y = range(0,step=hy,1);
z = range(0,step=hz,1);

x_list = 1 ./h_list_x;
y_list = 1 ./h_list_y;
z_list = 1 ./h_list_z;


# Matrix Size
N_x = Integer(x_list[i]);
N_y = Integer(y_list[j]);
N_z = Integer(z_list[k]);

(D1x, HIx, H1x, r1x) = diagonal_sbp_D1(p,N_x,xc=(0,1));
(D2x, S0x, SNx, HI2x, H2x, r2x) = diagonal_sbp_D2(p,N_x,xc=(0,1));


(D1y, HIy, H1y, r1y) = diagonal_sbp_D1(p,N_y,xc=(0,1));
(D2y, S0y, SNy, HI2y, H2y, r2y) = diagonal_sbp_D2(p,N_y,xc=(0,1));


(D1z, HIz, H1z, r1z) = diagonal_sbp_D1(p,N_y,xc=(0,1));
(D2z, S0z, SNz, HI2z, H2z, r2z) = diagonal_sbp_D2(p,N_y,xc=(0,1));

BSx = sparse(SNx - S0x);
BSy = sparse(SNy - S0y);
BSz = sparse(SNz - S0z);

# Forming 2d Operators
# e_1x = sparse(e(1,N_x+1));
# e_Nx = sparse(e(N_x+1,N_x+1));
# e_1y = sparse(e(1,N_y+1));
# e_Ny = sparse(e(N_y+1,N_y+1));
# e_1z = sparse(e(1,N_z+1));
# e_Nz = sparse(e(N_z+1,N_z+1));


# I_Nx = sparse(eyes(N_x+1));
# I_Ny = sparse(eyes(N_y+1));
# I_Nz = sparse(eyes(N_z+1));

# e_E = kron(e_Nx,I_Ny);
# e_W = kron(e_1x,I_Ny);
# e_S = kron(I_Nx,e_1y);
# e_N = kron(I_Nx,e_Ny);

# Forming 3D operators

D2_x = kron(D2x,I_Ny,I_Ny)
D2_y = kron(I_Ny,D2y,I_Ny)
D2_z = kron(I_Ny,I_Ny,D2z)

e_1x = sparse(e(1,N_x+1))
e_1y = sparse(e(1,N_y+1))
e_1z = sparse(e(1,N_z+1))

e_Nx = sparse(e(N_x+1,N_x+1))
e_Ny = sparse(e(N_y+1,N_y+1))
e_Nz = sparse(e(N_z+1,N_z+1))


# Penalty Parameters
tau_x = 13/hx
tau_y = 13/hy
tau_z = 13/hz

beta = -1

# SBP-SAT Terms
SAT_L = tau_x * 



# Analytical solutions
(analy_sol,x_ex,y_ex,z_ex) = form_analy_sol(;N=2^3)
plot(x_ex,y_ex,analy_sol[:,:,end],st=:surface)

# source Terms
source_terms = 3π^2 * analy_sol;
