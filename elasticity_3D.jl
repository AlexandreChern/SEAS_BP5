include("diagonal_sbp.jl")
include("3D_face.jl")
include("analy_sol.jl")
include("components.jl")
include("coefficients.jl")

using LinearAlgebra
using IterativeSolvers
using BenchmarkTools
using Plots

p = 2

level = 3

i = j = k = level
h_list_x = [1/2^1, 1/2^2, 1/2^3, 1/2^4, 1/2^5, 1/2^6, 1/2^7, 1/2^8,1/2^9,1/2^10]
h_list_y = [1/2^1, 1/2^2, 1/2^3, 1/2^4, 1/2^5, 1/2^6, 1/2^7, 1/2^8,1/2^9,1/2^10]
h_list_z = [1/2^1, 1/2^2, 1/2^3, 1/2^4, 1/2^5, 1/2^6, 1/2^7, 1/2^8,1/2^9,1/2^10]


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

Nx = N_x + 1
Ny = N_y + 1
Nz = N_z + 1

(D1x, HIx, H1x, r1x) = diagonal_sbp_D1(p,N_x,xc=(0,1));
(D2x, S0x, SNx, HI2x, H2x, r2x) = diagonal_sbp_D2(p,N_x,xc=(0,1));


(D1y, HIy, H1y, r1y) = diagonal_sbp_D1(p,N_y,xc=(0,1));
(D2y, S0y, SNy, HI2y, H2y, r2y) = diagonal_sbp_D2(p,N_y,xc=(0,1));


(D1z, HIz, H1z, r1z) = diagonal_sbp_D1(p,N_y,xc=(0,1));
(D2z, S0z, SNz, HI2z, H2z, r2z) = diagonal_sbp_D2(p,N_y,xc=(0,1));

BSx = sparse(SNx - S0x);
BSy = sparse(SNy - S0y);
BSz = sparse(SNz - S0z);

H_tilde = kron(H1x,H1y,H1z)

I_Nx = eyes(N_x+1)
I_Ny = eyes(N_y+1)
I_Nz = eyes(N_z+1)

D2_x = kron(I_Nz,I_Ny,D2x)
D2_y = kron(I_Nz,D2y,I_Nx)
D2_z = kron(D2z,I_Nx,I_Ny)

e_1x = e(1,N_x + 1)
e_1y = e(1,N_y + 1)
e_1z = e(1,N_z + 1)

e_Nx = e(N_x + 1,N_x + 1)
e_Ny = e(N_y + 1,N_y + 1)
e_Nz = e(N_z + 1,N_z + 1)

HI_Front = kron(I_Ny,I_Nz,HIx)
HI_End = kron(I_Ny,I_Nz,HIx)

HI_Left = kron(I_Nz,HIy,I_Nx)
HI_Right = kron(I_Nz,HIy,I_Nx)

HI_Bottom = kron(HIz,I_Nx,I_Ny)
HI_Top = kron(HIz,I_Nx,I_Ny)

BS_Front = kron(I_Ny,I_Nz,BSx)
BS_End = kron(I_Ny,I_Nz,BSx)

BS_Left = kron(I_Nz,BSy,I_Nx)
BS_Right = kron(I_Nz,BSy,I_Nx)

# BS_Bottom = kron(I_Nx,I_Ny,BSz)
# BS_Top = kron(I_Nx,I_Ny,BSz)
BS_Bottom = kron(BSz,I_Nx,I_Ny)
BS_Top = kron(BSz,I_Nx,I_Ny)

Front_operator = get_front_face(N_x+1,N_y+1,N_z+1)
End_operator = get_end_face(N_x+1,N_y+1,N_z+1)
Left_operator = get_left_face(N_x+1,N_y+1,N_z+1)
Right_operator = get_right_face(N_x+1,N_y+1,N_z+1)
Top_operator = get_top_face(N_x+1,N_y+1,N_z+1)
Bottom_operator = get_bottom_face(N_x+1,N_y+1,N_z+1)

# Penalty Parameters
tau_x = 13/hx
tau_y = 13/hy
tau_z = 13/hz


beta = -1



# A = -(D2_x + D2_y + D2_z) + SAT_Front + SAT_End + SAT_Left + SAT_Right + SAT_Bottom + SAT_Top

# b = source_terms[:] + SAT_Front_r * G_Front + SAT_End_r * G_End + SAT_Left_r * G_Left + SAT_Right_r * G_Right + SAT_Bottom_r * G_Bottom + SAT_Top_r * G_Top


# numerical_sol = A\b

# numerical_sol_3D = reshape(numerical_sol,N_x+1,N_y+1,N_z+1)

# analy_sol_3D u1 = u2 = u3 = sin(πx + πy + πz)

u1_filter = get_u1(Nx,Ny,Nz)
u2_filter = get_u2(Nx,Ny,Nz)
u3_filter = get_u3(Nx,Ny,Nz)


analy_sol = zeros(3*Nx*Ny*Nz)

for i in eachindex(numerical_sol)
    analy_sol[i] = rem(i,3)
end

# setting values for u1 u2 u3
analy_sol[1:3:end] = analy_sol_3D[:] 
analy_sol[2:3:end] = analy_sol_3D[:]
analy_sol[3:3:end] = analy_sol_3D[:] 

u1 = u1_filter * analy_sol
u2 = u2_filter * analy_sol
u3 = u3_filter * analy_sol


# First order derivatives
p_px = kron(I_Nz,I_Ny,D1x)
p_py = kron(I_Nz,D1y,I_Nx)
p_pz = kron(D1z,I_Ny,I_Nx)

# Second order derivatives
p2_px2 = kron(I_Nz,I_Ny,D2x)
p2_py2 = kron(I_Nz,D2y,I_Nx)
p2_pz2 = kron(D2z,I_Ny,I_Nx)

# crossterms

p2_pypx = kron(I_Nz,D1y,D1x) # equivalent to p_py * p_px
p2_pxpy = kron(I_Nz,D1y,D1x)

p2_pzpy = kron(D1z,D1y,I_Nx)
p2_pypz = kron(D1z,D1y,I_Nx)

p2_pxpz = kron(D1z,I_Ny,D1x)
p2_pzpx = kron(D1z,I_Ny,D1x)

# Express σ tensors as operators on u vector (u1, u2, u3 stacked)
sigma_11 = (K - 2/3*μ) * (p_px * u1_filter + p_py * u2_filter + p_pz * u3_filter) + 2 * μ * p_px * u1_filter
sigma_12 = μ*(p_py*u1_filter + p_px*u2_filter)
sigma_13 = μ*(p_pz*u1_filter + p_px*u3_filter)


sigma_21 = μ*(p_px*u2_filter + p_py*u1_filter) 
sigma_22 = (K - 2/3*μ) * (p_px * u1_filter + p_py * u2_filter + p_pz * u3_filter) + 2 * μ * p_py * u2_filter
sigma_23 = μ * (p_pz * u2_filter + p_py * u3_filter)

sigma_31 = μ * (p_px * u3_filter + p_pz * u1_filter)
sigma_32 = μ * (p_py * u3_filter + p_pz * u2_filter)
sigma_33 = (K - 2/3 * μ) * (p_px * u1_filter + p_py * u2_filter + p_pz * u3_filter) + 2 * μ * p_pz * u3_filter


# Deriving equation for u1 as operators on u vector (u1, u2, u3 stacked)
u1_operator = p_px * sigma_11 + p_py * sigma_12 + p_pz * sigma_13

# Deriving equation for u2 as operators on u vector (u1, u2, u3 stacked)
u2_operator = p_px * sigma_21 + p_py * sigma_22 + p_pz * sigma_13

# Deriving equation for u3 as operators on u vector (u1, u2, u3 stacked)
u3_operator = p_px * sigma_31 + p_py * sigma_32 + p_pz * sigma_33


# Assembling source source_terms

source_u1 = u1_filter' * u1_operator * analy_sol
source_u2 = u2_filter' * u2_operator * analy_sol
source_u3 = u3_filter' * u3_operator * analy_sol

# Assembling boundary conditions
# Getting boundary values

# u1
u1_Front = Front_operator' * u_Front(y,z)[:] # Dirichlet Conditions
u1_End = End_operator' * u_End(y,z)[:] # Dirichlet Conditions

u1_Top = Top_operator' * u_Top(x,y)[:] # Dirichlet Conditions
u1_Bottom = Bottom_operator' * u_Bottom(x,y)[:] # Dirichlet Conditions

u1_Left = Left_operator' * u_y_Left(x,z)[:] # Neumann Conditions
u1_Right = Right_operator' * u_y_Right(x,z)[:]

# u2
u2_Front = Front_operator' * u_Front(y,z)[:] # Dirichlet Conditions
u2_End = End_operator' * u_End(y,z)[:] # Dirichlet Conditions

u2_Top = Top_operator' * u_Top(x,y)[:] # Dirichlet Conditions
u2_Bottom = Bottom_operator' * u_Bottom(x,y)[:] # Dirichlet Conditions

u2_Left = Left_operator' * u_y_Left(x,z)[:] # Neumann Conditions
u2_Right = Right_operator' * u_y_Right(x,z)[:]

# u3
u3_Front = Front_operator' * u_Front(y,z)[:] # Dirichlet Conditions
u3_End = End_operator' * u_End(y,z)[:] # Dirichlet Conditions

u3_Top = Top_operator' * u_Top(x,y)[:] # Dirichlet Conditions
u3_Bottom = Bottom_operator' * u_Bottom(x,y)[:] # Dirichlet Conditions

u3_Left = Left_operator' * u_y_Left(x,z)[:] # Neumann Conditions
u3_Right = Right_operator' * u_y_Right(x,z)[:]



# Assembling left hand side

A1 = (u1_filter' * u1_operator)

A2 = (u2_filter' * u2_operator)

A3 = (u3_filter' * u3_operator)

A = A1 + A2 + A3

# Assembling right hand side
source = source_u1 + source_u2 + source_u3