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


# SBP-SAT Terms
SAT_Left = tau_y * HI_Left * Left_operator'* Left_operator + beta * HI_Left * BS_Left' * Left_operator'*Left_operator
SAT_Right = tau_y * HI_Right * Right_operator'* Right_operator + beta * HI_Right * BS_Right' * Right_operator'*Right_operator

SAT_Top = tau_z * HI_Top  * Top_operator' * Top_operator * BS_Top
SAT_Bottom = tau_z * HI_Bottom  * Bottom_operator' * Bottom_operator * BS_Bottom

SAT_Front = tau_x * HI_Top * Front_operator' * Front_operator * BS_Front
SAT_End = tau_x * HI_Bottom * End_operator' * End_operator * BS_End


## To Be Completed
SAT_Left_r = tau_y * HI_Left * Left_operator' + beta * HI_Left * BS_Left' * Left_operator'
SAT_Right_r = tau_y * HI_Right * Right_operator'+ beta * HI_Right * BS_Right' * Right_operator'

SAT_Top_r = tau_z * HI_Top  * Top_operator' 
SAT_Bottom_r = tau_z * HI_Bottom * Bottom_operator'

SAT_Front_r = tau_x * HI_Top * Front_operator'
SAT_End_r = tau_x * HI_Bottom * End_operator'



# Analytical solutions
(analy_sol_3D,x_ex,y_ex,z_ex) = form_analy_sol(;N=N_x)
# plot(x_ex,y_ex,analy_sol_3D[:,:,end],st=:surface)


G_Left = [sin(π*i  + π*k) for i ∈ x_ex, k ∈ z_ex][:] # Dirichlet
G_Right = [-sin(π*i  + π*k) for i ∈ x_ex, k ∈ z_ex][:] # Dirichlet
G_Bottom = [-π*cos(π*i + π*j) for i ∈ x_ex, j ∈ y_ex][:] # norm vector -1 Neumann
G_Top = [π*cos(π*i + π*j + π) for i ∈ x_ex, j ∈ y_ex][:] # sin(pi) = -1 Neumann
G_Front = [π*cos(π*j + π*k + π) for j ∈ y_ex, k ∈ z_ex][:] # Neumann
G_End = [-π*cos(π*j + π*k) for j ∈ y_ex, k ∈ z_ex][:] # Neumann

source_terms = 3π^2 * analy_sol_3D;


# A = -(D2_x + D2_y + D2_z) + SAT_Front + SAT_End + SAT_Left + SAT_Right + SAT_Bottom + SAT_Top

# b = source_terms[:] + SAT_Front_r * G_Front + SAT_End_r * G_End + SAT_Left_r * G_Left + SAT_Right_r * G_Right + SAT_Bottom_r * G_Bottom + SAT_Top_r * G_Top


# numerical_sol = A\b

# numerical_sol_3D = reshape(numerical_sol,N_x+1,N_y+1,N_z+1)
# analy_sol_3D

u1_filter = get_u1(Nx,Ny,Nz)
u2_filter = get_u2(Nx,Ny,Nz)
u3_filter = get_u3(Nx,Ny,Nz)


numerical_sol = zeros(3*Nx*Ny*Nz)

for i in eachindex(numerical_sol)
    numerical_sol[i] = rem(i,3)
end

u1 = u1_filter * numerical_sol
u2 = u2_filter * numerical_sol
u3 = u3_filter * numerical_sol


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


