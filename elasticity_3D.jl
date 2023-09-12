include("diagonal_sbp.jl")
include("3D_face.jl")
include("analy_sol.jl")
# include("analy_sol_2.jl")
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

x = collect(range(0,step=hx,1));
y = collect(range(0,step=hy,1));
z = collect(range(0,step=hz,1));

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
HI_tilde = kron(HIx, HIy, HIz)

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


# H and HI Operators for End and Front
# Size: N^2 by N^2
H_1 = kron(H1z,H1y)
HI_1 = kron(HIz,HIy)

H_2 = kron(H1z, H1y)
HI_2 = kron(HIz, HIy)

# H and HI operators for Left and Right
# Size: N^2 by N^2
H_3 = kron(H1z, H1x)
HI_3 = kron(HIz, HIx)

H_4 = kron(H1z, H1x)
HI_4 = kron(HIz, HIx)

# H and HI operators for Bottom and Top
# Size: N^2 by N^2

H_5 = kron(H1y, H1x)
HI_5 = kron(HIy,HIx)

H_6 = kron(H1y, H1x)
HI_6 = kron(HIy,HIx)

# BS operators for 6 faces
# Size: N^3 by N^3
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


beta = 1



# A = -(D2_x + D2_y + D2_z) + SAT_Front + SAT_End + SAT_Left + SAT_Right + SAT_Bottom + SAT_Top

# b = source_terms[:] + SAT_Front_r * G_Front + SAT_End_r * G_End + SAT_Left_r * G_Left + SAT_Right_r * G_Right + SAT_Bottom_r * G_Bottom + SAT_Top_r * G_Top


# numerical_sol = A\b

# numerical_sol_3D = reshape(numerical_sol,N_x+1,N_y+1,N_z+1)

# analy_sol_3D u1 = u2 = u3 = sin(πx + πy + πz)

u1_filter = get_u1(Nx,Ny,Nz)
u2_filter = get_u2(Nx,Ny,Nz)
u3_filter = get_u3(Nx,Ny,Nz)


# analy_sol = zeros(3*Nx*Ny*Nz)

# # for i in eachindex(numerical_sol)
# #     analy_sol[i] = rem(i,3)
# # end

# # # setting values for u1 u2 u3
# # analy_sol[1:3:end] = analy_sol_3D[:] 
# # analy_sol[2:3:end] = analy_sol_3D[:]
# # analy_sol[3:3:end] = analy_sol_3D[:] 

# u1 = u1_filter * analy_sol
# u2 = u2_filter * analy_sol
# u3 = u3_filter * analy_sol


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
sigma_11 = (K_v - 2/3*μ_v) * (p_px * u1_filter + p_py * u2_filter + p_pz * u3_filter) + 2 * μ_v * p_px * u1_filter
sigma_12 = μ_v*(p_py*u1_filter + p_px*u2_filter)
sigma_13 = μ_v*(p_pz*u1_filter + p_px*u3_filter)


sigma_21 = μ_v*(p_px*u2_filter + p_py*u1_filter) 
sigma_22 = (K_v - 2/3*μ_v) * (p_px * u1_filter + p_py * u2_filter + p_pz * u3_filter) + 2 * μ_v * p_py * u2_filter
sigma_23 = μ_v * (p_pz * u2_filter + p_py * u3_filter)

sigma_31 = μ_v * (p_px * u3_filter + p_pz * u1_filter)
sigma_32 = μ_v * (p_py * u3_filter + p_pz * u2_filter)
sigma_33 = (K_v - 2/3 * μ_v) * (p_px * u1_filter + p_py * u2_filter + p_pz * u3_filter) + 2 * μ_v * p_pz * u3_filter


# # # Deriving equation for u1 as operators on u vector (u1, u2, u3 stacked)
# u1_operator = p_px * sigma_11 + p_py * sigma_12 + p_pz * sigma_13 # should rewrite explicitly using p2_px2

# # # Deriving equation for u2 as operators on u vector (u1, u2, u3 stacked)
# u2_operator = p_px * sigma_21 + p_py * sigma_22 + p_pz * sigma_13

# # # Deriving equation for u3 as operators on u vector (u1, u2, u3 stacked)
# u3_operator = p_px * sigma_31 + p_py * sigma_32 + p_pz * sigma_33

# # TO DO: Need to rewrite this part 
u1_operator =  ( (K_v - 2/3 * μ_v) * (p2_px2 * u1_filter + p2_pxpy * u2_filter + p2_pxpz * u3_filter) 
            + 2 * μ_v * p2_px2 * u1_filter 
            + μ_v * (p2_py2 * u1_filter + p2_pxpy * u2_filter)
            + μ_v * (p2_pz2 * u1_filter + p2_pxpz * u3_filter)
)

u2_operator = ( μ_v * (p2_px2 * u2_filter + p2_pxpy * u1_filter)
            + (K_v - 2/3 * μ_v) * (p2_pxpy * u1_filter + p2_py2 * u2_filter + p2_pypz * u3_filter)
            + 2 * μ_v * p2_py2 * u2_filter
            + μ_v * (p2_pz2 * u2_filter + p2_pypz * u3_filter)
)

u3_operator = ( μ_v * (p2_px2 * u3_filter + p2_pxpz * u1_filter)
            + μ_v * (p2_py2 * u3_filter + p2_pypz * u2_filter)
            + (K_v - 2/3 * μ_v) * (p2_pxpz * u1_filter + p2_pypz * u2_filter + p2_pz2 * u3_filter)
            + 2 * μ_v * p2_pz2 * u3_filter
)






### Assembling SBP terms according to the note

### Assembling components of stress tensors and boundary operators
### Face 1
e_1 = End_operator'
e_1T = End_operator

T_11_1 = - (K_v + 4/3 * μ_v) * p_px #* u1_filter
T_12_1 = - (K_v - 2/3 * μ_v) * p_pz #* u2_filter # Not quite sure 
T_13_1 = - (K_v - 2/3 * μ_v) * p_py #* u3_filter

T_21_1 = - μ_v * p_px #* u1_filter
T_22_1 = - μ_v * p_py #* u2_filter
T_23_1 = 0

T_31_1 = - μ_v * p_pz #* u1_filter
T_32_1 = 0
T_33_1 = - μ_v * p_px #* u3_filter


## TO DO Fix Z values
Z_11_1 = (d * β / H1x[1]) * (K_v + 4/3 * μ_v + 2 * μ_v) * kron(I_Nz, I_Ny, I_Nx)#* u1_filter
Z_12_1 = (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, I_Nx) #* u2_filter
Z_13_1 = (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, I_Nx) #* u3_filter

Z_21_1 =  (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, I_Nx) #* u1_filter
Z_22_1 = (d * β / H1x[1]) * (K_v + 4/3 * μ_v + 2 * μ_v) * kron(I_Nz, I_Ny, I_Nx)#* u2_filter
Z_23_1 =  (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, I_Nx) #* u3_filter ## 0 ?

Z_31_1 = (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, I_Nx)#* u1_filter
Z_32_1 = (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, I_Nx)
Z_33_1 = (d * β / H1x[1]) * (K_v + 4/3 * μ_v + 2 * μ_v) * kron(I_Nz, I_Ny, I_Nx) #* u3_filter

# # Z version 2
# Z_11_1 = (d * β / H1x[1]) * (K_v + 4/3 * μ_v + 2 * μ_v) * kron(I_Nz, I_Ny, e_1x * e_1x')#* u1_filter
# Z_12_1 = (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, e_1x * e_1x') #* u2_filter
# Z_13_1 = (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, e_1x * e_1x') #* u3_filter

# Z_21_1 =  (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, e_1x * e_1x') #* u1_filter
# Z_22_1 = (d * β / H1x[1]) * (K_v + 4/3 * μ_v + 2 * μ_v) * kron(I_Nz, I_Ny, e_1x * e_1x')#* u2_filter
# Z_23_1 =  (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, e_1x * e_1x') #* u3_filter ## 0 ?

# Z_31_1 = (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, e_1x * e_1x')#* u1_filter
# Z_32_1 = (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, e_1x * e_1x')
# Z_33_1 = (d * β / H1x[1]) * (K_v + 4/3 * μ_v + 2 * μ_v) * kron(I_Nz, I_Ny, e_1x * e_1x') #* u3_filter


### Face 2
e_2 = Front_operator'
e_2T = Front_operator

T_11_2 = (K_v + 4/3)  * p_px #* u1_filter
T_12_2 = (K_v - 2/3 * μ_v) * p_pz #* u2_filter # Not quite sure 
T_13_2 = (K_v - 2/3 * μ_v) * p_py #* u3_filter

T_21_2 = μ_v * p_px #* u1_filter
T_22_2 = μ_v * p_py #* u2_filter
T_23_2 = 0

T_31_2 = μ_v * p_pz #* u1_filter
T_32_2 = 0
T_33_2 = μ_v * p_px #* u3_filter

## TO DO Fix Z values
Z_11_2 = (d * β / H1x[1]) * (K_v + 4/3 * μ_v + 2 * μ_v) * kron(I_Nz, I_Ny, I_Nx) #* u1_filter
Z_12_2 = (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, I_Nx) #* u2_filter
Z_13_2 = (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, I_Nx) #* u3_filter

Z_21_2 =  (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, I_Nx) #* u1_filter
Z_22_2 = (d * β / H1x[1]) * (K_v + 4/3 * μ_v + 2 * μ_v) * kron(I_Nz, I_Ny, I_Nx) #* u2_filter
Z_23_2 =  (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, I_Nx) #* u3_filter ## 0 ?

Z_31_2 = (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, I_Nx) #* u1_filter
Z_32_2 = (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, I_Nx)
Z_33_2 = (d * β / H1x[1]) * (K_v + 4/3 * μ_v + 2 * μ_v) * kron(I_Nz, I_Ny, I_Nx)#* u3_filter

# # Z version 2
# Z_11_2 = (d * β / H1x[1]) * (K_v + 4/3 * μ_v + 2 * μ_v) * kron(I_Nz, I_Ny, e_Nx * e_Nx') #* u1_filter
# Z_12_2 = (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, e_Nx * e_Nx') #* u2_filter
# Z_13_2 = (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, e_Nx * e_Nx') #* u3_filter

# Z_21_2 =  (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, e_Nx * e_Nx') #* u1_filter
# Z_22_2 = (d * β / H1x[1]) * (K_v + 4/3 * μ_v + 2 * μ_v) * kron(I_Nz, I_Ny, e_Nx * e_Nx') #* u2_filter
# Z_23_2 =  (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, e_Nx * e_Nx') #* u3_filter ## 0 ?

# Z_31_2 = (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, e_Nx * e_Nx') #* u1_filter
# Z_32_2 = (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, e_Nx * e_Nx')
# Z_33_2 = (d * β / H1x[1]) * (K_v + 4/3 * μ_v + 2 * μ_v) * kron(I_Nz, I_Ny, e_Nx * e_Nx')#* u3_filter




### Face 3
e_3 = Left_operator'
e_3T = Left_operator

T_11_3 = - μ_v * p_py #* u1_filter
T_12_3 = - μ_v * p_px #* u2_filter
T_13_3 = 0 # Face 1

T_21_3 = - (K_v - 2/3 * μ_v) * p_px #* u1_filter
T_22_3 = - (K_v + 4/3 * μ_v) * p_py #* u2_filter
T_23_3 = - (K_v - 2/3 * μ_v) * p_pz #* u3_filter

T_31_3 = 0
T_32_3 = - μ_v * p_pz #* u2_filter
T_33_3 = - μ_v * p_py #* u3_filter

### Face 4
e_4 = Right_operator'
e_4T = Right_operator

T_11_4 = μ_v * p_py #* u1_filter
T_12_4 = μ_v * p_px #* u2_filter
T_13_4 = 0

T_21_4 = (K_v - 2/3 * μ_v) * p_px #* u1_filter
T_22_4 = (K_v + 4/3 * μ_v) * p_py #* u2_filter
T_23_4 = (K_v - 2/3 * μ_v) * p_pz # * u3_filter

T_31_4 = 0
T_32_4 = μ_v * p_pz #* u2_filter
T_33_4 = μ_v * p_py #* u3_filter

### Face 5
e_5 = Bottom_operator'
e_5T = Bottom_operator

T_11_5 = - μ_v * p_pz #* u1_filter
T_12_5 = 0
T_13_5 = - μ_v * p_px #* u3_filter

T_21_5 = 0
T_22_5 = - μ_v * p_pz #* u2_filter
T_23_5 = - μ_v * p_py #* u3_filter

T_31_5 = - (K_v - 2/3 * μ_v) * p_px #* u1_filter
T_32_5 = - (K_v - 2/3 * μ_v) * p_py #* u2_filter
T_33_5 = - (K_v + 4/3 * μ_v) * p_pz #* u3_filter

### Face 6
e_6 = Top_operator'
e_6T = Top_operator

T_11_6 = μ_v * p_pz #* u1_filter
T_12_6 = 0
T_13_6 = μ_v * p_px #* u3_filter

T_21_6 = 0
T_22_6 = μ_v * p_pz #* u2_filter
T_23_6 = μ_v * p_py #* u3_filter

T_31_6 = (K_v - 2/3 * μ_v) * p_px #* u1_filter
T_32_6 = (K_v - 2/3 * μ_v) * p_py #* u2_filter
T_33_6 = (K_v + 4/3 * μ_v) * p_pz #* u3_filter

### Assembling SBP terms for left-hand-side (LHS) traction condition

SAT_1_LHS = - HI_tilde * (
        e_3 * H_3 * e_3T * (T_11_3 * u1_filter .+ T_12_3 * u2_filter .+ T_13_3 * u3_filter)
    +   e_4 * H_4 * e_4T * (T_11_4 * u1_filter .+ T_12_4 * u2_filter .+ T_13_4 * u3_filter)
    +   e_5 * H_5 * e_5T * (T_11_5 * u1_filter .+ T_12_5 * u2_filter .+ T_13_5 * u3_filter)
    +   e_6 * H_6 * e_6T * (T_11_6 * u1_filter .+ T_12_6 * u2_filter .+ T_13_6 * u3_filter)

    #     e_3 * e_3T * (T_11_3 * u1_filter .+ T_12_3 * u2_filter .+ T_13_3 * u3_filter)
    # +   e_4 * e_4T * (T_11_4 * u1_filter .+ T_12_4 * u2_filter .+ T_13_4 * u3_filter)
    # +   e_5 * e_5T * (T_11_5 * u1_filter .+ T_12_5 * u2_filter .+ T_13_5 * u3_filter)
    # +   e_6 * e_6T * (T_11_6 * u1_filter .+ T_12_6 * u2_filter .+ T_13_6 * u3_filter)

) 

SAT_2_LHS = - HI_tilde * (
        e_3 * H_3 * e_3T * (T_21_3 * u1_filter .+ T_22_3 * u2_filter .+ T_23_3 * u3_filter)
    +   e_4 * H_4 * e_4T * (T_21_4 * u1_filter .+ T_22_4 * u2_filter .+ T_23_4 * u3_filter)
    +   e_5 * H_5 * e_5T * (T_21_5 * u1_filter .+ T_22_5 * u2_filter .+ T_23_5 * u3_filter)
    +   e_6 * H_5 * e_6T * (T_21_6 * u1_filter .+ T_22_6 * u2_filter .+ T_23_6 * u3_filter)

    #     e_3 * e_3T * (T_21_3 * u1_filter .+ T_22_3 * u2_filter .+ T_23_3 * u3_filter)
    # +   e_4 * e_4T * (T_21_4 * u1_filter .+ T_22_4 * u2_filter .+ T_23_4 * u3_filter)
    # +   e_5 * e_5T * (T_21_5 * u1_filter .+ T_22_5 * u2_filter .+ T_23_5 * u3_filter)
    # +   e_6 * e_6T * (T_21_6 * u1_filter .+ T_22_6 * u2_filter .+ T_23_6 * u3_filter)

) 


SAT_3_LHS = - HI_tilde * (
        e_3 * H_3 * e_3T * (T_31_3 * u1_filter .+ T_32_3 * u2_filter .+ T_33_3 * u3_filter)
    +   e_4 * H_4 * e_4T * (T_31_4 * u1_filter .+ T_32_4 * u2_filter .+ T_33_4 * u3_filter)
    +   e_5 * H_5 * e_5T * (T_31_5 * u1_filter .+ T_32_5 * u2_filter .+ T_33_5 * u3_filter)
    +   e_6 * H_6 * e_6T * (T_31_6 * u1_filter .+ T_32_6 * u2_filter .+ T_33_6 * u3_filter)

    #     e_3 * e_3T * (T_31_3 * u1_filter .+ T_32_3 * u2_filter .+ T_33_3 * u3_filter)
    # +   e_4 * e_4T * (T_31_4 * u1_filter .+ T_32_4 * u2_filter .+ T_33_4 * u3_filter)
    # +   e_5 * e_5T * (T_31_5 * u1_filter .+ T_32_5 * u2_filter .+ T_33_5 * u3_filter)
    # +   e_6 * e_6T * (T_31_6 * u1_filter .+ T_32_6 * u2_filter .+ T_33_6 * u3_filter)
)



### Assembling SBP terms for left-hand-side (LHS) Dirichlet condition

SAT_tilde_1_LHS =  HI_tilde * (
        (T_11_1' .- Z_11_1') * (e_1 * H_1 * (e_1T)) * u1_filter
    +   (T_21_1' .- Z_21_1') * (e_1 * H_1 * (e_1T)) * u2_filter
    +   (T_31_1' .- Z_31_1') * (e_1 * H_1 * (e_1T)) * u3_filter
    +   (T_11_2' .- Z_11_2') * (e_2 * H_2 * (e_2T)) * u1_filter
    +   (T_21_2' .- Z_21_2') * (e_2 * H_2 * (e_2T)) * u2_filter
    +   (T_31_2' .- Z_31_2') * (e_2 * H_2 * (e_2T)) * u3_filter
)

SAT_tilde_2_LHS =  HI_tilde * (
        (T_12_1' .- Z_12_1') * (e_1 * H_1 * (e_1T)) * u1_filter
    +   (T_22_1' .- Z_22_1') * (e_1 * H_1 * (e_1T)) * u2_filter
    +   (T_32_1' .- Z_32_1') * (e_1 * H_1 * (e_1T)) * u3_filter
    +   (T_12_2' .- Z_12_2') * (e_2 * H_2 * (e_2T)) * u1_filter
    +   (T_22_2' .- Z_22_2') * (e_2 * H_2 * (e_2T)) * u2_filter
    +   (T_32_2' .- Z_32_2') * (e_2 * H_2 * (e_2T)) * u3_filter
)

SAT_tilde_3_LHS =  HI_tilde * (
        (T_13_1' .- Z_13_1') * (e_1 * H_1 * (e_1T)) * u1_filter
    +   (T_23_1' .- Z_23_1') * (e_1 * H_1 * (e_1T)) * u2_filter
    +   (T_33_1' .- Z_33_1') * (e_1 * H_1 * (e_1T)) * u3_filter
    +   (T_13_2' .- Z_13_2') * (e_2 * H_2 * (e_2T)) * u1_filter
    +   (T_23_2' .- Z_23_2') * (e_2 * H_2 * (e_2T)) * u2_filter
    +   (T_33_2' .- Z_33_2') * (e_2 * H_2 * (e_2T)) * u3_filter
)


# New formulation of SAT_tilde_LHS version 2

# SAT_tilde_1_LHS =  HI_tilde * (
#         (e_1 * H_1 * (e_1T)) * (T_11_1' .- Z_11_1') * u1_filter
#     +   (e_1 * H_1 * (e_1T)) * (T_21_1' .- Z_21_1') * u2_filter
#     +   (e_1 * H_1 * (e_1T)) * (T_31_1' .- Z_31_1') * u3_filter
#     +   (e_2 * H_2 * (e_2T)) * (T_11_2' .- Z_11_2') * u1_filter
#     +   (e_2 * H_2 * (e_2T)) * (T_21_2' .- Z_21_2') * u2_filter
#     +   (e_2 * H_2 * (e_2T)) * (T_31_2' .- Z_31_2') * u3_filter
# )

# SAT_tilde_2_LHS =  HI_tilde * (
#         (e_1 * H_1 * (e_1T)) * (T_12_1' .- Z_12_1') * u1_filter
#     +   (e_1 * H_1 * (e_1T)) * (T_22_1' .- Z_22_1') * u2_filter
#     +   (e_1 * H_1 * (e_1T)) * (T_32_1' .- Z_32_1') * u3_filter
#     +   (e_2 * H_2 * (e_2T)) * (T_12_2' .- Z_12_2') * u1_filter
#     +   (e_2 * H_2 * (e_2T)) * (T_22_2' .- Z_22_2') * u2_filter
#     +   (e_2 * H_2 * (e_2T)) * (T_32_2' .- Z_32_2') * u3_filter
# )

# SAT_tilde_3_LHS =  HI_tilde * (
#         (e_1 * H_1 * (e_1T)) * (T_13_1' .- Z_13_1') * u1_filter
#     +   (e_1 * H_1 * (e_1T)) * (T_23_1' .- Z_23_1') * u2_filter
#     +   (e_1 * H_1 * (e_1T)) * (T_33_1' .- Z_33_1') * u3_filter
#     +   (e_2 * H_2 * (e_2T)) * (T_13_2' .- Z_13_2') * u1_filter
#     +   (e_2 * H_2 * (e_2T)) * (T_23_2' .- Z_23_2') * u2_filter
#     +   (e_2 * H_2 * (e_2T)) * (T_33_2' .- Z_33_2') * u3_filter
# )

# New formulation for SAT_tilde_LHS version 3
# SAT_tilde_1_LHS =  HI_tilde * (
#         (T_11_1' .- Z_11_1') * (e_1 * (e_1T)) * u1_filter
#     +   (T_21_1' .- Z_21_1') * (e_1 * (e_1T)) * u2_filter
#     +   (T_31_1' .- Z_31_1') * (e_1 * (e_1T)) * u3_filter
#     +   (T_11_2' .- Z_11_2') * (e_2 * (e_2T)) * u1_filter
#     +   (T_21_2' .- Z_21_2') * (e_2 * (e_2T)) * u2_filter
#     +   (T_31_2' .- Z_31_2') * (e_2 * (e_2T)) * u3_filter
# )

# SAT_tilde_2_LHS =  HI_tilde * (
#         (T_12_1' .- Z_12_1') * (e_1 * (e_1T)) * u1_filter
#     +   (T_22_1' .- Z_22_1') * (e_1 * (e_1T)) * u2_filter
#     +   (T_32_1' .- Z_32_1') * (e_1 * (e_1T)) * u3_filter
#     +   (T_12_2' .- Z_12_2') * (e_2 * (e_2T)) * u1_filter
#     +   (T_22_2' .- Z_22_2') * (e_2 * (e_2T)) * u2_filter
#     +   (T_32_2' .- Z_32_2') * (e_2 * (e_2T)) * u3_filter
# )

# SAT_tilde_3_LHS =  HI_tilde * (
#         (T_13_1' .- Z_13_1') * (e_1 * (e_1T)) * u1_filter
#     +   (T_23_1' .- Z_23_1') * (e_1 * (e_1T)) * u2_filter
#     +   (T_33_1' .- Z_33_1') * (e_1 * (e_1T)) * u3_filter
#     +   (T_13_2' .- Z_13_2') * (e_2 * (e_2T)) * u1_filter
#     +   (T_23_2' .- Z_23_2') * (e_2 * (e_2T)) * u2_filter
#     +   (T_33_2' .- Z_33_2') * (e_2 * (e_2T)) * u3_filter
# )

# Forming analytical solutions
u1 = form_analy_sol(;N = N_x)[1][:] # u1 is the only non-zero component
u2 = form_analy_sol(;N = N_x)[2][:] # u2 = 0 for the test case
u3 = form_analy_sol(;N = N_x)[3][:] # u3 = 0 for the test case

analy_sol = u1_filter' * u1 + u2_filter' * u2 + u3_filter' * u3


# SAT_1_LHS * u_analy
# SAT_2_LHS * u_analy
# SAT_3_LHS * u_analy

# SAT_tilde_1_LHS * u_analy
# SAT_tilde_2_LHS * u_analy
# SAT_tilde_3_LHS * u_analy




# Assembling boundary conditions
# Getting boundary values

# u1
u1_Front_value = Front_operator' * u1_Front(y,z)[:] # Dirichlet Conditions
u1_End_value = End_operator' * u1_End(y,z)[:] # Dirichlet Conditions

u1_Top_value = Top_operator' * u1_Top(x,y)[:] # Dirichlet Conditions
u1_Bottom_value = Bottom_operator' * u1_Top(x,y)[:] # Dirichlet Conditions

u1_Left_value = Left_operator' * u1_y_Left(x,z)[:] # Neumann Conditions
u1_Right_value = Right_operator' * u1_y_Right(x,z)[:]

# u2
u2_Front_value = Front_operator' * u1_Front(y,z)[:] # Dirichlet Conditions
u2_End_value = End_operator' * u1_End(y,z)[:] # Dirichlet Conditions

u2_Top_value = Top_operator' * u1_Top(x,y)[:] # Dirichlet Conditions
u2_Bottom_value = Bottom_operator' * u1_Top(x,y)[:] # Dirichlet Conditions

u2_Left_value = Left_operator' * u1_y_Left(x,z)[:] # Neumann Conditions
u2_Right_value = Right_operator' * u1_y_Right(x,z)[:]

# u3
u3_Front_value = Front_operator' * u1_Front(y,z)[:] # Dirichlet Conditions
u3_End_value = End_operator' * u1_End(y,z)[:] # Dirichlet Conditions

u3_Top_value = Top_operator' * u1_Top(x,y)[:] # Dirichlet Conditions
u3_Bottom_value = Bottom_operator' * u1_Top(x,y)[:] # Dirichlet Conditions

u3_Left_value = Left_operator' * u1_y_Left(x,z)[:] # Neumann Conditions
u3_Right_value = Right_operator' * u1_y_Right(x,z)[:]



# # Assembling left hand side

E1 = (u1_filter' * H_tilde * u1_operator)

E2 = (u2_filter' * H_tilde * u2_operator)

E3 = (u3_filter' * H_tilde * u3_operator)

E = ( E1 + E2 + E3)

# Assembling right hand side
# Assembling source source_terms
source_u1 = u1_filter' * H_tilde * ((K_v - 2/3 * μ_v) * (-π^2 * u1_analy(x,y,z)[:] -π^2 * u2_analy(x,y,z)[:]) 
            + 2 * μ_v * (-π^2 * u1_analy(x,y,z)[:])
            + μ_v * (-π^2 * u1_analy(x,y,z)[:] -π^2 * u2_analy(x,y,z)[:]) 
            + μ_v * (-π^2 * u1_analy(x,y,z)[:])
            )

source_u2 = u2_filter' * H_tilde * (μ_v * (-π^2 * u1_analy(x,y,z)[:] - π^2 * u2_analy(x,y,z)[:])
            + (K_v - 2/3 * μ_v) * (-π^2 * u1_analy(x,y,z)[:] -π^2 * u2_analy(x,y,z)[:] ) 
            + 2 * μ_v * (-π^2 * u2_analy(x,y,z)[:])
            + μ_v * (-π^2 * u2_analy(x,y,z)[:])
            )

source_u3 = u3_filter' * H_tilde * (μ_v * (-π^2 * u1_analy(x,y,z)[:])
            + μ_v * (-π^2 * u2_analy(x,y,z)[:])
            + (K_v - 2/3 * μ_v) * (-π^2 * u1_analy(x,y,z)[:] + -π^2 * u2_analy(x,y,z)[:])
            # u3 is set to be zero 
            )
source = ( source_u1 + source_u2 + source_u3)

# Assembling boundary data
# Face 1: Dirichlet
g₁¹ = u1_End(y,z)
g₂¹ = u2_End(y,z)
g₃¹ = zeros(Ny,Nz)

# Face 2: Dirichlet
g₁² = u1_Front(y,z)
g₂² = u2_Front(y,z)
g₃² = zeros(Ny,Nz)

# Face 3: Neumann
g₁³ = u1_y_Left(x,z) + u2_x_Left(x,z)
g₂³ = (K_v - 2/3 * μ_v) * u1_x_Left(x,z) + (K_v + 4/3 * μ_v) * u2_y_Left(x,z)
g₃³ = zeros(Nx,Nz)

# Face 4: Neumann
g₁⁴ = u1_y_Right(x,z) + u2_x_Right(x,z)
g₂⁴ = (K_v - 2/3 * μ_v) * u1_x_Right(x,z) + (K_v + 4/3 * μ_v) * u2_y_Right(x,z)
g₃⁴ = zeros(Nx,Nz)

# Face 5: Neumann
g₁⁵ = u1_z_Bottom(x,y)
g₂⁵ = u2_z_Bottom(x,y)
g₃⁵ = (K_v - 2/3 * μ_v) * (u1_x_Bottom(x,y) + u2_y_Bottom(x,y))

# Face 6: Neumann
g₁⁶ = u1_z_Top(x,y)
g₂⁶ = u2_z_Top(x,y)
g₃⁶ = (K_v - 2/3 * μ_v) * (u1_x_Top(x,y) + u2_y_Top(x,y))


### Assembling SBP terms for right-hand-side (RHS) traction condition
SAT_1_RHS = - HI_tilde * (
        e_3 * H_3 * g₁³[:]
    +   e_4 * H_4 * g₁⁴[:]
    +   e_5 * H_5 * g₁⁵[:]
    +   e_6 * H_6 * g₁⁶[:]
    #     e_3 * g₁³[:]
    # +   e_4 * g₁⁴[:]
    # +   e_5 * g₁⁵[:]
    # +   e_6 * g₁⁶[:]
)

SAT_2_RHS = - HI_tilde * (
        e_3 * H_3 * g₂³[:]
    +   e_4 * H_4 * g₂⁴[:]
    +   e_5 * H_5 * g₂⁵[:]
    +   e_6 * H_6 * g₂⁶[:]
    #     e_3 * g₂³[:]
    # +   e_4 * g₂⁴[:]
    # +   e_5 * g₂⁵[:]
    # +   e_6 * g₂⁶[:]
)

SAT_3_RHS = - HI_tilde * (
        e_3 * H_3 * g₃³[:]
    +   e_4 * H_4 * g₃⁴[:]
    +   e_5 * H_5 * g₃⁵[:]
    +   e_6 * H_6 * g₃⁶[:]
    #     e_3 * g₃³[:]
    # +   e_4 * g₃⁴[:]
    # +   e_5 * g₃⁵[:]
    # +   e_6 * g₃⁶[:]
)



### Assembling SBP terms for right-hand-side (RHS) Dirichlet condition
SAT_tilde_1_RHS =  HI_tilde * (
        (T_11_1' .- Z_11_1') * (e_1 * H_1 * g₁¹[:])
    +   (T_21_1' .- Z_21_1') * (e_1 * H_1 * g₂¹[:])
    +   (T_31_1' .- Z_31_1') * (e_1 * H_1 * g₃¹[:])
    +   (T_11_2' .- Z_11_2') * (e_2 * H_2 * g₁²[:])
    +   (T_21_2' .- Z_21_2') * (e_2 * H_2 * g₂²[:])
    +   (T_31_2' .- Z_31_2') * (e_2 * H_2 * g₃²[:])
)

SAT_tilde_2_RHS =  HI_tilde * (
        (T_12_1' .- Z_12_1') * (e_1 * H_1 * g₁¹[:])
    +   (T_22_1' .- Z_22_1') * (e_1 * H_1 * g₂¹[:])
    +   (T_32_1' .- Z_32_1') * (e_1 * H_1 * g₃¹[:])
    +   (T_12_2' .- Z_12_2') * (e_2 * H_2 * g₁²[:])
    +   (T_22_2' .- Z_22_2') * (e_2 * H_2 * g₂²[:])
    +   (T_32_2' .- Z_31_2') * (e_2 * H_2 * g₃²[:])
)

SAT_tilde_3_RHS =  HI_tilde * (
        (T_13_1' .- Z_13_1') * (e_1 * H_1 * g₁¹[:])
    +   (T_23_1' .- Z_23_1') * (e_1 * H_1 * g₂¹[:])
    +   (T_33_1' .- Z_33_1') * (e_1 * H_1 * g₃¹[:])
    +   (T_13_2' .- Z_13_2') * (e_2 * H_2 * g₁²[:])
    +   (T_23_2' .- Z_23_2') * (e_2 * H_2 * g₂²[:])
    +   (T_33_2' .- Z_33_2') * (e_2 * H_2 * g₃²[:])
)

# # New formulation for SAT_tilde_RHS version 2
# SAT_tilde_1_RHS =  HI_tilde * (
#         (e_1 * H_1 * (e_1T)) * (T_11_1' .- Z_11_1') * (e_1 * g₁¹[:])
#     +   (e_1 * H_1 * (e_1T)) * (T_21_1' .- Z_21_1') * (e_1 * g₂¹[:])
#     +   (e_1 * H_1 * (e_1T)) * (T_31_1' .- Z_31_1') * (e_1 * g₃¹[:])
#     +   (e_2 * H_2 * (e_2T)) * (T_11_2' .- Z_11_2') * (e_2 * g₁²[:])
#     +   (e_2 * H_2 * (e_2T)) * (T_21_2' .- Z_21_2') * (e_2 * g₂²[:])
#     +   (e_2 * H_2 * (e_2T)) * (T_31_2' .- Z_31_2') * (e_2 * g₃²[:])
# )

# SAT_tilde_2_RHS =  HI_tilde * (
#         (e_1 * H_1 * (e_1T)) * (T_12_1' .- Z_12_1') * (e_1 * g₁¹[:])
#     +   (e_1 * H_1 * (e_1T)) * (T_22_1' .- Z_22_1') * (e_1 * g₂¹[:])
#     +   (e_1 * H_1 * (e_1T)) * (T_32_1' .- Z_32_1') * (e_1 * g₃¹[:])
#     +   (e_2 * H_2 * (e_2T)) * (T_12_2' .- Z_12_2') * (e_2 * g₁²[:])
#     +   (e_2 * H_2 * (e_2T)) * (T_22_2' .- Z_22_2') * (e_2 * g₂²[:])
#     +   (e_2 * H_2 * (e_2T)) * (T_32_2' .- Z_31_2') * (e_2 * g₃²[:])
# )

# SAT_tilde_3_RHS =  HI_tilde * (
#         (e_1 * H_1 * (e_1T)) * (T_13_1' .- Z_13_1') * (e_1 * g₁¹[:])
#     +   (e_1 * H_1 * (e_1T)) * (T_23_1' .- Z_23_1') * (e_1 * g₂¹[:])
#     +   (e_1 * H_1 * (e_1T)) * (T_33_1' .- Z_33_1') * (e_1 * g₃¹[:])
#     +   (e_2 * H_2 * (e_2T)) * (T_13_2' .- Z_13_2') * (e_2 * g₁²[:])
#     +   (e_2 * H_2 * (e_2T)) * (T_23_2' .- Z_23_2') * (e_2 * g₂²[:])
#     +   (e_2 * H_2 * (e_2T)) * (T_33_2' .- Z_33_2') * (e_2 * g₃²[:])
# )

# # New formulation for SAT_tilde_RHS version 3
# SAT_tilde_1_RHS =  HI_tilde * (
#         (T_11_1' .- Z_11_1') * (e_1 * g₁¹[:])
#     +   (T_21_1' .- Z_21_1') * (e_1 * g₂¹[:])
#     +   (T_31_1' .- Z_31_1') * (e_1 * g₃¹[:])
#     +   (T_11_2' .- Z_11_2') * (e_2 * g₁²[:])
#     +   (T_21_2' .- Z_21_2') * (e_2 * g₂²[:])
#     +   (T_31_2' .- Z_31_2') * (e_2 * g₃²[:])
# )

# SAT_tilde_2_RHS =  HI_tilde * (
#         (T_12_1' .- Z_12_1') * (e_1 * g₁¹[:])
#     +   (T_22_1' .- Z_22_1') * (e_1 * g₂¹[:])
#     +   (T_32_1' .- Z_32_1') * (e_1 * g₃¹[:])
#     +   (T_12_2' .- Z_12_2') * (e_2 * g₁²[:])
#     +   (T_22_2' .- Z_22_2') * (e_2 * g₂²[:])
#     +   (T_32_2' .- Z_31_2') * (e_2 * g₃²[:])
# )

# SAT_tilde_3_RHS =  HI_tilde * (
#         (T_13_1' .- Z_13_1') * (e_1 * g₁¹[:])
#     +   (T_23_1' .- Z_23_1') * (e_1 * g₂¹[:])
#     +   (T_33_1' .- Z_33_1') * (e_1 * g₃¹[:])
#     +   (T_13_2' .- Z_13_2') * (e_2 * g₁²[:])
#     +   (T_23_2' .- Z_23_2') * (e_2 * g₂²[:])
#     +   (T_33_2' .- Z_33_2') * (e_2 * g₃²[:])
# )


# Assembling LHS of the linear system

M = (E + u1_filter' * H_tilde * SAT_1_LHS
        + u2_filter' * H_tilde * SAT_2_LHS 
        + u3_filter' * H_tilde * SAT_3_LHS 
        + u1_filter' * H_tilde * SAT_tilde_1_LHS 
        + u2_filter' * H_tilde * SAT_tilde_2_LHS 
        + u3_filter' * H_tilde * SAT_tilde_3_LHS
    )

RHS = (source + u1_filter' * H_tilde * SAT_1_RHS 
        + u2_filter' * H_tilde * SAT_2_RHS 
        + u3_filter' * H_tilde * SAT_3_RHS
        + u1_filter' * H_tilde * SAT_tilde_1_RHS 
        + u2_filter' * H_tilde * SAT_tilde_2_RHS 
        + u3_filter' * H_tilde * SAT_tilde_3_RHS)


u_direct = M \ RHS

error_direct = sqrt((u1_filter * u_direct - u1)' * H_tilde * (u1_filter * u_direct - u1) 
            + (u2_filter * u_direct - u2)' * H_tilde * (u2_filter * u_direct - u2) 
            + (u3_filter * u_direct - u3)' * H_tilde * (u3_filter * u_direct - u3) 
)

reshape((u1_filter * u_direct - u1), Nx, Ny, Nz)
reshape((u2_filter * u_direct - u2), Nx, Ny, Nz)
reshape((u3_filter * u_direct - u3), Nx, Ny, Nz)

sqrt((u1_filter * u_direct - u1)' * H_tilde * (u1_filter * u_direct - u1))
sqrt((u2_filter * u_direct - u2)' * H_tilde * (u2_filter * u_direct - u2))
sqrt((u3_filter * u_direct - u3)' * H_tilde * (u3_filter * u_direct - u3))