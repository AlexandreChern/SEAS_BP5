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

error_norms = []
for level = 2:6
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
    I_Nx = eyes(N_x+1)
    I_Ny = eyes(N_y+1)
    I_Nz = eyes(N_z+1)

    # e_E = kron(e_Nx,I_Ny);
    # e_W = kron(e_1x,I_Ny);
    # e_S = kron(I_Nx,e_1y);
    # e_N = kron(I_Nx,e_Ny);

    # Forming 3D operators

    # D2_x = kron(D2x,I_Ny,I_Ny)
    # D2_y = kron(I_Ny,D2y,I_Ny)
    # D2_z = kron(I_Ny,I_Ny,D2z)
    D2_x = kron(I_Ny,I_Nz,D2x)
    D2_y = kron(I_Nz,D2y,I_Nx)
    D2_z = kron(D2z,I_Nx,I_Ny)

    # e_1x = sparse(e(1,N_x+1))
    # e_1y = sparse(e(1,N_y+1))
    # e_1z = sparse(e(1,N_z+1))

    # e_Nx = sparse(e(N_x+1,N_x+1))
    # e_Ny = sparse(e(N_y+1,N_y+1))
    # e_Nz = sparse(e(N_z+1,N_z+1))

    e_1x = e(1,N_x + 1)
    e_1y = e(1,N_y + 1)
    e_1z = e(1,N_z + 1)

    e_Nx = e(N_x + 1,N_x + 1)
    e_Ny = e(N_y + 1,N_y + 1)
    e_Nz = e(N_z + 1,N_z + 1)


    # HI_Front = kron(I_Nx,HIy,HIz)
    # HI_End = kron(I_Nx,HIy,HIz)
    # HI_Front = kron(HIy,HIz,I_Nx)
    # HI_End = kron(HIy,HIz,I_Nx)

    HI_Front = kron(I_Ny,I_Nz,HIx)
    HI_End = kron(I_Ny,I_Nz,HIx)

    # HI_Left = kron(HIz,I_Ny,HIx)
    # HI_Right = kron(HIz,I_Ny,HIx)
    # HI_Left = kron(HIz,I_Ny,HIx)
    # HI_Right = kron(HIz,I_Ny,HIx)
    HI_Left = kron(I_Nz,HIy,I_Nx)
    HI_Right = kron(I_Nz,HIy,I_Nx)

    # HI_Bottom = kron(HIx,HIy,I_Nz)
    # HI_Top = kron(HIx,HIy,I_Nz)
    # HI_Bottom = kron(I_Nz,HIx,HIy)
    # HI_Top = kron(I_Nz,HIx,HIy)

    HI_Bottom = kron(HIz,I_Nx,I_Ny)
    HI_Top = kron(HIz,I_Nx,I_Ny)

    # BS_Front = kron(BSx,I_Ny,I_Nz)
    # BS_End = kron(BSx,I_Ny,I_Nz)

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

    # source Terms
    source_terms = 3π^2 * analy_sol_3D;


    A = -(D2_x + D2_y + D2_z) + SAT_Front + SAT_End + SAT_Left + SAT_Right + SAT_Bottom + SAT_Top

    b = source_terms[:] + SAT_Front_r * G_Front + SAT_End_r * G_End + SAT_Left_r * G_Left + SAT_Right_r * G_Right + SAT_Bottom_r * G_Bottom + SAT_Top_r * G_Top


    numerical_sol = A\b

    numerical_sol_3D = reshape(numerical_sol,N_x+1,N_y+1,N_z+1)
    analy_sol_3D

    # plot(x_ex,y_ex,numerical_sol_3D[:,:,end],st=:surface)


    error_norm = norm((numerical_sol-analy_sol_3D[:])'*H_tilde*(numerical_sol-analy_sol_3D[:]))
    append!(error_norms,error_norm)

    # numerical_laplacian = -(D2_x + D2_y + D2_z)*analy_sol_3D[:]

    # extrema(numerical_sol_3D[:] - analy_sol_3D[:])

    # Front_LHS = SAT_Front * analy_sol_3D[:]
    # Front_RHS = SAT_Front_r * G_Front


    # Bottom_LHS = SAT_Bottom * analy_sol_3D[:]
    # Bottom_RHS = SAT_Bottom_r * G_Bottom

    # HI_Bottom * Bottom_operator'* Bottom_operator * BS_Bottom *  analy_sol_3D[:]
    # HI_Bottom * Bottom_operator' * G_Bottom
end

@show error_norms
@show log2.(error_norms)