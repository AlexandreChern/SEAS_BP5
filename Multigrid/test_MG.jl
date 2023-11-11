include("Assembling_3D_matrices.jl")
include("utils.jl")

clear_mg_struct_CUDA(mg_struct_CUDA)
initialize_mg_struct_CUDA(mg_struct_CUDA, 32, 32, 32, 5)

# u_direct_1 = mg_struct_CUDA.A_CPU_mg[1] \ Array(mg_struct_CUDA.b_mg[1])
# extrema((u_direct_1 - mg_struct_CUDA.u_exact[1]))


get_lams(mg_struct_CUDA)


f_in = mg_struct_CUDA.b_mg[1]

mg_solver_CUDA(mg_struct_CUDA, f_in; max_mg_iterations=5, n_levels=3)

mg_struct_CUDA.x_CUDA[1] .= 0
mgcg_CUDA(mg_struct_CUDA,nx=32,ny=32,nz=32,n_levels=4,precond=true,max_cg_iter=10)

dot(mg_struct_CUDA.r_CUDA[1],mg_struct_CUDA.r_CUDA[1]) / dot(mg_struct_CUDA.r_CUDA[1], mg_struct_CUDA.A_mg[1] * mg_struct_CUDA.r_CUDA[1])


dot(mg_struct_CUDA.b_mg[1],mg_struct_CUDA.b_mg[1]) / dot(mg_struct_CUDA.b_mg[1], mg_struct_CUDA.A_mg[1] * mg_struct_CUDA.b_mg[1])

mg_struct_CUDA





# exploring interpolation operators

N = 16
N_xh = N_yh = N_zh = N
hx_h = 1 / N_xh
hy_h = 1 / N_yh
hz_h = 1 / N_zh

N_2h = div(N,2)
N_x2h = N_y2h = N_z2h = N_2h
hx_2h = 1 / N_x2h
hy_2h = 1 / N_y2h
hz_2h = 1 / N_z2h



M_h, RHS_h, H_tilde_h, HI_tilde_h, analy_sol_h = Assembling_3D_matrices(N_xh, N_yh, N_zh;p=2)
M_2h, RHS_2h, H_tilde_2h, HI_tilde_2h, analy_sol_2h = Assembling_3D_matrices(N_x2h, N_y2h, N_z2h;p=2)


RHS_h_1_reshaped = reshape(u1_filter_MF(RHS_h), N_xh + 1, N_yh + 1, N_zh + 1)
RHS_2h_1_reshaped = reshape(u1_filter_MF(RHS_2h), N_x2h + 1, N_y2h + 1, N_z2h + 1)


RHS_h_1_reshaped = reshape(HI_tilde_h * u1_filter_MF(RHS_h), N_xh + 1, N_yh + 1, N_zh + 1)
RHS_2h_1_reshaped = reshape(HI_tilde_2h * u1_filter_MF(RHS_2h), N_x2h + 1, N_y2h + 1, N_z2h + 1)

RHS_h_1_reshaped[1,:,:]
RHS_2h_1_reshaped[1,:,:]

rest_h = restriction_matrix_v0(N_xh,N_yh,N_zh,N_x2h,N_y2h,N_z2h) 
RHS_restricted = H_tilde_2h * rest_h * HI_tilde_h * u1_filter_MF(RHS_h)

RHS_resricted_reshaped = reshape(RHS_restricted,N_x2h + 1, N_y2h + 1, N_z2h + 1) / 2
reshape(u1_filter_MF(RHS_2h), N_x2h + 1, N_y2h + 1, N_z2h + 1)

plot(0:hx_h:1,0:hy_h:1, RHS_h_1_reshaped[1,:,:], st=:surface)
plot(0:hx_2h:1,0:hy_2h:1, RHS_2h_1_reshaped[1,:,:], st=:surface)