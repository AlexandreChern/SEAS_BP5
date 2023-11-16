include("Assembling_3D_matrices.jl")
include("utils.jl")

clear_mg_struct_CUDA(mg_struct_CUDA)
initialize_mg_struct_CUDA(mg_struct_CUDA, 128, 128, 128, 6)

# u_direct_1 = mg_struct_CUDA.A_CPU_mg[1] \ Array(mg_struct_CUDA.b_mg[1])
# extrema((u_direct_1 - mg_struct_CUDA.u_exact[1]))


get_lams(mg_struct_CUDA)


f_in = mg_struct_CUDA.b_mg[1]

mg_solver_CUDA(mg_struct_CUDA, nx = 64, ny = 64, nz=64, f_in; max_mg_iterations=10, n_levels=2, v1=10, v2 = 100, v3 = 10, print_results=true, scaling_factor=1, iter_algo_num=1)

mg_struct_CUDA.x_CUDA[1] .= 0
mgcg_CUDA(mg_struct_CUDA,nx=64,ny=64,nz=64,n_levels=6,precond=true,max_mg_iterations=1, v1=5, v2=100, v3=5, max_cg_iter=30,scaling_factor=1) # check mgcg implementation! precond=false should give good convergence
x_out, history = cg(mg_struct_CUDA.A_mg[1], mg_struct_CUDA.b_mg[1], log=true)
history.data

dot(mg_struct_CUDA.r_CUDA[1],mg_struct_CUDA.r_CUDA[1]) / dot(mg_struct_CUDA.r_CUDA[1], mg_struct_CUDA.A_mg[1] * mg_struct_CUDA.r_CUDA[1])


dot(mg_struct_CUDA.b_mg[1],mg_struct_CUDA.b_mg[1]) / dot(mg_struct_CUDA.b_mg[1], mg_struct_CUDA.A_mg[1] * mg_struct_CUDA.b_mg[1])

mg_struct_CUDA



############################### N = 128 ######################################
clear_mg_struct_CUDA(mg_struct_CUDA)
initialize_mg_struct_CUDA(mg_struct_CUDA, 128, 128, 128, 7)
get_lams(mg_struct_CUDA)

f_in = mg_struct_CUDA.b_mg[1]
mg_solver_CUDA(mg_struct_CUDA, nx = 128, ny = 128, nz=128, f_in; max_mg_iterations=10, n_levels=7, v1=10, v2 = 100, v3 = 10, print_results=true, scaling_factor=1, iter_algo_num=1)
mg_struct_CUDA.x_CUDA[1] .= 0
mgcg_CUDA(mg_struct_CUDA,nx=128,ny=128,nz=128,n_levels=7,precond=true,max_mg_iterations=1, v1=5, v2=100, v3=5, max_cg_iter=30,scaling_factor=1) 

mg_struct_CUDA.A_CPU_mg[1]


clear_mg_struct_CUDA(mg_struct_CUDA)
initialize_mg_struct_CUDA(mg_struct_CUDA, 256, 256, 256, 8)
get_lams(mg_struct_CUDA)

f_in = mg_struct_CUDA.b_mg[1]
mg_solver_CUDA(mg_struct_CUDA, nx = 256, ny = 256, nz=256, f_in; max_mg_iterations=10, n_levels=8, v1=10, v2 = 10, v3 = 10, print_results=true, scaling_factor=1, iter_algo_num=1)
mg_struct_CUDA.x_CUDA[1] .= 0
mgcg_CUDA(mg_struct_CUDA,nx=256,ny=256,nz=256,n_levels=8,precond=false,max_mg_iterations=1, v1=5, v2=5, v3=5, max_cg_iter=3000,scaling_factor=1) 

mg_struct_CUDA.x_CUDA[1] .= 0
mgcg_CUDA(mg_struct_CUDA,nx=256,ny=256,nz=256,n_levels=8,precond=true,max_mg_iterations=1, v1=5, v2=5, v3=5, max_cg_iter=30,scaling_factor=1) 


mg_struct_CUDA.A_CPU_mg[1]




# exploring interpolation operators

N = 32
N_xh = N_yh = N_zh = N
hx_h = 1 / N_xh
hy_h = 1 / N_yh
hz_h = 1 / N_zh

N_2h = div(N,2)
N_x2h = N_y2h = N_z2h = N_2h
hx_2h = 1 / N_x2h
hy_2h = 1 / N_y2h
hz_2h = 1 / N_z2h



M_h, RHS_h, H_tilde_h, HI_tilde_h, analy_sol_h, source_h = Assembling_3D_matrices(N_xh, N_yh, N_zh;p=2)
M_2h, RHS_2h, H_tilde_2h, HI_tilde_2h, analy_sol_2h, source_2h = Assembling_3D_matrices(N_x2h, N_y2h, N_z2h;p=2)


RHS_h_1_reshaped = reshape(u1_filter_MF(RHS_h), N_xh + 1, N_yh + 1, N_zh + 1)
RHS_2h_1_reshaped = reshape(u1_filter_MF(RHS_2h), N_x2h + 1, N_y2h + 1, N_z2h + 1)


RHS_h_1_reshaped = reshape(HI_tilde_h * u1_filter_MF(RHS_h), N_xh + 1, N_yh + 1, N_zh + 1)
RHS_2h_1_reshaped = reshape(HI_tilde_2h * u1_filter_MF(RHS_2h), N_x2h + 1, N_y2h + 1, N_z2h + 1)
[]
RHS_h_1_reshaped[1,:,:]
RHS_2h_1_reshaped[1,:,:]

rest_h = restriction_matrix_v0(N_xh,N_yh,N_zh,N_x2h,N_y2h,N_z2h) 
RHS_restricted = H_tilde_2h * rest_h * HI_tilde_h * u1_filter_MF(RHS_h)

RHS_resricted_reshaped = reshape(RHS_restricted,N_x2h + 1, N_y2h + 1, N_z2h + 1) / 2
reshape(u1_filter_MF(RHS_2h), N_x2h + 1, N_y2h + 1, N_z2h + 1)

plot(0:hx_h:1,0:hy_h:1, RHS_h_1_reshaped[:,:,1], st=:surface)
plot(0:hx_2h:1,0:hy_2h:1, RHS_2h_1_reshaped[:,:,1], st=:surface)



source_h_u1_reshaped = reshape(HI_tilde_h * u1_filter_MF(source_h), N_xh + 1, N_yh + 1, N_zh + 1)
source_2h_u1_reshaped = reshape(HI_tilde_2h * u1_filter_MF(source_2h), N_x2h + 1, N_y2h + 1, N_z2h + 1)


plot(0:hx_h:1,0:hy_h:1, source_h_u1_reshaped[1,:,:], st=:surface)
plot(0:hx_2h:1,0:hy_2h:1, source_2h_u1_reshaped[1,:,:], st=:surface)



mg_struct_CUDA.λ_mins
mg_struct_CUDA.λ_maxs

extrema(eigvals(Matrix(mg_struct_CUDA.A_CPU_mg[end])))
extrema(eigvals(Matrix(mg_struct_CUDA.A_CPU_mg[end-1])))
extrema(eigvals(Matrix(mg_struct_CUDA.A_CPU_mg[end-2])))
extrema(eigvals(Matrix(mg_struct_CUDA.A_CPU_mg[end-3])))