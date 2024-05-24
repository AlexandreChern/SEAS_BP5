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


############################### N = 64 ######################################
clear_mg_struct_CUDA(mg_struct_CUDA)
initialize_mg_struct_CUDA(mg_struct_CUDA, 64, 64, 64, 6)
get_lams(mg_struct_CUDA)

f_in = mg_struct_CUDA.b_mg[1]
mg_solver_CUDA(mg_struct_CUDA, nx = 64, ny = 64, nz=64, f_in; max_mg_iterations=10, n_levels=6, v1=10, v2 = 10, v3 = 10, print_results=true, scaling_factor=1, iter_algo_num=1)
mg_struct_CUDA.x_CUDA[1] .= 0
mgcg_CUDA(mg_struct_CUDA,nx=64,ny=64,nz=64,n_levels=6,precond=true,max_mg_iterations=1, v1=5, v2=10, v3=5, max_cg_iter=30,scaling_factor=1,print_results=true) 

mg_struct_CUDA.A_CPU_mg[1]

@benchmark for _ in 1:1
    mg_struct_CUDA.x_CUDA[1] .= 0
    mgcg_CUDA(mg_struct_CUDA,nx=64,ny=64,nz=64,n_levels=6,precond=true,max_mg_iterations=1, v1=5, v2=100, v3=5, max_cg_iter=30,scaling_factor=1, rel_tol=1e-6) 
end

@benchmark for _ in 1:1
    mg_struct_CUDA.x_CUDA[1] .= 0
    mgcg_CUDA(mg_struct_CUDA,nx=64,ny=64,nz=64,n_levels=6,precond=true,max_mg_iterations=1, v1=5, v2=100, v3=5, max_cg_iter=30,scaling_factor=1, rel_tol=1e-7) 
end



############################### N = 128 ######################################
clear_mg_struct_CUDA(mg_struct_CUDA)
initialize_mg_struct_CUDA(mg_struct_CUDA, 128, 128, 128, 7)
get_lams(mg_struct_CUDA)

f_in = mg_struct_CUDA.b_mg[1]
mg_solver_CUDA(mg_struct_CUDA, nx = 128, ny = 128, nz=128, f_in; max_mg_iterations=10, n_levels=7, v1=10, v2 = 100, v3 = 10, print_results=true, scaling_factor=1, iter_algo_num=1)
mg_struct_CUDA.x_CUDA[1] .= 0
mgcg_CUDA(mg_struct_CUDA,nx=128,ny=128,nz=128,n_levels=7,precond=true,max_mg_iterations=1, v1=5, v2=100, v3=5, max_cg_iter=30,scaling_factor=1) 

mg_struct_CUDA.A_CPU_mg[1]

@benchmark for _ in 1:1
    mg_struct_CUDA.x_CUDA[1] .= 0
    mgcg_CUDA(mg_struct_CUDA,nx=128,ny=128,nz=128,n_levels=7,precond=true,max_mg_iterations=1, v1=5, v2=100, v3=5, max_cg_iter=30,scaling_factor=1, rel_tol=1e-6) 
end

@benchmark for _ in 1:1
    mg_struct_CUDA.x_CUDA[1] .= 0
    mgcg_CUDA(mg_struct_CUDA,nx=128,ny=128,nz=128,n_levels=7,precond=true,max_mg_iterations=1, v1=5, v2=100, v3=5, max_cg_iter=30,scaling_factor=1, rel_tol=1e-7) 
end



############################### N = 256 ######################################
clear_mg_struct_CUDA(mg_struct_CUDA)
initialize_mg_struct_CUDA(mg_struct_CUDA, 256, 256, 256, 8)
get_lams(mg_struct_CUDA)

f_in = mg_struct_CUDA.b_mg[1]
mg_solver_CUDA(mg_struct_CUDA, nx = 256, ny = 256, nz=256, f_in; max_mg_iterations=10, n_levels=8, v1=10, v2 = 10, v3 = 10, print_results=true, scaling_factor=1, iter_algo_num=1)
mg_struct_CUDA.x_CUDA[1] .= 0
mgcg_CUDA(mg_struct_CUDA,nx=256,ny=256,nz=256,n_levels=8,precond=false,max_mg_iterations=1, v1=5, v2=5, v3=5, max_cg_iter=3000,scaling_factor=1) 

mg_struct_CUDA.x_CUDA[1] .= 0
mgcg_CUDA(mg_struct_CUDA,nx=256,ny=256,nz=256,n_levels=8,precond=true,max_mg_iterations=1, v1=5, v2=5, v3=5, max_cg_iter=30,scaling_factor=1) 


# julia> mgcg_CUDA(mg_struct_CUDA,nx=256,ny=256,nz=256,n_levels=8,precond=true,max_mg_iterations=1, v1=5, v2=5, v3=5, max_cg_iter=30,scaling_factor=1)        
# using H_tilde_2h * H_tilde_h^-1
# (k, norm_v_initial_norm) = (1, 0.006436951926747335)
# (k, norm_v_initial_norm) = (2, 0.0005254223510206546)
# (k, norm_v_initial_norm) = (3, 0.00012181443910793622)
# (k, norm_v_initial_norm) = (4, 4.9241457128124536e-5)
# (k, norm_v_initial_norm) = (5, 1.6505258274511717e-5)
# (k, norm_v_initial_norm) = (6, 5.1670865098935386e-6)
# (k, norm_v_initial_norm) = (7, 1.6681474624175163e-6)
# (k, norm_v_initial_norm) = (8, 5.167277942923694e-7)
# (k, norm_v_initial_norm) = (9, 1.4793698693336772e-7)
# (k, norm_v_initial_norm) = (10, 4.557867719399852e-8)
# (k, norm_v_initial_norm) = (11, 1.425339878856736e-8)

# julia> mgcg_CUDA(mg_struct_CUDA,nx=256,ny=256,nz=256,n_levels=8,precond=true,max_mg_iterations=1, v1=5, v2=5, v3=5, max_cg_iter=30,scaling_factor=1)     
# using constant scaling of 8
# (k, norm_v_initial_norm) = (1, 0.005335368218840409)
# (k, norm_v_initial_norm) = (2, 0.00194702294752371)
# (k, norm_v_initial_norm) = (3, 0.00018189868319323977)
# (k, norm_v_initial_norm) = (4, 6.832784342162066e-5) 
# (k, norm_v_initial_norm) = (5, 2.829042376293409e-5)
# (k, norm_v_initial_norm) = (6, 4.4056196113501565e-6) 
# (k, norm_v_initial_norm) = (7, 2.460823014646125e-6) 
# (k, norm_v_initial_norm) = (8, 7.026687034884572e-7) 
# (k, norm_v_initial_norm) = (9, 1.609102149369374e-7) 
# (k, norm_v_initial_norm) = (10, 9.305216028187453e-8)
# (k, norm_v_initial_norm) = (11, 3.1782700282913165e-8)
# (k, norm_v_initial_norm) = (12, 9.446818769670756e-9)  

mg_struct_CUDA.A_CPU_mg[1]




# exploring interpolation operators

N = 64
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


source_h_u1_reshaped = reshape(u1_filter_MF(source_h), N_xh + 1, N_yh + 1, N_zh + 1)
source_2h_u1_reshaped = reshape(u1_filter_MF(source_2h), N_x2h + 1, N_y2h + 1, N_z2h + 1)
source_h_u1_reshaped[9,:,:]
source_2h_u1_reshaped[5,:,:]

plot(0:hx_h:1,0:hy_h:1, source_h_u1_reshaped[1,:,:], st=:surface)
plot(0:hx_2h:1,0:hy_2h:1, source_2h_u1_reshaped[1,:,:], st=:surface)



mg_struct_CUDA.λ_mins
mg_struct_CUDA.λ_maxs

extrema(eigvals(Matrix(mg_struct_CUDA.A_CPU_mg[end])))
extrema(eigvals(Matrix(mg_struct_CUDA.A_CPU_mg[end-1])))
extrema(eigvals(Matrix(mg_struct_CUDA.A_CPU_mg[end-2])))
extrema(eigvals(Matrix(mg_struct_CUDA.A_CPU_mg[end-3])))