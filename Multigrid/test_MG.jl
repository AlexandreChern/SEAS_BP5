include("Assembling_3D_matrices.jl")
include("utils.jl")

clear_mg_struct_CUDA(mg_struct_CUDA)
initialize_mg_struct_CUDA(mg_struct_CUDA, 32, 32, 32, 5)

# u_direct_1 = mg_struct_CUDA.A_CPU_mg[1] \ Array(mg_struct_CUDA.b_mg[1])
# extrema((u_direct_1 - mg_struct_CUDA.u_exact[1]))


get_lams(mg_struct_CUDA)


f_in = mg_struct_CUDA.b_mg[1]

mg_solver_CUDA(mg_struct_CUDA, f_in; max_mg_iterations=2, n_levels=3)

mgcg_CUDA(mg_struct_CUDA,nx=32,ny=32,nz=32,n_levels=3,precond=true)


dot(mg_struct_CUDA.r_CUDA[1],mg_struct_CUDA.r_CUDA[1]) / dot(mg_struct_CUDA.r_CUDA[1], mg_struct_CUDA.A_mg[1] * mg_struct_CUDA.r_CUDA[1])


dot(mg_struct_CUDA.b_mg[1],mg_struct_CUDA.b_mg[1]) / dot(mg_struct_CUDA.b_mg[1], mg_struct_CUDA.A_mg[1] * mg_struct_CUDA.b_mg[1])