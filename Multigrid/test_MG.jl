include("Assembling_3D_matrices.jl")
include("utils.jl")

clear_mg_struct_CUDA(mg_struct_CUDA)
initialize_mg_struct_CUDA(mg_struct_CUDA, 16, 16, 16, 4)

u_direct_1 = mg_struct_CUDA.A_CPU_mg[1] \ Array(mg_struct_CUDA.b_mg[1])
extrema((u_direct_1 - mg_struct_CUDA.u_exact[1]))

sqrt(
    (CuArray(u1_filter(u_direct_1)) - CuArray(u1_filter(mg_struct_CUDA.u_exact[1])))' 
    * mg_struct_CUDA.H_mg[1] * (CuArray(u1_filter(u_direct_1)) - CuArray(u1_filter(mg_struct_CUDA.u_exact[1])))
+   (CuArray(u2_filter(u_direct_1)) - CuArray(u2_filter(mg_struct_CUDA.u_exact[1])))' 
    * mg_struct_CUDA.H_mg[1] * (CuArray(u2_filter(u_direct_1)) - CuArray(u2_filter(mg_struct_CUDA.u_exact[1])))
+   (CuArray(u3_filter(u_direct_1)) - CuArray(u3_filter(mg_struct_CUDA.u_exact[1])))' 
    * mg_struct_CUDA.H_mg[1] * (CuArray(u3_filter(u_direct_1)) - CuArray(u3_filter(mg_struct_CUDA.u_exact[1])))
)

get_lams(mg_struct_CUDA)