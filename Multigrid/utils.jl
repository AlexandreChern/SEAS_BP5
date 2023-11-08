using CUDA

mutable struct MG_CUDA
    A_mg
    b_mg
    A_CPU_mg
    λ_mins
    λ_maxs
    odata_mg
    H_mg
    H_inv_mg
    f_mg
    u_mg
    r_mg
    prol_fine_mg
    rest_mg
    prol_mg
    lnx_mg
    lny_mg
    u_exact
    discretization_error

    # struct for mgcg_CUDA
    x_CUDA
    r_CUDA
    r_new_CUDA
    z_CUDA
    z_new_CUDA
    p_CUDA
end

mg_struct_CUDA = MG_CUDA([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[])

function clear_mg_struct_CUDA(mg_struct_CUDA)
    println("Clearing mg_struct")
    mg_struct_CUDA.A_mg = []
    mg_struct_CUDA.b_mg = []
    mg_struct_CUDA.A_CPU_mg = []
    mg_struct_CUDA.λ_mins = []
    mg_struct_CUDA.λ_maxs = []
    mg_struct_CUDA.odata_mg = []
    mg_struct_CUDA.H_mg = []
    mg_struct_CUDA.H_inv_mg = []
    mg_struct_CUDA.f_mg = []
    mg_struct_CUDA.u_mg = []
    mg_struct_CUDA.r_mg = []
    mg_struct_CUDA.prol_fine_mg = []
    mg_struct_CUDA.rest_mg = []
    mg_struct_CUDA.prol_mg = []
    mg_struct_CUDA.lnx_mg = []
    mg_struct_CUDA.lny_mg = []
    mg_struct_CUDA.u_exact = []
    mg_struct_CUDA.discretization_error = []

    # struct for mgcg_CUDA
    mg_struct_CUDA.x_CUDA = []
    mg_struct_CUDA.r_CUDA = []
    mg_struct_CUDA.r_new_CUDA = []
    mg_struct_CUDA.z_CUDA = []
    mg_struct_CUDA.z_new_CUDA = []
    mg_struct_CUDA.p_CUDA = []
end


function clear_urf_CUDA(mg_struct_CUDA)
    for i in 1:length(mg_struct_CUDA.u_mg)
        mg_struct_CUDA.u_mg[i] .= 0
        mg_struct_CUDA.r_mg[i] .= 0
        mg_struct_CUDA.f_mg[i] .= 0
        mg_struct_CUDA.prol_fine_mg[i] .= 0
    end
end

function initialize_mg_struct_CUDA(mg_struct_CUDA, nx, ny, nz, n_level)
    A_mg = mg_struct_CUDA.A_mg
    b_mg = mg_struct_CUDA.b_mg
    A_CPU_mg = mg_struct_CUDA.A_CPU_mg
    odata_mg = mg_struct_CUDA.odata_mg
    H_mg = mg_struct_CUDA.H_mg
    H_inv_mg = mg_struct_CUDA.H_inv_mg
    f_mg = mg_struct_CUDA.f_mg
    u_mg = mg_struct_CUDA.u_mg
    r_mg = mg_struct_CUDA.r_mg
    prol_fine_mg = mg_struct_CUDA.prol_fine_mg
    rest_mg = mg_struct_CUDA.rest_mg
    prol_mg = mg_struct_CUDA.prol_mg
    lnx_mg = mg_struct_CUDA.lnx_mg
    lny_mg = mg_struct_CUDA.lny_mg
    u_exact = mg_struct_CUDA.u_exact

    if isempty(A_mg)
        for k in 1:n_level
            hx = 1/nx
            hy = 1/ny
            hz = 1/nz
            if k == 1
                A, b, H_tilde, HI_tilde = Assembling_3D_matrices(nx,ny,nz)
                push!(A_CPU_mg, A)
                push!(A_mg, CUDA.CUSPARSE.CuSparseMatrixCSR(A))
                push!(b_mg, CuArray(b))
                push!(H_mg, CUDA.CUSPARSE.CuSparseMatrixCSR(H_tilde))
                push!(H_inv_mg, CUDA.CUSPARSE.CuSparseMatrixCSR(HI_tilde))
            else
                A, b, H_tilde, HI_tilde = Assembling_3D_matrices(nx,ny,nz)
                push!(A_CPU_mg, A)
                push!(A_mg, CUDA.CUSPARSE.CuSparseMatrixCSR(A))
                push!(b_mg, CuArray(b))
                push!(H_mg, CUDA.CUSPARSE.CuSparseMatrixCSR(H_tilde))
                push!(H_inv_mg, CUDA.CUSPARSE.CuSparseMatrixCSR(HI_tilde))
            end
        end 
    end
    # # struct for mgcg_CUDA
    # mg_struct_CUDA.x_CUDA
    # mg_struct_CUDA.r_CUDA
    # mg_struct_CUDA.r_new_CUDA
    # mg_struct_CUDA.z_CUDA
    # mg_struct_CUDA.z_new_CUDA
    # mg_struct_CUDA.p_CUDA
end