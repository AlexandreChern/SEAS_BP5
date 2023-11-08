using CUDA
using Arpack

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
                A, b, H_tilde, HI_tilde, analy_sol = Assembling_3D_matrices(nx,ny,nz)
                push!(A_CPU_mg, A)
                push!(A_mg, CUDA.CUSPARSE.CuSparseMatrixCSR(A))
                push!(b_mg, CuArray(b))
                push!(H_mg, CUDA.CUSPARSE.CuSparseMatrixCSR(H_tilde))
                push!(H_inv_mg, CUDA.CUSPARSE.CuSparseMatrixCSR(HI_tilde))
                push!(u_exact, analy_sol)
            else
                A, b, H_tilde, HI_tilde, analy_sol = Assembling_3D_matrices(nx,ny,nz)
                push!(A_CPU_mg, A)
                push!(A_mg, CUDA.CUSPARSE.CuSparseMatrixCSR(A))
                push!(b_mg, CuArray(b))
                push!(H_mg, CUDA.CUSPARSE.CuSparseMatrixCSR(H_tilde))
                push!(H_inv_mg, CUDA.CUSPARSE.CuSparseMatrixCSR(HI_tilde))
                push!(u_exact, analy_sol)
            end
            nx, ny, nz = div(nx,2), div(ny,2), div(nz,2)
            hx, hy = 2*hx, 2*hy, 2*hz
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


function u1_filter(u)
    return u[1:3:end]
end

function u2_filter(u)
    return u[2:3:end]
end

function u3_filter(u)
    return u[3:3:end]
end


function discretization_error(k)
    u_direct_k = mg_struct_CUDA.A_CPU_mg[k] \ Array(mg_struct_CUDA.b_mg[k])
    # extrema((u_direct_k - mg_struct_CUDA.u_exact[k]))

    sqrt(
        (CuArray(u1_filter(u_direct_k)) - CuArray(u1_filter(mg_struct_CUDA.u_exact[k])))' 
        * mg_struct_CUDA.H_mg[k] * (CuArray(u1_filter(u_direct_k)) - CuArray(u1_filter(mg_struct_CUDA.u_exact[k])))
    +   (CuArray(u2_filter(u_direct_k)) - CuArray(u2_filter(mg_struct_CUDA.u_exact[k])))' 
        * mg_struct_CUDA.H_mg[k] * (CuArray(u2_filter(u_direct_k)) - CuArray(u2_filter(mg_struct_CUDA.u_exact[k])))
    +   (CuArray(u3_filter(u_direct_k)) - CuArray(u3_filter(mg_struct_CUDA.u_exact[k])))' 
        * mg_struct_CUDA.H_mg[k] * (CuArray(u3_filter(u_direct_k)) - CuArray(u3_filter(mg_struct_CUDA.u_exact[k])))
    )
end



function get_lams(mg_struct_CUDA)
    # TO DO, get 
    empty!(mg_struct_CUDA.λ_mins)
    empty!(mg_struct_CUDA.λ_maxs)
    # empty!(mg_struct.αs)
    # empty!(mg_struct.δs)
    reverse_Amg = reverse(mg_struct_CUDA.A_CPU_mg)
    if size(reverse_Amg[1])[1] > 289
        println("The minimal A matrix is too large for λ_min calculation")
        return 0
    end
    for k in eachindex(reverse_Amg)
        # lam_max, v_max = eigs(reverse_Amg[k], nev=1, which=:LR)
        # lam_max = real(lam_max[1]) # try different formulations
        if size(reverse_Amg[k])[1] <= 1089 # nx <= 32
        # if size(reverse_Amg[k])[1] <= 4225 # nx <= 64 not consistent, sometimes work
            lam_min, v_min = eigs(reverse_Amg[k], nev=1, which=:SR)
            lam_min = real(lam_min[1])
            # @show lam_min

            lam_max, v_max = eigs(reverse_Amg[k], nev=1, which=:LR)
            lam_max = real(lam_max[1]) # try different formulations
        else
            lam_min = mg_struct_CUDA.λ_mins[1] / 4
            # @show mg_struct_CUDA.λ_mins[1]
            # @show lam_min
            if size(reverse_Amg[k])[1] <= 16641
                lam_max, v_max = eigs(reverse_Amg[k], nev=1, which=:LR)
                lam_max = real(lam_max[1]) # try different formulations
            else
                lam_max = mg_struct_CUDA.λ_maxs[1] + (mg_struct_CUDA.λ_maxs[1] - mg_struct_CUDA.λ_maxs[2]) * 0.6
            end
            # @show mg_struct_CUDA.λ_maxs[1]
        end
        pushfirst!(mg_struct_CUDA.λ_mins, lam_min)
        pushfirst!(mg_struct_CUDA.λ_maxs, lam_max)
        # α = (lam_min + lam_max) / 2
        # δ = α - 0.99 * lam_min
        # pushfirst!(mg_struct.αs, α)
        # pushfirst!(mg_struct.δs, δ)
    end 
    # return lam_mins, lam_maxs, αs, δs
end