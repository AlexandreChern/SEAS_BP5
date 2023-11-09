using CUDA
using Arpack
include("interpolations.jl")

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
    lnz_mg
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

mg_struct_CUDA = MG_CUDA([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[])

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
    mg_struct_CUDA.lnz_mg = []
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
    lnz_mg = mg_struct_CUDA.lnz_mg
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
                push!(H_mg, CUDA.CUSPARSE.CuSparseMatrixCSR(kron(H_tilde, sparse(I,3,3)))) # kron(H_tilde, sparse(I,3,3))
                push!(H_inv_mg, CUDA.CUSPARSE.CuSparseMatrixCSR(kron(HI_tilde, sparse(I,3,3))))
                push!(f_mg, CuArray(zeros(size(b))))
                push!(r_mg, CuArray(zeros(size(b))))
                push!(u_mg, CuArray(zeros(size(b))))
                push!(u_exact, analy_sol)
                push!(prol_fine_mg, CuArray(zeros(size(b))))
                push!(rest_mg, CUDA.CUSPARSE.CuSparseMatrixCSR(kron(restriction_matrix_v0(nx,ny,nz,div(nx,2),div(ny,2),div(nz,2)),sparse(I,3,3))))
                push!(prol_mg, CUDA.CUSPARSE.CuSparseMatrixCSR(kron(prolongation_matrix_v0(nx,ny,nz,div(nx,2),div(ny,2),div(nz,2)),sparse(I,3,3))))
            else
                A, b, H_tilde, HI_tilde, analy_sol = Assembling_3D_matrices(nx,ny,nz)
                push!(A_CPU_mg, A)
                push!(A_mg, CUDA.CUSPARSE.CuSparseMatrixCSR(A))
                push!(b_mg, CuArray(b))
                push!(H_mg, CUDA.CUSPARSE.CuSparseMatrixCSR(kron(H_tilde, sparse(I,3,3))))
                push!(H_inv_mg, CUDA.CUSPARSE.CuSparseMatrixCSR(kron(HI_tilde, sparse(I,3,3))))
                push!(f_mg, CuArray(zeros(size(b))))
                push!(r_mg, CuArray(zeros(size(b))))
                push!(u_mg, CuArray(zeros(size(b))))
                push!(u_exact, analy_sol)
                push!(prol_fine_mg, CuArray(zeros(size(b))))
                push!(rest_mg, CUDA.CUSPARSE.CuSparseMatrixCSR(kron(restriction_matrix_v0(nx,ny,nz,div(nx,2),div(ny,2),div(nz,2)),sparse(I,3,3))))
                push!(prol_mg, CUDA.CUSPARSE.CuSparseMatrixCSR(kron(prolongation_matrix_v0(nx,ny,nz,div(nx,2),div(ny,2),div(nz,2)),sparse(I,3,3))))
            end
            push!(lnx_mg, nx)
            push!(lny_mg, ny)
            push!(lnz_mg, nz)
            nx, ny, nz = div(nx,2), div(ny,2), div(nz,2)
            hx, hy = 2*hx, 2*hy, 2*hz
        end 
    end
    get_lams(mg_struct_CUDA)
    # # struct for mgcg_CUDA
    # mg_struct_CUDA.x_CUDA
    # mg_struct_CUDA.r_CUDA
    # mg_struct_CUDA.r_new_CUDA
    # mg_struct_CUDA.z_CUDA
    # mg_struct_CUDA.z_new_CUDA
    # mg_struct_CUDA.p_CUDA
    push!(mg_struct_CUDA.x_CUDA, CuArray(zeros(size(mg_struct_CUDA.b_mg[1]))))
    push!(mg_struct_CUDA.r_CUDA, CuArray(zeros(size(mg_struct_CUDA.b_mg[1]))))
    push!(mg_struct_CUDA.r_new_CUDA, CuArray(zeros(size(mg_struct_CUDA.b_mg[1]))))
    push!(mg_struct_CUDA.z_CUDA, CuArray(zeros(size(mg_struct_CUDA.b_mg[1]))))
    push!(mg_struct_CUDA.z_new_CUDA, CuArray(zeros(size(mg_struct_CUDA.b_mg[1]))))
    push!(mg_struct_CUDA.p_CUDA, CuArray(zeros(size(mg_struct_CUDA.b_mg[1]))))
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
    # TODO rewrite the function to get correct interpolation for eigenvalues
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
        if size(reverse_Amg[k])[1] <= 2187 # nx <= 8
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
                # lam_max = mg_struct_CUDA.λ_maxs[1] + (mg_struct_CUDA.λ_maxs[1] - mg_struct_CUDA.λ_maxs[2]) * 0.6
                lam_max = mg_struct_CUDA.λ_maxs[1] / 2
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


function mg_solver_CUDA(mg_struct_CUDA, f_in; 
                        nx=16, ny=16, nz=16, 
                        n_levels=3, v1=5, v2=5, v3=5,
                        tolerance=1e-10,
                        max_mg_iterations=1, 
                        use_direct_sol=false,
                        dynamic_richardson_ω=false)
    if isempty(mg_struct_CUDA.A_mg)
        initialize_mg_struct_CUDA(mg_struct_CUDA,nx,ny,nz,n_levels)
    end
    clear_urf_CUDA(mg_struct_CUDA)
    println("Starting Multigrid V-cycle")

    mg_struct_CUDA.f_mg[1][:] .= copy(f_in)[:]
    mg_struct_CUDA.r_mg[1][:] .= mg_struct_CUDA.f_mg[1][:] .- mg_struct_CUDA.A_mg[1] * mg_struct_CUDA.u_mg[1]
    @show norm(mg_struct_CUDA.r_mg[1])
    mg_iter_count = 0

    if nx < (2^n_levels)
        println("Number of levels exceeds the possible number.")
        return 0
    end

    for iteration_count in 1:max_mg_iterations
        mg_iter_count += 1
        @show mg_iter_count

        ω_richardson = 2 / (mg_struct_CUDA.λ_mins[1] + mg_struct_CUDA.λ_maxs[1])
        for i in 1:v1
            mg_struct_CUDA.u_mg[1][:] .+= ω_richardson * (mg_struct_CUDA.f_mg[1][:] .- mg_struct_CUDA.A_mg[1] * mg_struct_CUDA.u_mg[1][:])
            mg_struct_CUDA.r_mg[1][:] .= mg_struct_CUDA.f_mg[1][:] .- mg_struct_CUDA.A_mg[1] * mg_struct_CUDA.u_mg[1]
            @show norm(mg_struct_CUDA.r_mg[1][:])
        end

        mg_struct_CUDA.r_mg[1][:] .= mg_struct_CUDA.f_mg[1][:] .- mg_struct_CUDA.A_mg[1] * mg_struct_CUDA.u_mg[1]

        for k in 2:n_levels
            ω_richardson = 2 / (mg_struct_CUDA.λ_mins[k] + mg_struct_CUDA.λ_maxs[k])
            if k == 2
                mg_struct_CUDA.r_mg[k-1] = mg_struct_CUDA.r_mg[1]
            else
                mg_struct_CUDA.r_mg[k-1][:] .= mg_struct_CUDA.f_mg[k-1][:] .- mg_struct_CUDA.A_mg[k-1] * mg_struct_CUDA.u_mg[k-1]
            end

            mg_struct_CUDA.f_mg[k] .= mg_struct_CUDA.H_mg[k] * mg_struct_CUDA.rest_mg[k-1] * mg_struct_CUDA.H_inv_mg[k-1] * mg_struct_CUDA.r_mg[k-1]

            if k < n_levels
                println("pre-smoothing")
                for i in 1:v1
                    mg_struct_CUDA.u_mg[k][:] .+= ω_richardson * (mg_struct_CUDA.f_mg[k][:] .- mg_struct_CUDA.A_mg[k] * mg_struct_CUDA.u_mg[k][:])
                    mg_struct_CUDA.r_mg[k][:] .= mg_struct_CUDA.f_mg[k][:] .- mg_struct_CUDA.A_mg[k] * mg_struct_CUDA.u_mg[k]
                    @show norm(mg_struct_CUDA.r_mg[k][:])     
                end
                mg_struct_CUDA.r_mg[k][:] .= mg_struct_CUDA.f_mg[k][:] .- mg_struct_CUDA.A_mg[k] * mg_struct_CUDA.u_mg[k]
            elseif k == n_levels
                println("coarsest grid smoothing")
                for i in 1:v2
                    mg_struct_CUDA.u_mg[k][:] .+= ω_richardson * (mg_struct_CUDA.f_mg[k][:] .- mg_struct_CUDA.A_mg[k] * mg_struct_CUDA.u_mg[k][:])
                    mg_struct_CUDA.r_mg[k][:] .= mg_struct_CUDA.f_mg[k][:] .- mg_struct_CUDA.A_mg[k] * mg_struct_CUDA.u_mg[k]
                    @show norm(mg_struct_CUDA.r_mg[k][:])     
                end
            end
        end

        println("Post smoothing")

        for k = n_levels:-1:2
            ω_richardson = 2 / (mg_struct_CUDA.λ_mins[k] + mg_struct_CUDA.λ_maxs[k])
            mg_struct_CUDA.prol_fine_mg[k-1] = mg_struct_CUDA.prol_mg[k-1] * mg_struct_CUDA.u_mg[k]
            mg_struct_CUDA.u_mg[k-1] .+= mg_struct_CUDA.prol_fine_mg[k-1]
            ω_richardson = 2 / (mg_struct_CUDA.λ_mins[k-1] + mg_struct_CUDA.λ_maxs[k-1])
            for i in 1:v3
                mg_struct_CUDA.u_mg[k][:] .+= ω_richardson * (mg_struct_CUDA.f_mg[k][:] .- mg_struct_CUDA.A_mg[k] * mg_struct_CUDA.u_mg[k][:])    
                mg_struct_CUDA.r_mg[k][:] .= mg_struct_CUDA.f_mg[k][:] .- mg_struct_CUDA.A_mg[k] * mg_struct_CUDA.u_mg[k]
                @show norm(mg_struct_CUDA.r_mg[k][:])        
            end
        end
        mg_struct_CUDA.r_mg[1][:] .= mg_struct_CUDA.f_mg[1][:] .- mg_struct_CUDA.A_mg[1] * mg_struct_CUDA.u_mg[1][:]

        println("Finish V-cycle")
        @show norm(mg_struct_CUDA.r_mg[1][:])
    end
end


function mgcg_CUDA(mg_struct_CUDA;nx=64,ny=64,n_levels=3,v1=5,v2=5,v3=5, ω=1.0, ω_richardson=2/1000, max_cg_iter=10, max_mg_iterations=1,iter_algo_num=3, precond=true,dynamic_richardson_ω=false)
    if nx != mg_struct_CUDA.lnx_mg[1]
        clear_mg_struct_CUDA(mg_struct_CUDA)
        initialize_mg_struct_CUDA(mg_struct_CUDA, nx, ny, nz, n_level)
    end
   
    mg_struct_CUDA.r_CUDA[1][:] .= mg_struct_CUDA.b_mg[1][:] .- mg_struct_CUDA.A_mg[1] * mg_struct_CUDA.x_CUDA[1][:]

    init_rms = norm(mg_struct_CUDA.r_CUDA[1])
    # z_CUDA = CuArray(zeros(nx+1,ny+1))

    if precond == true
        mg_solver_CUDA(mg_struct_CUDA, mg_struct_CUDA.r_CUDA[1], n_levels = n_levels)
        mg_struct_CUDA.z_CUDA[1] .= mg_struct_CUDA.u_mg[1]
    else
        mg_struct_CUDA.z_CUDA[1][:] .= mg_struct_CUDA.r_CUDA[1][:]
    end

    # p_CUDA = CuArray(zeros(size(r_CUDA)))
    # p_CUDA .= z_CUDA
    mg_struct_CUDA.p_CUDA[1] .= mg_struct_CUDA.z_CUDA[1]

    counter = 0
    for k in 1:max_cg_iter
        mfA_CUDA(mg_struct_CUDA.p_CUDA[1],mg_struct_CUDA,1)

        # α = dot(r_CUDA[:],z_CUDA[:]) / (dot(p_CUDA[:],A_CUDA * p_CUDA[:]))
        α = dot(mg_struct_CUDA.r_CUDA[1][:], mg_struct_CUDA.z_CUDA[1][:]) / (dot(mg_struct_CUDA.p_CUDA[1][:],mg_struct_CUDA.odata_mg[1]))

        mg_struct_CUDA.x_CUDA[1] .+= α .* mg_struct_CUDA.p_CUDA[1]


        # r_new_CUDA = r_CUDA[:] .- α * A_CUDA * p_CUDA[:]
        mg_struct_CUDA.r_new_CUDA[1][:] = mg_struct_CUDA.r_CUDA[1][:] .- α * mg_struct_CUDA.odata_mg[1][:]

        norm_v_initial_norm = norm(mg_struct_CUDA.r_new_CUDA[1]) / init_rms
        # @show k, norm(x_CUDA[:] - mg_struct_CUDA.u_exact[1][:])
        
        # L2_error = sqrt((x_CUDA[:] - mg_struct_CUDA.u_exact[1][:])' * CUDA.CUSPARSE.CuSparseMatrixCSR(mg_struct_CUDA.H_mg[1]) * (x_CUDA[:] - mg_struct_CUDA.u_exact[1][:]) )

        mfA_H((mg_struct_CUDA.x_CUDA[1][:] - mg_struct_CUDA.u_exact[1][:]),mg_struct_CUDA,1)
        L2_error = sqrt(dot((mg_struct_CUDA.x_CUDA[1][:] - mg_struct_CUDA.u_exact[1][:])', mg_struct_CUDA.odata_mg[1]))

        # @show k, norm_v_initial_norm, L2_error

        # println("")

        if norm(mg_struct_CUDA.r_new_CUDA[1]) < 1e-8 * init_rms
            break
        end

        if precond == true
            mg_solver_CUDA(mg_struct_CUDA, mg_struct_CUDA.r_new_CUDA[1], n_level = n_level, v1=v1,v2=v2,v3=v3, max_mg_iterations=max_mg_iterations, nx = nx, ny = ny, iter_algo_num=iter_algo_num, dynamic_richardson_ω=dynamic_richardson_ω)
            mg_struct_CUDA.z_new_CUDA[1] .= mg_struct_CUDA.u_mg[1]
        else
            mg_struct_CUDA.z_new_CUDA[1] .= copy(mg_struct_CUDA.r_new_CUDA[1])
        end

        β = dot(mg_struct_CUDA.r_new_CUDA[1][:], mg_struct_CUDA.z_new_CUDA[1][:]) / (dot(mg_struct_CUDA.r_CUDA[1][:],mg_struct_CUDA.z_CUDA[1][:]))

        mg_struct_CUDA.p_CUDA[1][:] .= mg_struct_CUDA.z_new_CUDA[1][:] .+ β * mg_struct_CUDA.p_CUDA[1][:]
        mg_struct_CUDA.z_CUDA[1][:] .= mg_struct_CUDA.z_new_CUDA[1][:]
        mg_struct_CUDA.r_CUDA[1][:] .= mg_struct_CUDA.r_new_CUDA[1][:]
        counter += 1

        # η_alg = dot((x_CUDA[:]), mg_struct_CUDA.A_mg[1] * (x_CUDA[:]))
        # η_tot = dot((mg_struct_CUDA.u_exact[1])[:], mg_struct_CUDA.A_mg[1] * mg_struct_CUDA.u_exact[1][:])
        # @show k, η_alg, η_tot, η_alg / η_tot
    end
    return mg_struct_CUDA.x_CUDA[1], counter     
end