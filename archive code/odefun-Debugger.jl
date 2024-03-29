using DifferentialEquations
using Printf
using DelimitedFiles
using IterativeSolvers

include("coefficients.jl")
include("helper.jl")
include("domain.jl")

global const ctr = Ref{Int64}(1)

odeparam = (
    reject_step = [false],                          # to reject a step or not
    Vp = BP5_coeff.Vp,                              # plate rate
    M = M,                                          # LHS of the linear system M*u=RHS
    M_GPU = M_GPU,                                  # GPU array of the LHS system
    u = zeros(size(RHS)),                           # solution for the linear system 
    u_old = zeros(size(RHS)),                       # solution from the previous step
    Δτb = spzeros(2 * (N_x + 1) * (N_y + 1)),       # store the traction computed
    τb = spzeros(2 * (N_x + 1) * (N_y + 1)),        # shear stress vector \boldsymbol{τ} = [τ; τ_z]
    τfb = spzeros(2 * (N_x + 1) * (N_y + 1)),
    counter = [],                                   # counter for slip with Vmax >= threshold
    RHS = RHS,                                      # RHS of the linear system
    μshear = BP5_coeff.cs^2 * BP5_coeff.ρ ,         # constant?
    RSa = BP5_coeff.a0,                             # rate-and-state distance a0
    RSb = BP5_coeff.b0,                             # rate-and-state distance b0
    σn = BP5_coeff.σn,                              # effective normal stress
    η = BP5_coeff.cs * BP5_coeff.ρ / 2,             # constant?
    RSV0 = BP5_coeff.V0,                            # rate-and-state reference slip rate
    τ0 = zeros((N_x + 1) * (N_y + 1)),              # pre-stress                                     # 
    RSL = BP5_coeff.L,                              # rate-and-state critical slip distance L
    RSf0 = BP5_coeff.f0,                            # rate-and-state reference friction coefficient 
    N = N_x,                                        # number of grids in each direction, assuming idential of grid in x,y,z directions
    δNp = N_x + 1,                                  # number of grid points in each direction, assuming idential of grid in x,y,z directions
    Face_operators,                                 # getting face values from 3D SparseArrays
    updators,                                       # updating RHS values using SBP-SAT operators for Dirichlet Operations
    u_filters,                                      # filtering u1, u2, u3 from stacked u
    stride_time = 5
);

struct odeparam_struct
    #=
    # reject_step
    # Vp
    # M
    # u
    # u_old
    # Δτ
    # τ
    # counter
    # ge
    # μsher
    # RSa
    # RSb
    # σn
    # η
    # RSV0
    # τz0
    # RSDc
    # RSf0
    # δNp
    # N
    =#
end 


# ODE function
function odefun(dψV, ψδ, odeparam, t)

    @unpack_namedtuple odeparam;

    # If reject return
    if reject_step[1]
        return
    end
    # End if reject block

    @show t



    dψ, V, ψ, δ = create_view(dψV, ψδ) # creating "views" to get dψ, V, ψ, δ

    dψ .= 0;
    V .= 0;
    ## End setting up dψV and ψδ

    # Updating RHS using δ
    RHS .+= updators[1] * δ[1:2:end]
    RHS .+= updators[2] * δ[2:2:end]

    # Updating RHS using remote loading for face 2 for V2
    RHS .+= updators[3] * (fill(t .* Vp, div(length(δ),2)))

 
    u_iterative, history = cg(M_GPU, CuArray(RHS), log=true);    # solving with non preconditioned cg
    u_iterative, history = cg(M_GPU, CuArray(RHS_new), log=true);    # solving with non preconditioned cg
    # this can be replaced with MGCG in future 
    # @show history.iters
    if typeof(u_iterative) == CuArray{Float64, 1, CUDA.Mem.DeviceBuffer}
    u_iterative = Array(u_iterative)
    end
    u[:] .= u_iterative;

    u1 = u1_filter_matrix * u
    u2 = u2_filter_matrix * u
    u3 = u3_filter_matrix * u

    get_front_face(reshape(u1, Nx, Ny, Nz))
    get_front_face(reshape(u2, Nx, Ny, Nz))
    get_front_face(reshape(u3, Nx, Ny, Nz))


    # End of solving 

    # updating values TODO
    # computate tractions using u: similar to calculation in 
    # assembling 3D matrices
  

    # do non-linear solve for corresponding slip rate TODO
    # given traction on the fault and state variables   
    # check https://github.com/Thrase/Thrase.jl/blob/main/src/odefun.jl
    
    # update rate-and-state region
    # update rate-and-state region using slip rate 
    # and update the region out of rate-and-state 
    # using steady state slip rate

    # len_τ = div(length(Δτb), 2)
    Δτ = @view Δτb[1:2:length(Δτb)] 
    Δτz = @view Δτb[2:2:length(Δτb)]

    τ0 = @view τb[1:2:length(τb)]
    τz0 =  @view τb[2:2:length(τb)]

    # Δτz .= compute_traction_τz() # TODO
    Δτ .= Face_operators[1] * sigma_21 * u_iterative
    Δτz .= Face_operators[1] * sigma_31 * u_iterative
    
    Vn1 = 1e-10 # use very small values
    Vn2 = 1e-10 # use very small values

    Vn1_v = fill(Vn1, fN2 * fN3)
    Vn2_v = fill(Vn2, fN2 * fN3)

    τfb = τb + Δτb 

    τ2 = (τ0 + Δτ)[RS_filter_2D_nzind]
    τ3 = (τz0 + Δτz)[RS_filter_2D_nzind]

    
    (V2_v, V3_v, _, _, iter) = newtbndv_vectorized(rateandstate_vectorized, Vn1_v, Vn2_v, ψ, σn, Vector(τ2), Vector(τ3), RSas, η, RSV0; ftol=1e-12, maxiter=100, atolx=1e-4, rtolx=1e-4)
    @show iter

    # out side of RS, V2 = Vp, V3 = 0
    V[1:2:end] .= Vp
    V[2:2:end] .= 0

    # inside RS, set V2 and V3 with solutions
    V[2 .* RS_filter_2D_nzind .- 1] .= V2_v
    V[2 .* RS_filter_2D_nzind] .= V3_v

    # Update ψ
    # dψ[n] = (RSb * RSV0 / RSDc) * (exp((RSf0 - ψn) / RSb) - abs(Vn) / RSV0) # BP1
    dψ .= (RSb * RSV0 / RSL) .* (exp.((RSf0 .- ψ) ./ RSb) .- sqrt.(V2_v.^2 .+ V3_v.^2) ./ RSV0)
    nothing
end
