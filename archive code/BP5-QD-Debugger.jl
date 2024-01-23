# loading helper functions for calculations
include("helper.jl")

# loading coefficients for BP5 problem
include("coefficients.jl")

# loading computational domain and linear system
include("domain.jl")

# loading odefun defined for ODEProblem
include("odefun.jl")


u1_filter_matrix_GPU =CUDA.CUSPARSE.CuSparseMatrixCSR(u1_filter_matrix)
u2_filter_matrix_GPU =CUDA.CUSPARSE.CuSparseMatrixCSR(u2_filter_matrix)
u3_filter_matrix_GPU =CUDA.CUSPARSE.CuSparseMatrixCSR(u3_filter_matrix)


b = BP5_coeff.b0
Vzero = 1e-20

@unpack_namedtuple odeparam;


# calling create_ψVδ to create variables used for odefun, create_ψVδ defined in domain.jl
dψV, ψδ = create_ψVδ()
# calling create_view to create "views" for dψ, V, ψ, δ, create_view defined in domain.jl
dψ, V, ψ, δ = create_view(dψV, ψδ)

# adding temporate θ
θ = zeros(size(ψ))

# 
RSas = zeros(fN2 * fN3)
for i in 1:fN2
    for j in 1:fN3
        index = i + (j - 1) * fN2
        x2 = (i - 1) * BP5_coeff.Δz/1000 - BP5_coeff.lf/2
        x3 = (j - 1) * BP5_coeff.Δz/1000
        RSas[index] = a_func(x2, x3, BP5_coeff)
    end
end
# RSa_reshaped = reshape(RSa, fN2, fN3)'

# initializing \boldsymbol{V} over the entire region
for i in 1:Ny
    for j in 1:Nz
        index = i + (j - 1) * Ny
        V[2*index - 1] = BP5_coeff.Vinit
        V[2*index] = Vzero
    end
end

# initializing \boldsymbol{τ}^0 for the entire domain
V_norm = norm([BP5_coeff.Vinit, Vzero])
τ = @view τb[1:2:length(τb)]
τz = @view τb[2:2:length(τb)]

# only using \tau values for the RS zone

for i in 1:fN2
    for j in 1:fN3
        index = i + (j - 1) * fN2
        tau_index = RS_filter_2D_nzind[index]
        # τ[tau_index] = BP5_coeff.Vinit
        # τz[tau_index] = Vzero
        τ0 = BP5_coeff.σn * RSas[index] * asinh( (BP5_coeff.Vinit / (2*BP5_coeff.V0) * 
                        exp((BP5_coeff.f0 + BP5_coeff.b0 * log(BP5_coeff.V0 / BP5_coeff.Vinit)) / 
                        RSas[index]))  + η * BP5_coeff.Vinit) 
        τ[tau_index] = τ0 * BP5_coeff.Vinit / V_norm
        τz[tau_index] = τ0 * Vzero / V_norm

        θ0 = BP5_coeff.L / BP5_coeff.V0 * exp(RSas[index] / BP5_coeff.b0 * 
        log(2 * BP5_coeff.V0 / BP5_coeff.Vinit * sinh((τ0 - η * BP5_coeff.Vinit ) / (RSas[index] * BP5_coeff.σn)) )
        - BP5_coeff.f0 / BP5_coeff.b0)
        ψ0 = BP5_coeff.f0 + BP5_coeff.b0 * log(BP5_coeff.V0 * θ0 / BP5_coeff.L)
        ψ[index] = ψ0
        θ[index] = θ0
    end
end

try
    mkdir(path)
catch
    # folder already exists
end

# Creating output
callback_func = SavingCallback(
    (ψδ, t, i)->write_to_file(path, ψδ, t, i, odeparam, station_strings, station_indices)
                ,SavedValues(Float64, Float64))
global ctr[] = 1
create_text_files(path, station_strings, station_indices,0)


tspan = (0, sim_years * year_seconds)
prob = ODEProblem(odefun, ψδ, tspan, odeparam)

sol = solve(prob, Tsit5(); dt=0.2, abstol = 1e-5, reltol = 1e-5, gamma=0.2,save_everystep=true,
    callback=callback_func)
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

    #= Unpacking named tuple p 
    # reject_step = p.reject_step
    # Vp = p.Vp
    # Δτ = p.Δτ
    # τ = p.τ
    # τ0 = p.τ0

    # counter = p.counter
    # μshear = p.μshear
    # RSa = p.RSa
    # RSb = p.RSb
    # σn = p.σn
    # η = p.η
    # RSV0 = p.RSV0
    # RSf0 = p.RSf0

    # # domain information
    # N = p.N
    # δNp = p.δNp
    
    # # linear system
    # M = p.M
    # RHS = p.RHS
    # u = p.u
    # u_old = p.u_old

    # # End of unpacking
    =#

    # automatically unpacking named tuples p
    # which is a variable for function odefun

    @unpack_namedtuple odeparam;

    # If reject return
    if reject_step[1]
        return
    end
    # End if reject block

    @show t

    

    ## Setting up ratees of change for state and slip

    # ψ = @view ψδ[1:(fN2 + 1)*(fN3 + 1)]
    # δ = @view ψδ[(fN2 + 1)*(fN3 + 1) + 1:end]

    # dψ = @view dψV[1:(fN2 + 1)*(fN3 + 1)]
    # V = @view dψV[(fN2 + 1)*(fN3 + 1) + 1:end]

    dψ, V, ψ, δ = create_view(dψV, ψδ) # creating "views" to get dψ, V, ψ, δ

    dψ .= 0;
    V .= 0;
    ## End setting up dψV and ψδ

    # Updating RHS using δ
    RHS .+= updators[1] * δ[1:2:end]
    RHS .+= updators[2] * δ[2:2:end]

    # Updating RHS using remote loading for face 2 for V2
    RHS .+= updators[3] * (fill(t .* Vp, div(length(δ),2)))
    # Updating RHS using remote loading for face 2 for V3
    # RSH .+= updators[3] * (fill(0, div(length(δ),2)))

    # End updating RHS using δ

    # Solving linear system using iterative methods
    u_iterative, history = cg(M_GPU, CuArray(RHS), log=true);    # solving with non preconditioned cg
    # this can be replaced with MGCG in future 
    # @show history.iters
    if typeof(u_iterative) == CuArray{Float64, 1, CUDA.Mem.DeviceBuffer}
    u_iterative = Array(u_iterative)
    end
    u[:] .= u_iterative;
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
