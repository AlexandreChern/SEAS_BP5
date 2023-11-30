using DifferentialEquations
using Printf
using DelimitedFiles
using IterativeSolvers

include("coefficients.jl")
include("helper.jl")
include("domain.jl")

odeparam = (
    reject_step = [false],                          # to reject a step or not
    Vp = BP5_coeff.Vp,                              # plate rate
    M = M,                                          # LHS of the linear system M*u=RHS
    u = zeros(size(RHS)),                           # solution for the linear system 
    u_old = zeros(size(RHS)),                       # solution from the previous step
    Δτ = zeros(2 * (N_x + 1) * (N_y + 1)),          # store the traction computed
    τ = zeros(2 * (N_x + 1) * (N_y + 1)),           # shear stress vector \boldsymbol{τ} = [τ; τ_z]
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
);

struct odeparam_struct
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
end 


# ODE function
function odefun(dψV, ψδ, p, t)

    # # Unpacking named tuple p 
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

    # automatically unpacking named tuples
    @unpack_namedtuple p

    # If reject return
    if reject_step[1]
        return
    end
    # End if reject block

    # Solving linear system using iterative methods
    u_iterative, history = cg(M, RHS, log=true);    # solving with non preconditioned cg
                                                    # this can be replaced with MGCG in future 
    @show history.iters
    u[:] .= u_iterative
    # End of solving 

    # Setting up ratees of change for state and slip
    dψ = @view ψδ[1:δNp]
    δ = @view ψδ[(1:δNp) .+ δNp]

    dψ = @view dψV[1:δNp]
    V = @view dψV[(1:δNp) .+ δNp]

    dψ .= 0
    V .= 0

    # End setting up dψV and ψδ



end
