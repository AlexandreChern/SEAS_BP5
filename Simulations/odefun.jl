using DifferentialEquations
using Printf
using DelimitedFiles
using IterativeSolvers

include("coefficients.jl")
include("helper.jl")
include("domain.jl")

odeparam = (
    reject_step = [false],
    Vp = BP5_coeff.Vp,
    M = M,
    u = zeros(size(RHS)),
    u_old = zeros(size(RHS)),
    Δτ = Δτ,
    τ = τ,
    counter = [],
    ge = ge,
    μshear = BP5_coeff.cs^2 * BP5_coeff.ρ ,
    RSa = BP5_coeff.a0,
    RSb = BP5_coeff.b0,
    σn = BP5_coeff.σn,
    η = BP5_coeff.cs * BP5_coeff.ρ / 2,
    RSV0 = BP5_coeff.V0,
    τz0 = τz0,
    RSDc = RSDc,
    RSf0 = RSf0,
    N = N_x,
    δNp = N_x + 1,
)

struct odeparam_struct
    reject_step
    Vp
    M
    u
    u_old
    Δτ
    τ
    counter
    ge
    μsher
    RSa
    RSb
    σn
    η
    RSV0
    τz0
    RSDc
    RSf0
    δNp
    N
end 


# ODE function
function odefun(dψV, ψδ, p, t)
    # TO DO
end
