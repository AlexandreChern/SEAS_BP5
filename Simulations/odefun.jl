using DifferentialEquations
using Printf
using DelimitedFiles
using IterativeSolvers


odeparam = (
    reject_step = [false],
    Vp = Vp,
    M = M,
    u = u,
    u_old = u_old,
    Δτ = Δτ,
    τ = τ,
    counter = [],
    ge = ge,
    μshear = μshear,
    RSa = RSa,
    RSb = RSb,
    σn = σn,
    η = η,
    RSV0 = RSV0,
    τz0 = τz0,
    RSDc = RSDc,
    RSf0 = RSf0,
    σNP = σNp,
    N = N
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
