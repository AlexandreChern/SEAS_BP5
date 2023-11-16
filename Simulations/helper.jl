# Helper functions for BP5-QD and BP5-FD simulations
const year_seconds = 31556926

using DifferentialEquations
using Printf
using Plots

function odefun(dψV, ψδ, p, t)
    
end

struct coefficients
    #       # variable names                # recommended values
    Vp      # plate rate                    # recommended values
    ρ       # density
    cs      # shear wave speed
    a0      # rate-and-state parameter
    amax    # rate-and-state parameter
    b0      # rate-and-state parameter
end