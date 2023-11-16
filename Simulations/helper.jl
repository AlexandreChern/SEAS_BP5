# Helper functions for BP5-QD and BP5-FD simulations
const year_seconds = 31556926

using DifferentialEquations
using Printf
using Plots

function odefun(dψV, ψδ, p, t)
    
end

struct coefficients
    # Table 1: Parameter values used in this benchmark problem 
            # variable names                        # recommended values                
    ρ       # density                               2670 kg/m³
    cs      # shear wave speed                      3.464 km/s
    ν       # Poisson's ratio                       0.25
    a0      # rate-and-state parameter              0.004
    amax    # rate-and-state parameter              0.04
    b0      # rate-and-state parameter              0.03
    σn      # effective normal stress               25 MPa
    L       # critical slip distance                0.14 m/0.13 m 
    Vp      # plate rate                            10⁻⁹ m/s
    Vinit   # initial slip rate                     10⁻⁹ m/s
    V0      # reference slip rate                   10⁻⁶ m/s 
    f0      # reference friction coefficient        0.6
    hs      # width of shallow VS zone              2 km
    ht      # width of VW-VS transition zone        2 km
    H       # width of uniform VW region            12 km
    l       # length of uniform VW region           60 km
    Wf      # width of rate-and-state fault         40 km
    lf      # length of rate-and-state fault        100 km
    w       # width of favorible nucleation zone    12 km
    Δz      # suggested cell size                   1000m
    tf      # final simulation time                 1800 years
end



function f_func(V, θ, a, BP5_coeff::coefficients)
    return a * asinh(
        V / (2 * BP5_coeff.V0) * exp((BP5_coeff.f0 + BP5_coeff.b0 * ln(BP5_coeff.V0 / BP5_coeff.L) ) / a)
    )
end

function a_func(x2, x3, BP5_coeff::coefficients)
    if (BP5_coeff.hs + BP5_coeff.ht ≤ x3 ≤ BP5_coeff.hs + BP5_coeff.ht + BP5_coeff.H) && (abs(x2) ≤ BP5_coeff.l / 2)
        return BP5_coeff.a0
    elseif (0 ≤ x3 ≤ BP5_coeff.hs) || (BP5_coeff.hs + 2 * BP5_coeff.ht + BP5_coeff.H ≤ x3 ≤ BP5_coeff.Wf)
        return BP5_coeff.amax
    else
        r = max(abs(x3 - BP5_coeff.hs - BP5_coeff.ht - BP5_coeff.H/2) - BP5_coeff.H/2, abs(x2) - BP5_coeff.l/2) / BP5_coeff.ht
        return BP5_coeff.a0 + r * (BP5_coeff.amax - BP5_coeff.a0)
    end
end


a_func(20,20,BP5_coeff)