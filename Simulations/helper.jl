# Helper functions for BP5-QD and BP5-FD simulations
using DifferentialEquations
using Printf
using Plots

function _unpack(x::NamedTuple)
    kk = keys(x)
    vv = values(x)
    i = 1
    for k in kk
        @eval $k = $vv[$i]
        i += 1
    end
end


macro unpack_namedtuple(arg)
    quote
        _unpack($arg)
    end |> esc
end


struct coefficients
    # Table 1: Parameter values used in this benchmark problem 
    #        # variable names                        # recommended values                
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

# default constructor
coefficients() = coefficients(
    2670,                   # ρ
    3.464,                  # cs
    0.25,                   # ν
    0.004,                  # a0
    0.04,                   # amax
    0.03,                   # b0            value for b in this problem
    25,                     # σn
    0.14,                   # L
    1E-9,                   # Vp
    1E-9,                   # Vinit
    1E-6,                   # V0
    0.6,                    # f0
    2,                      # hs
    2,                      # ht
    12,                     # H
    60,                     # l
    40,                     # Wf
    100,                    # lf
    12,                     # w
    1000,                   # Δz in meter, 
    1800                    # tf
)

# initial state variable over the entire fault
# x2 
function θ0_func(Nx, Ny, BP5_coeff::coefficients)
    return fill(BP5_coeff.L / BP5_coeff.Vinit, Nx * Ny)
end

# for BP5, b is set to be the constant value b0
# using b in the variables of the functions below 
# to write more modular code easier to maintain

function f_func(V, θ, a, b, BP5_coeff::coefficients)
    return a * asinh(
        V / (2 * BP5_coeff.V0) * exp((BP5_coeff.f0 + b * ln(BP5_coeff.V0 / BP5_coeff.L)) / a) # is b the value b0 in coefficients?
    )
end

function a_func(x2, x3, BP5_coeff::coefficients)
    if (BP5_coeff.hs + BP5_coeff.ht ≤ x3 ≤ BP5_coeff.hs + BP5_coeff.ht + BP5_coeff.H) && (abs(x2) ≤ BP5_coeff.l / 2)
        return BP5_coeff.a0
    elseif (0 ≤ x3 ≤ BP5_coeff.hs) || (BP5_coeff.hs + 2 * BP5_coeff.ht + BP5_coeff.H ≤ x3 ≤ BP5_coeff.Wf) || (BP5_coeff.l / 2 + BP5_coeff.ht ≤ abs(x2) ≤ BP5_coeff.lf / 2)
        return BP5_coeff.amax
    else
        r = max(abs(x3 - BP5_coeff.hs - BP5_coeff.ht - BP5_coeff.H / 2) - BP5_coeff.H / 2, abs(x2) - BP5_coeff.l / 2) / BP5_coeff.ht
        return BP5_coeff.a0 + r * (BP5_coeff.amax - BP5_coeff.a0)
    end
end

# auxiliary function to determine which region a belongs to
function a_func_region(x2, x3, BP5_coeff::coefficients)
    if (BP5_coeff.hs + BP5_coeff.ht ≤ x3 ≤ BP5_coeff.hs + BP5_coeff.ht + BP5_coeff.H) && (abs(x2) ≤ BP5_coeff.l / 2)
        return 0
    elseif (0 ≤ x3 ≤ BP5_coeff.hs) || (BP5_coeff.hs + 2 * BP5_coeff.ht + BP5_coeff.H ≤ x3 ≤ BP5_coeff.Wf) || (BP5_coeff.l / 2 + BP5_coeff.ht ≤ abs(x2) ≤ BP5_coeff.lf / 2)
        return 2
    else
        r = max(abs(x3 - BP5_coeff.hs - BP5_coeff.ht - BP5_coeff.H / 2) - BP5_coeff.H / 2, abs(x2) - BP5_coeff.l / 2) / BP5_coeff.ht
        return 1
    end
end

# For BP5-QD, the scalar pre-stress τ⁰ is chosen as the steady-state stress
function τ0_QD_func(a, b, η, BP5_coeff::coefficients)
    return BP5_coeff.σn * a * asinh(BP5_coeff.Vinit / (2 * BP5_coeff.V0) * exp((BP5_coeff.f0 + b) / a)) + η * BP5_coeff.Vinit
end

# For BP5-QD, the scalar pre-stress τ⁰ is chosen as the steady-state stress
function τ0_FD_func(a, b, BP5_coeff::coefficients)
    return BP5_coeff.σn * asinh(BP5_coeff.Vinit / (2 * BP5_coeff.V0) * exp((BP5_coeff.f0 + b) / a))
end

# a higher pre-stress along the x2-direction
function τi0_FD_func(a, b, δ, τ, BP5_coeff::coefficients)
    return BP5_coeff.σn * asinh(BP5_coeff.Vinit / (2 * BP5_coeff.V0) * exp((BP5_coeff.f0 + b) / a)) + δ * τ
end


# quasi-static process zone
function Λ0_func(C, b, μ, BP5_coeff::coefficients)
    return C * μ * BP5_coeff.L / (b * BP5_coeff.σn)
end

# nucleation zone
function h_func(a, b, μ, BP5_coeff::coefficients)
    return π / 2 * (μ * b * BP5_coeff.L) / ((b - a)^2 * BP5_coeff.σn)^2
end

# fault strength
function F_func(f, Vbold, BP5_coeff::coefficients)
    return BP5_coeff.σn * f * Vbold / norm(V)
end


# boundary functions
# Dirichlet boundary conditions
function bc_Dirichlet(face, Vp, δ, x, t)
    if face == 1 # check this
        return (δ ./ 2)
    elseif face == 2 # check the 
        return fill(t * Vp / 2, size(x))
    end
end

# Neumann boundary conditions
function bc_Neumann(face, Vp, δ, x, t)
    return zeros(size(x))
end


# update boundary conditions
function boundary_update!(RHS, bc_Dirichlet, bc_Neumann)
    # TODO
end

# test functions 
# let    
#     a_func(20,20,BP5_coeff)
# end

# # rate and state function scalar values
# function rateandstate(V, psi, σn, ϕ, η, a, V0)
#     Y = (1 ./ (2 .* V0)) .* exp.(psi ./ a)
#     f = a .* asinh.(V .* Y)
#     dfdV = a .* (1 ./ sqrt.(1 + (V .* Y) .^ 2)) .* Y

#     g = σn .* f + η .* V - ϕ
#     dgdV = σn .* dfdV + η
#     (g, dgdV)
# end

# # newtom method for solving V
# function newtbndv(func, xL, xR, x; ftol=1e-6, maxiter=500, minchange=0,
#     atolx=1e-4, rtolx=1e-4)
#     (fL, _) = func(xL)
#     (fR, _) = func(xR)
#     if fL .* fR > 0
#         return (typeof(x)(NaN), typeof(x)(NaN), -maxiter)
#     end

#     (f, df) = func(x)
#     dxlr = xR - xL

#     for iter = 1:maxiter
#         dx = -f / df
#         x = x + dx

#         if x < xL || x > xR || abs(dx) / dxlr < minchange
#             x = (xR + xL) / 2
#             dx = (xR - xL) / 2
#         end

#         (f, df) = func(x)

#         if f * fL > 0
#             (fL, xL) = (f, x)
#         else
#             (fR, xR) = (f, x)
#         end
#         dxlr = xR - xL

#         if abs(f) < ftol && abs(dx) < atolx + rtolx * (abs(dx) + abs(x))
#             return (x, f, iter)
#         end
#     end
#     return (x, f, -maxiter)
# end


function rateandstate(V2, V3, psi, σn, τ2, τ3, η, a, V0)

    V = sqrt(V2^2 + V3^2)
    dV_dV2 = 0.5*(V2^2 + V3^2)^(-0.5) .* 2 .* V2
    dV_dV3 = 0.5*(V2^2 + V3^2)^(-0.5) .* 2 .* V3
    
    Y = (1 ./ (2 .* V0)) .* exp.(psi ./ a)
    f = a .* asinh.(V .* Y)  # compute friction coefficient 
    df_dV2  = a .* (1 ./ sqrt.(1 + (V .* Y).^2)) .* (dV_dV2 .* Y)  # derivative wrt V_2
    df_dV3  = a .* (1 ./ sqrt.(1 + (V .* Y).^2)) .* (dV_dV3 .* Y)  # derivative wrt V_2

    g1 = σn .* f .* V2 / V   + η .* V2 - τ2
    g2 = σn .* f .* V3 / V   + η .* V3 - τ3
   

    A2 = V2/V
    A3 = V3/V

    dA2_dV2 = (V - V2*dV_dV2)/V^2
    dA2_dV3 = (-V2*dV_dV3)/V^2
    dA3_dV2 = (-V3*dV_dV2)/V^2
    dA3_dV3 = (V - V3*dV_dV3)/V^2

    dg1_dV2 = σn .* (df_dV2 .* V2 / V + f .* dA2_dV2) + η
    dg1_dV3 = σn .* (df_dV3 .* V2 / V + f .* dA2_dV3)

    dg2_dV2 = σn .* (df_dV2 .* V3 / V + f .* dA3_dV2) 
    dg2_dV3 = σn .* (df_dV3 .* V3 / V + f .* dA3_dV3) + η 

    return (g1, g2, dg1_dV2, dg1_dV3, dg2_dV2, dg2_dV3)
end
  

function newtbndv(func, x, y; ftol = 1e-12, maxiter = 500, 
                    atolx = 1e-4, rtolx = 1e-4)

    (f, g, dfx, dfy, dgx, dgy) = func(x, y)
    for iter = 1:maxiter

        z = [x; y] 
        (f, g, dfx, dfy, dgx, dgy) = func(x, y)
      
        J = [dfx dfy; dgx dgy] 
        dx, dy = -J\[f; g]
      
        x = x + dx
        y = y + dy
  
       if abs(f) < ftol && abs(dx) < atolx + rtolx * (abs(dx) + abs(x)) && abs(g) < ftol && abs(dy) < atolx + rtolx * (abs(dy) + abs(y))
            
            return (x, y, f, g, iter)
       end
    end
    return (x, y, f, g, -maxiter)
end


# Plot the slip in 2D from BP1 problem
function plot_slip(S, δNp, yf, stride_time)

    m = length(yf)
    no_time_steps = size(S.t)
    slip_final = S.u[end][end]

    for i = 1:stride_time:no_time_steps[1]

        slip_t = S.u[i][δNp+1:end] # slip at time t
        #pyplot()
        display(plot(slip_t, -yf, xtickfont=font(18),
            ytickfont=font(18),
            guidefont=font(18),
            legendfont=font(18), ylabel="Depth (km)", xlabel="Slip (m)", xlims=(0, slip_final)))
        sleep(0.1)
    end

    #nothing
end



########################## Test functions, subjected to changes ################################


function test_unpack(p)
    @unpack_namedtuple p
    @show Vp
    @show RHS
    # for i in keys(p)
    #     @show @eval $(i)
    # end
end