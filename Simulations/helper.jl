# Helper functions for BP5-QD and BP5-FD simulations
using DifferentialEquations
using Printf
using Plots
using CUDA
using DelimitedFiles
using JSON

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


function rateandstate_vectorized(V2_v, V3_v, ψ, σn, τ2_v, τ3_v, η, RSas, RSV0)
    V  = sqrt.(V2_v.^2 .+ V3_v.^2) 
    dV_dV2_v = 0.5 ./ sqrt.(V2_v.^2 .+ V3_v.^2) .* 2 .* V2_v
    dV_dV3_v = 0.5 ./ sqrt.(V2_v.^2 .+ V3_v.^2) .* 2 .* V3_v

    Y = (1 ./ (2 .* RSV0)) .* exp.(ψ ./ RSas)
    f = RSas .* asinh.(V .* Y)
    df_dV2_v = RSas .* (1 ./ sqrt.(1 .+ (V .* Y).^2)) .* (dV_dV2_v .* Y)
    df_dV3_v = RSas .* (1 ./ sqrt.(1 .+ (V .* Y).^2)) .* (dV_dV3_v .* Y)

    g1 = σn .* f .* V2_v ./ V .+ η .* V2_v .- τ2_v
    g2 = σn .* f .* V3_v ./ V .+ η .* V3_v .- τ3_v

    dA2_dV2_v = (V .- V2_v .* dV_dV2_v) ./ V.^2
    dA2_dV3_v = (-V2_v .* dV_dV3_v) ./ V.^2
    dA3_dV2_v = (-V3_v .* dV_dV2_v) ./ V.^2
    dA3_dV3_v = (V .- V3_v .* dV_dV3_v) ./ V.^2

    dg1_dV2_v = σn .* (df_dV2_v .* V2_v ./ V .+ f .* dA2_dV2_v) .+ η
    dg1_dV3_v = σn .* (df_dV3_v .* V2_v ./ V .+ f .* dA2_dV3_v)

    dg2_dV2_v = σn .* (df_dV2_v .* V3_v ./ V .+ f .* dA3_dV2_v)
    dg2_dV3_v = σn .* (df_dV3_v .* V3_v ./ V .+ f .* dA3_dV3_v) .+ η

    return (g1, g2, dg1_dV2_v, dg1_dV3_v, dg2_dV2_v, dg2_dV3_v)
end

function newtbndv_vectorized(rateandstate_vectorized, V2_v, V3_v, ψ, σn, τ2_v, τ3_v, η, RSas, RSV0;  
                            ftol=1e-12, maxiter=100, α=1.0, atolx = 1e-4, rtolx=1e-4) # change atolx to 1e-8 and rtolx to 1e-8 for better stability
    (f_v, g_v, dfx_v, dfy_v, dgx_v, dgy_v) = rateandstate_vectorized(V2_v, V3_v,  
                        ψ, σn, τ2_v, τ3_v, η, RSas, RSV0)
    # @show V2_v[1], V3_v[1], ψ[1], τ2_v[1], τ3_v[1], RSV0
    for iter = 1:maxiter
        (f_v, g_v, dfx_v, dfy_v, dgx_v, dgy_v) = rateandstate_vectorized(V2_v, V3_v,  
                                        ψ, σn, τ2_v, τ3_v, η, RSas, RSV0)
        inv_J = map_jacobian_inv.(dfx_v, dfy_v, dgx_v, dgy_v)
        # @show dfx_v[1], dfy_v[1], dgx_v[1], dgy_v[1]
        # @show inv_J[1]
        # println()
        dV2_vV3_v =  -inv_J .* map_cat.(f_v, g_v)
        dV2_v = get_first.(dV2_vV3_v)
        dV3_v = get_second.(dV2_vV3_v)
        V2_v = V2_v + α * dV2_v
        V3_v = V3_v + α * dV3_v
        # @show dV2_v[1], dV3_v[1]
        
        # TODO vectorized control flow
        if all(abs.(f_v) .< ftol) && all(abs.(dV2_v) .< atolx .+ rtolx .* (abs.(dV2_v) .+ abs.(V2_v))) && all(abs.(g_v) .< ftol) && all(abs.(dV2_v) .< atolx .+ rtolx .* (abs.(dV3_v) .+ abs.(V3_v)))
            return (V2_v, V3_v, f_v, g_v, iter)
        end
    end
    return (V2_v, V3_v, f_v, g_v, -maxiter)
end


function map_jacobian(x, y, a, b)
    return [x y; a b]
end

function map_jacobian_inv(x, y, a, b)
    return inv([x y; a b])
end

function map_cat(x,y)
    return [x;y]
end

function get_first(a)
    return a[1]
end

function get_second(a)
    return a[2]
end

# Auxiliary functions for debugging
# Plut RS region
function plot_RS(vector)
    spy(reshape(vector, fN2, fN3)')
end
# Plot favorible region
function plot_f(vector)
    spy(reshape(vector, fN2_VW_favorable, fN3_VW)')
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

# Plot slipt against time 
function plot_time_series(path, filenm, region, fieldname, col; headerlines=16, figuresize=(800,600))
    default(fontfamily = "Times", xtickfont = font(20), ytickfont = font(20), legendfont = font(20), guidefont = font(20))

    filename = string(path, filenm)
    A = readdlm(filename, skipstart=headerlines)

    if region == "fault"
        T = A[:, 1] ./ (60*60*24*30) # seconds to years
        slip = A[:, 2]
        V = A[:, 3]
        tau = A[:, 4]
        state = A[:, 5]

        if fieldname == "slip"
            plot(T, slip, color=col, linewidth=2, size=figuresize)
        elseif fieldname == "V"
            plot(T, V, color=col, linewidth=2, size=figuresize)
        elseif fieldname == "10toV"
            plot(T, 10 .^ V, color=col, linewidth=2, size=figuresize)
        elseif fieldname == "tau"
            plot(T, tau, color=col, linewidth=2, size=figuresize)
        else
            plot(T, state, color=col, linewidth=2, size=figuresize)
        end

    elseif region == "surface"
        T = A[:, 1] ./ (60*60*24*30) # seconds to years
        horiz = A[:, 2]
        vert = A[:, 3]
        Vhoriz = A[:, 4]
        Vvert = A[:, 5]

        if fieldname == "horiz"
            plot(T, horiz, color=col)
        elseif fieldname == "vert"
            plot(T, vert, color=col)
        elseif fieldname == "Vhoriz"
            plot(T, Vhoriz, color=col)
        else
            plot(T, Vvert, color=col)
        end
    end

    xlabel!("Time (months)")
    ylabel!(fieldname)
    output_figure_name = split(filename,'.')[1]
    savefig("$output_figure_name.png")
    return
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


function create_text_files(path, station_strings, station_indices, δ, τb, θ, t)

    # path_to_slip = path * "slip.dat"
    # io = open(path_to_slip, "w")
    # write(io, "0.0 0.0 ")
    # for  i in 1:length(flt_loc)
    #     write(io, "$(station_strings[i])")
    # end
    # write(io,"\n")
    # close(io)

    # vv = Array{Float64}(undef, 1, 2+length(flt_loc))
    # vv[1] = t
    # vv[2] = log10(RSVinit)
    # vv[3:end] = δ[flt_loc_indices]

    # io = open(path_to_slip, "a")
    # writedlm(io, vv)
    # close(io)
    Vzero = 1e-20
    for n = 1:length(station_strings)
        XXX = path * "fltst_strk" * station_strings[n] * ".txt"
        ww = Array{Float64}(undef, 1, 8)
        ww[1] = t
        RS_index = RS_filter_2D_nzind[station_indices[n]]
        ww[2] = δ[2 * RS_index-1]
        ww[3] = δ[2 * RS_index]
        # ww[4] = log10(BP5_coeff.Vinit)
        ww[4] = log10(abs(V2_v[station_indices[n]]))
        ww[5] = log10(Vzero)
        ww[6] = τb[2 * RS_index - 1] # need to define this
        ww[7] = τb[2 * RS_index]
        ww[8] = log10(θ[station_indices[n]])  # 
        open(XXX, "w") do io
            write(io, "# This is the file header")
            write(io, "# problem=SEAS Benchmark BP5-QD\n")  # 
            write(io, "# code=Thrase\n")
            write(io, "# modeler=B. A. Erickson\n")
            write(io, "# date=2023/01/09\n")
            write(io, "# element size=1000 m\n")
            write(io, "#location=on fault, 0km along strike, 8km away from the fault, 0km depth")
            write(io, "# minimum_time_step=0.1\n")
            write(io, "# maximum_time_step=3.157e6\n")
            write(io, "# num_time_steps=2400\n")
            write(io, "# Column #1 = Time(s)\n")
            write(io, "# Column #2 = Slip_2(m)\n")
            write(io, "# Column #3 = Slip_3(m)\n")
            write(io, "# Column #4 = Slip_rate_2(log10 m/s)\n")
            write(io, "# Column #5 = Slip_rate_3(log10 m/s)\n")
            write(io, "# Column #6 = Shear_stress_2 (MPa)\n")
            write(io, "# Column #7 = Shear_stress_3 (MPa)\n")
            write(io, "# Column #8 = State (log10 s)\n")
            write(io, "t slip_2 slip_3 slip_rate_2 slip_rate_3 shear_stress_2 shear_stress_3 state\n")
            write(io, "# Here is the time-series data.\n")
            writedlm(io, ww)
        end
    end

end



function write_to_file(path, ψδ, t, i, odeparam, station_strings, station_indices)
    @unpack_namedtuple odeparam;
    Vmax = 0.0
    # if t == (sim_years ./ 31556926)
   
    if isdefined(i, :fsallast)
        dψV = i.fsallast
        dψ, V, ψ, δ = create_view(dψV, ψδ)
        if mod(ctr[], odeparam.stride_time) == 0
            for n = 1:length(station_strings)
                XXX = path * "fltst_strk" * station_strings[n] * ".txt"
                ww = Array{Float64}(undef, 1, 8)
                ww[1] = t
                RS_index = RS_filter_2D_nzind[station_indices[n]] # Check if should be just station_indices[n]
                ww[2] = δ[2 * RS_index-1]
                ww[3] = δ[2 * RS_index]
                V2_real = abs((V[2 * RS_index-1]))
                V3_real = abs((V[2 * RS_index]))
                ww[4] = log10(V2_real)
                ww[5] = log10(V3_real)
                ww[6] = τfb[2 * RS_index - 1] # need to define this
                ww[7] = τfb[2 * RS_index]
                θ = BP5_coeff.L / BP5_coeff.V0 * exp((ψ[station_indices[n]] - BP5_coeff.f0)/BP5_coeff.b0)
                # ww[8] = ψ[station_indices[n]]  # 
                ww[8] = log10(θ)
                open(XXX, "a") do io
                    writedlm(io, ww)
                end
            end
        end
        if mod(ctr[], 1000) == 0
            restart_values = (t = t, 
            ψ = ψ,
            δ = δ,
            V = V)
            open(path * "restart_values.json", "w") do file
                JSON.print(file, restart_values)
            end
        end
        global ctr[] += 1
        @show ctr[]
    end
    Vmax
end


function find_value_in_file(filename, value)
    # Open to file
    open(filename) do file
        for line in each line(file)
            columns = split(line)
            if columns[i] == value
        end
    end
end