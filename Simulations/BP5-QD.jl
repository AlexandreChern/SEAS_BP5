# loading helper functions for calculations
include("helper.jl")


# loading coefficients for BP5 problem
# include("coefficients.jl") # now included in domain.jl

# loading computational domain and linear system
# domain_file = "domain_4km.jl"
domain_file = "domain.jl"
include(domain_file)

# loading odefun defined for ODEProblem
include("odefun.jl")

CUDA.allowscalar(false)

function main()
    # loading coefficients for BP5 problem
    # create coefficents if not defined

    # if @isdefined BP5_coeff
    #     println("BP5 coefficients defined, using defined values")
    # else
    #     println("BP5 coefficients not defined, defining it now")
    #     BP5_coeff = coefficients() # calling default constructor
    # end

    @show BP5_coeff.Δz

    # setting b = b0
    b = BP5_coeff.b0
    Vzero = 1e-20

    @unpack_namedtuple odeparam


    # calling create_ψVδ to create variables used for odefun, create_ψVδ defined in domain.jl
    dψV, ψδ = create_ψVδ()
    # calling create_view to create "views" for dψ, V, ψ, δ, create_view defined in domain.jl
    dψ, V, ψ, δ = create_view(dψV, ψδ)

    # adding temporate θ
    θ = zeros(size(ψ))

    # 
    # RSas = zeros(fN2 * fN3)
    for i in 1:fN2
        for j in 1:fN3
            index = i + (j - 1) * fN2
            x2 = (i - 1) * BP5_coeff.Δz / 1000 - BP5_coeff.lf / 2
            x3 = (j - 1) * BP5_coeff.Δz / 1000
            RSas[index] = a_func(x2, x3, BP5_coeff)
        end
    end
    # RSa_reshaped = reshape(RSa, fN2, fN3)'

    # initializing \boldsymbol{V} over the entire region
    for i in 1:Ny
        for j in 1:Nz
            index = i + (j - 1) * Ny
            V[2*index-1] = BP5_coeff.Vinit
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
            τ0 = BP5_coeff.σn * RSas[index] * asinh((BP5_coeff.Vinit / (2 * BP5_coeff.V0) *
                                                     exp((BP5_coeff.f0 + BP5_coeff.b0 * log(BP5_coeff.V0 / BP5_coeff.Vinit)) /
                                                         RSas[index])) + η * BP5_coeff.Vinit)
            τ[tau_index] = τ0 * BP5_coeff.Vinit / V_norm
            τz[tau_index] = τ0 * Vzero / V_norm

            θ0 = BP5_coeff.L / BP5_coeff.V0 * exp(RSas[index] / BP5_coeff.b0 *
                                                  log(2 * BP5_coeff.V0 / BP5_coeff.Vinit * sinh((τ0 - η * BP5_coeff.Vinit) / (RSas[index] * BP5_coeff.σn)))
                                                  -
                                                  BP5_coeff.f0 / BP5_coeff.b0)
            ψ0 = BP5_coeff.f0 + BP5_coeff.b0 * log(BP5_coeff.V0 * θ0 / BP5_coeff.L)
            ψ[index] = ψ0
            θ[index] = θ0
        end
    end

    # initial velocity for favorable region
    # for index in VW_favorable_filter_RS_nzind
    #     τ_index = RS_filter_2D_nzind[index] 
    #     V[2* τ_index - 1] = 0.03
    #     δ[2* τ_index - 1] = BP5_coeff.L
    #     τ0 = BP5_coeff.σn * RSas[index] * asinh( (0.03 / (2*BP5_coeff.V0) * 
    #                         exp((BP5_coeff.f0 + BP5_coeff.b0 * log(BP5_coeff.V0 / BP5_coeff.Vinit)) / 
    #                         RSas[index]))  + η * 0.03) 
    #     τ[2*τ_index - 1] = τ0                            
    # end
    for i in 1:fN2
        for j in 1:fN3
            index = i + (j - 1) * fN2
            tau_index = RS_filter_2D_nzind[index]
            if index in VW_favorable_filter_RS_nzind
                τ0 = BP5_coeff.σn * RSas[index] * asinh((0.03 / (2 * BP5_coeff.V0) *
                                                         exp((BP5_coeff.f0 + BP5_coeff.b0 * log(BP5_coeff.V0 / BP5_coeff.Vinit)) /
                                                             RSas[index])) + η * 0.03)
                τ[tau_index] = τ0
                V2_v[index] = 0.03
                RSLs[index] = 0.13
                # δ[index] ?
                # δ[2*tau_index - 1] = 0.13? not this meaning, it means L in the rate_and_state_vectorized use a different L for VW region
            end
        end
    end
    # for index in VW_favorable_filter_RS_nzind
    #     V2_v[index] = 0.03
    #     # V3_v[index]
    # end


    # τ0 = BP5_coeff.σn * BP5_coeff.amax * asinh(BP5_coeff.Vinit / (2 * BP5_coeff.V0) *
    #         exp.((BP5_coeff.f0 + BP5_coeff.b0 * log.(BP5_coeff.V0 / BP5_coeff.Vinit)) /
    #         BP5_coeff.amax)) + η * BP5_coeff.Vinit
    # τb = τ0 * V / V_norm

    time_string = Dates.format(now(), "yyyymmddHHMM")
    path_time = "$path/$time_string/"
    try
        mkdir(path_time)
    catch
        # folder already exists
    end

    global ctr[] = 1
    create_text_files(path_time, station_strings, station_indices, δ, τb, θ, 0)
    # Creating output
    callback_func = SavingCallback(
        (ψδ, t, i) -> write_to_file(path_time, ψδ, t, i, odeparam, station_strings, station_indices), SavedValues(Float64, Float64))


    tspan = (0, sim_years * year_seconds)
    prob = ODEProblem(odefun, ψδ, tspan, odeparam)

    function stepcheck(_, odeparam, _)
        if odeparam.reject_step[1]
            odeparam.reject_step[1] = false
            println("reject")
            return true
        end
        return false
    end

    # sol = solve(prob, Tsit5(); dt=0.00001, abstol = 1e-6, reltol = 1e-6, gamma=0.2,save_everystep=true,
    #     callback=callback_func)

    # check atol = 1e-12, rtol = 1e-12
    # check save_everystep = false
    sol = solve(prob, Tsit5(); isoutofdomain=stepcheck, gamma=0.2, dt=0.001, dtmin=1e-8, abstol=1e-12, reltol=1e-12, save_everystep=true,
        callback=callback_func)

end


main()