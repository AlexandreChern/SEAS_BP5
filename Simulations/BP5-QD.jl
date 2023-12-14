# loading helper functions for calculations
include("helper.jl")

# loading coefficients for BP5 problem
include("coefficients.jl")

# loading computational domain and linear system
include("domain.jl")

# loading odefun defined for ODEProblem
include("odefun.jl")

function main()
    # loading coefficients for BP5 problem
    # create coefficents if not defined
    if @isdefined BP5_coeff
        println("BP5 coefficients defined, using defined values")
    else
        println("BP5 coefficients not defined, defining it now")
        BP5_coeff = coefficients() # calling default constructor
    end

    # setting b = b0
    b = BP5_coeff.b0
    Vzero = 1e-20

   

    # calling create_ψVδ to create variables used for odefun, create_ψVδ defined in domain.jl
    dψV, ψδ = create_ψVδ()
    # calling create_view to create "views" for dψ, V, ψ, δ, create_view defined in domain.jl
    dψ, V, ψ, δ = create_view(dψV, ψδ)

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
        end
    end

    τ0 = BP5_coeff.σn * BP5_coeff.amax * asinh(BP5_coeff.Vinit / (2 * BP5_coeff.V0) *
            exp.((BP5_coeff.f0 + BP5_coeff.b0 * log.(BP5_coeff.V0 / BP5_coeff.Vinit)) /
            BP5_coeff.amax)) + η * BP5_coeff.Vinit
    # τb = τ0 * V / V_norm

    τz0 = 



    
end