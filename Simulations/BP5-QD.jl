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
        BP5_coeff = coefficients(
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
            1000,                   # Δz
            1800                    # tf
        )
    end

    # setting b = b0
    b = BP5_coeff.b0

    # setting up ψδ, odefun for ODEProblem()
    dψV = zeros(2 * odeparam.δNp * odeparam.δNp)
    ψδ = zeros(2 * odeparam.δNp * odeparam.δNp)

    
end