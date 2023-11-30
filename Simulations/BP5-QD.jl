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

    # setting up ψδ, odefun for ODEProblem()
    dψV = zeros(2 * odeparam.δNp * odeparam.δNp)
    ψδ = zeros(2 * odeparam.δNp * odeparam.δNp)

    
end