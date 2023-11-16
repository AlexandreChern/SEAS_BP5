include("helper.jl")
include("coefficients.jl")


function main()
    if @isdefined BP5_coeff
        println("BP5 coefficients defined, using defined values")
    else
        println("BP5 coefficients not defined, defining it now")
        BP5_coeff = coefficients(
            2670,
            3.464,
            0.25,
            0.004,
            0.04,
            0.03,
            25,
            0.14,
            1E-9,
            1E-9,
            1E-6,
            0.6,
            2,
            2,
            12,
            60,
            40,
            100,
            12,
            1000,
            1800
        )
    end

    

    
end