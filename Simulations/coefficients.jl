include("helper.jl")
const year_seconds = 31556926
sim_years = 10

# calling parameterized constructor to set values for BP5
# if BP5_coeff is not defined here, the BP5-QD.jl will 
# call default constructor to construct BP5_coeff using
# values in helper.jl (line 51)
BP5_coeff = coefficients(
    2670,                   # ρ
    3.464,                  # cs
    0.25,                   # ν
    0.004,                  # a0
    0.04,                   # amax
    0.03,                   # b0
    25,                     # σn
    0.14,                   # L or Dc 
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


