using CUDA

mutable struct MG_CUDA
    A_mg
    b_mg
    A_CPU_mg
    λ_mins
    λ_maxs
    odata_mg
    H_mg
    H_inv_mg
    f_mg
    u_mg
    r_mg
    prol_fine_mg
    rest_mg
    prol_mg
    lnx_mg
    lny_mg
    u_exact
    discretization_error

    # struct for mgcg_CUDA
    x_CUDA
    r_CUDA
    r_new_CUDA
    z_CUDA
    z_new_CUDA
    p_CUDA
end