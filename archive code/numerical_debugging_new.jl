# s1 = reshape(HI_3 * e_3 * e_3T * (T_11_3 * u1_filter .+ T_12_3 * u2_filter .+ T_13_3 * u3_filter) * analy_sol, Nx, Ny, Nz)
# s2 = reshape(HI_3 * e_3 * g₁³[:],Nx,Ny,Nz)

# s3 = reshape(HI_4 * e_4 * e_4T * (T_11_4 * u1_filter .+ T_12_4 * u2_filter .+ T_13_4 * u3_filter) * analy_sol, Nx, Ny, Nz)
# s4 = reshape(HI_4 * e_4 * g₁⁴[:], Nx, Ny, Nz)

# s5 = reshape(HI_5 * e_5 * e_5T * (T_11_5 * u1_filter .+ T_12_5 * u2_filter .+ T_13_5 * u3_filter) * analy_sol, Nx, Ny, Nz)
# s6 = reshape(HI_5 * e_5 * g₁⁵[:], Nx, Ny, Nz)

# s7 = reshape(HI_6 * e_6 * e_6T * (T_11_6 * u1_filter .+ T_12_6 * u2_filter .+ T_13_6 * u3_filter) * analy_sol, Nx, Ny, Nz)
# s8 = reshape(HI_6 * e_6 * g₁⁶[:], Nx, Ny, Nz)

# s9 = reshape(HI_3 * e_3 * e_3T * (T_21_3 * u1_filter .+ T_22_3 * u2_filter .+ T_23_3 * u3_filter) * analy_sol, Nx, Ny, Nz)
# s10 = reshape(HI_3 * e_3 * g₂³[:], Nx, Ny, Nz)

# s11 = reshape(HI_4 * e_4 * e_4T * (T_21_4 * u1_filter .+ T_22_4 * u2_filter .+ T_23_4 * u3_filter) * analy_sol, Nx, Ny, Nz)
# s12 = reshape(HI_4 * e_4 * g₂⁴[:], Nx, Ny, Nz)

# s13 = reshape(HI_5 * e_5 * e_5T * (T_21_5 * u1_filter .+ T_22_5 * u2_filter .+ T_23_5 * u3_filter) * analy_sol, Nx, Ny, Nz)
# s14 = reshape(HI_5 * e_5 * g₂⁵[:], Nx, Ny, Nz)

# s15 = reshape(HI_6 * e_5 * e_5T * (T_21_5 * u1_filter .+ T_22_5 * u2_filter .+ T_23_5 * u3_filter) * analy_sol, Nx, Ny, Nz)
# s16 = reshape( HI_6 * e_6 * g₂⁶[:], Nx, Ny, Nz)

# s17 = reshape(HI_3 * e_3 * e_3T * (T_31_3 * u1_filter .+ T_32_3 * u2_filter .+ T_33_3 * u3_filter) * analy_sol, Nx, Ny, Nz)
# s18 = reshape(HI_3 * e_3 * g₃³[:], Nx, Ny, Nz)

# s19 = reshape(HI_4 * e_4 * e_4T * (T_31_4 * u1_filter .+ T_32_4 * u2_filter .+ T_33_4 * u3_filter) * analy_sol, Nx, Ny, Nz)
# s20 = reshape(HI_4 * e_4 * g₃⁴[:], Nx, Ny, Nz)

# s21 = reshape(HI_5 * e_5 * e_5T * (T_31_5 * u1_filter .+ T_32_5 * u2_filter .+ T_33_5 * u3_filter) * analy_sol, Nx, Ny, Nz)
# s22 = reshape(HI_5 * e_5 * g₃⁵[:], Nx, Ny, Nz)

# s23 = reshape(HI_6 * e_6 * e_6T * (T_31_6 * u1_filter .+ T_32_6 * u2_filter .+ T_33_6 * u3_filter) * analy_sol, Nx, Ny, Nz)
# s24 = reshape(HI_6 * e_6 * g₃⁶[:], Nx, Ny, Nz)


# s25 = reshape((T_11_1' .- Z_11_1') * (HI_1 * e_1 * (e_1T)) * u1_filter * analy_sol, Nx, Ny, Nz)
# s26 = reshape((T_11_1' .- Z_11_1') * (HI_1 * e_1 * g₁¹[:]), Nx, Ny, Nz)

# s27 = reshape((T_21_1' .- Z_21_1') * (HI_1 * e_1 * (e_1T)) * u2_filter * analy_sol, Nx, Ny, Nz)
# s28 = reshape((T_21_1' .- Z_21_1') * (HI_1 * e_1 * g₂¹[:]), Nx, Ny, Nz)

# s29 = reshape((T_31_1' .- Z_31_1') * (HI_1 * e_1 * (e_1T)) * u3_filter * analy_sol, Nx, Ny, Nz)
# s30 = reshape((T_31_1' .- Z_31_1') * (HI_1 * e_1 * g₃¹[:]), Nx, Ny, Nz)

# s31 = reshape((T_11_2' .- Z_11_2') * (HI_2 * e_2 * (e_2T)) * u1_filter * analy_sol, Nx, Ny, Nz)
# s32 = reshape((T_11_2' .- Z_11_2') * (HI_2 * e_2 * g₁²[:]), Nx, Ny, Nz)

# s33 = reshape((T_21_2' .- Z_21_2') * (HI_2 * e_2 * (e_2T)) * u2_filter * analy_sol, Nx, Ny, Nz)
# s34 = reshape((T_21_2' .- Z_21_2') * (HI_2 * e_2 * g₂²[:]), Nx, Ny, Nz)

# s35 = reshape((T_31_2' .- Z_31_2') * (HI_2 * e_2 * (e_2T)) * u3_filter * analy_sol, Nx, Ny, Nz)
# s36 = reshape((T_31_2' .- Z_31_2') * (HI_2 * e_2 * g₃²[:]), Nx, Ny, Nz)

# s37 = reshape((T_12_1' .- Z_12_1') * (HI_1 * e_1 * (e_1T)) * u1_filter * analy_sol, Nx, Ny, Nz)
# s38 = reshape((T_12_1' .- Z_12_1') * (HI_1 * e_1 * g₁¹[:]), Nx, Ny, Nz)

# s39 = reshape((T_22_1' .- Z_22_1') * (HI_1 * e_1 * (e_1T)) * u2_filter * analy_sol, Nx, Ny, Nz)
# s40 = reshape((T_22_1' .- Z_22_1') * (HI_1 * e_1 * g₂¹[:]), Nx, Ny, Nz)

# s41 = reshape((T_32_1' .- Z_32_1') * (HI_1 * e_1 * (e_1T)) * u3_filter * analy_sol, Nx, Ny, Nz)
# s42 = reshape((T_32_1' .- Z_32_1') * (HI_1 * e_1 * g₃¹[:]), Nx, Ny, Nz)

# s43 = reshape((T_12_2' .- Z_12_2') * (HI_2 * e_2 * (e_2T)) * u1_filter * analy_sol, Nx, Ny, Nz)
# s44 = reshape((T_12_2' .- Z_12_2') * (HI_2 * e_2 * g₁²[:]), Nx, Ny, Nz)

# s45 = reshape((T_22_2' .- Z_22_2') * (HI_2 * e_2 * (e_2T)) * u2_filter * analy_sol, Nx, Ny, Nz)
# s46 = reshape((T_22_2' .- Z_22_2') * (HI_2 * e_2 * g₂²[:]), Nx, Ny, Nz)

# s47 = reshape((T_32_2' .- Z_32_2') * (HI_2 * e_2 * (e_2T)) * u3_filter * analy_sol, Nx, Ny, Nz)
# s48 = reshape((T_32_2' .- Z_31_2') * (HI_2 * e_2 * g₃²[:]), Nx, Ny, Nz)

# s49 = reshape((T_13_1' .- Z_13_1') * (HI_1 * e_1 * (e_1T)) * u1_filter * analy_sol, Nx, Ny, Nz)
# s50 = reshape((T_13_1' .- Z_13_1') * (HI_1 * e_1 * g₁¹[:]), Nx, Ny, Nz)

# s51 = reshape((T_23_1' .- Z_23_1') * (HI_1 * e_1 * (e_1T)) * u2_filter * analy_sol, Nx, Ny, Nz)
# s52 = reshape((T_23_1' .- Z_23_1') * (HI_1 * e_1 * g₂¹[:]), Nx, Ny, Nz)

# s53 = reshape((T_33_1' .- Z_33_1') * (HI_1 * e_1 * (e_1T)) * u3_filter * analy_sol, Nx, Ny, Nz)
# s54 = reshape((T_33_1' .- Z_33_1') * (HI_1 * e_1 * g₃¹[:]), Nx, Ny, Nz)


# s55 = reshape((T_13_2' .- Z_13_2') * (HI_2 * e_2 * (e_2T)) * u1_filter * analy_sol, Nx, Ny, Nz)
# s56 = reshape((T_13_2' .- Z_13_2') * (HI_2 * e_2 * g₁²[:]), Nx, Ny, Nz)

# s57 = reshape((T_23_2' .- Z_23_2') * (HI_2 * e_2 * (e_2T)) * u2_filter * analy_sol, Nx, Ny, Nz)
# s58 = reshape((T_23_2' .- Z_23_2') * (HI_2 * e_2 * g₂²[:]), Nx, Ny, Nz)

# s59 = reshape((T_33_2' .- Z_33_2') * (HI_2 * e_2 * (e_2T)) * u3_filter * analy_sol, Nx, Ny, Nz)
# s60 = reshape((T_33_2' .- Z_33_2') * (HI_2 * e_2 * g₃²[:]), Nx, Ny, Nz)

function show_last(A;n=9)
    return A[end-n:end, end-n:end]
end


SAT_tilde_1_LHS_new = (
        HI_1 * (T_11_1_new' .- Z_11_1_new') * (e_1 * (e_1T)) * u1_filter
    +   HI_1 * (T_21_1_new' .- Z_21_1_new') * (e_1 * (e_1T)) * u2_filter
    +   HI_1 * (T_31_1_new' .- Z_31_1_new') * (e_1 * (e_1T)) * u3_filter
    +   HI_2 * (T_11_2_new' .- Z_11_2_new') * (e_2 * (e_2T)) * u1_filter
    +   HI_2 * (T_21_2_new' .- Z_21_2_new') * (e_2 * (e_2T)) * u2_filter
    +   HI_2 * (T_31_2_new' .- Z_31_2_new') * (e_2 * (e_2T)) * u3_filter
)

SAT_tilde_1_LHS_new_v2 = HI_tilde * (
        (T_11_1_new' .- Z_11_1_new') * (e_1 * H_1 * (e_1T)) * u1_filter
    +   (T_21_1_new' .- Z_21_1_new') * (e_1 * H_1 *(e_1T)) * u2_filter
    +   (T_31_1_new' .- Z_31_1_new') * (e_1 * H_1 *(e_1T)) * u3_filter
    +   (T_11_2_new' .- Z_11_2_new') * (e_2 * H_2 *(e_2T)) * u1_filter
    +   (T_21_2_new' .- Z_21_2_new') * (e_2 * H_2 *(e_2T)) * u2_filter
    +   (T_31_2_new' .- Z_31_2_new') * (e_2 * H_2 *(e_2T)) * u3_filter
)


# Checking SPD property

# Checking for u_1
# ∂^/∂x^2
E_11 = H_tilde * (K_v + 4/3 * μ_v) * p2_px2_new

SBP_SAT_11 = H_tilde * (
        HI_1 * T_11_1_new' * (e_1 * e_1T)
    +   HI_2 * T_11_2_new' * (e_2 * e_2T)
)

(E_11 + SBP_SAT_1) - (E_11 + SBP_SAT_1)'


# ∂^2/∂y^2
E_12 = H_tilde * (p2_py2_new)


SBP_SAT_12 = H_tilde * beta * ( # beta = -1
        HI_3 * (e_3 * e_3T) * T_11_3_new
    +   HI_4 * (e_4 * e_4T) * T_11_4_new 
)

(E_12 + SBP_SAT_12)
(E_12 + SBP_SAT_12) - (E_12 + SBP_SAT_12)'

# ∂^2/∂z^2
E_13 = H_tilde * (p2_pz2_new)

SBP_SAT_13 = H_tilde * beta * (
        HI_5 * (e_5 * e_5T) * T_11_5_new
    +   HI_6 * (e_6 * e_6T) * T_11_6_new
)

(E_13 + SBP_SAT_13)
(E_13 + SBP_SAT_13) - (E_13 + SBP_SAT_13)'

# ∂^2 / ∂x∂y

E_1xy = H_tilde * p2_pxpy_new;
E_1xz = H_tilde * p2_pxpz_new;
E_1yz = H_tilde * p2_pypz_new;

E_1xy[end-9:end, end-9:end]

SBP_SAT_1xy = H_tilde * beta * (
    HI_1 * p_px_new #* (e_1 * (e_1T)) * u1_filter
+   HI_3 * p_py_new #*  (e_1 * (e_1T)) * u1_filter
# +   HI_2 * p_px_new
# +   HI_2 * p_py_new
# +   HI_3 * p_px_new
# +   HI_3 * p_py_new
);

SBP_SAT_1xy[end-9:end, end-9:end]
# Checking for u2
# ∂^2/∂x^2


p2_pxy_ver1 = p_px_new' * p_py_new

p2_pyx_ver1 = p_py_new' * p_px_new

p2_pxpy_ver1 = H_tilde * (p2_pxy_ver1 + p2_pyx_ver1)


# temp1
temp1 = HI_3 * e_3 * e_3T * (T_11_3_new * u1_filter .+ T_12_3_new * u2_filter .+ T_13_3_new * u3_filter)
temp2 = HI_tilde * e_3 * H_3 * e_3T * (T_11_3_new * u1_filter .+ T_12_3_new * u2_filter .+ T_13_3_new * u3_filter)
@assert temp1 == temp2


temp3 = HI_1 * (T_11_1_new' .- Z_11_1_new') * (e_1 * (e_1T)) * u1_filter
temp4 = HI_tilde * (T_11_1_new' .- Z_11_1_new') * (e_1 * H_1 * (e_1T)) * u1_filter
@assert temp3 == temp4


M_partial_11 = u1_filter' * (
# p2_px2_new
    (
        H_tilde * (K_v + 4/3 * μ_v) * p2_px2_new 
# p_px_hat_new
    +   (T_11_1_new') * (e_1 * H_1 * e_1T)
    +   (T_11_2_new') * (e_2 * H_2 * e_2T)
    ) * u1_filter

# p2_py2_new
+   (
        H_tilde * μ_v * (p2_py2_new ) 
# p_py_hat_new
    +   beta * e_3 * H_3 * e_3T * T_11_3_new
    +   beta * e_4 * H_4 * e_4T * T_11_4_new
    ) * u1_filter

# p2_pz2_new
+   (
        H_tilde * μ_v * (p2_pz2_new)
# p_pz_hat_new
    +   beta * e_5 * H_5 * e_5T * T_11_5_new
    +   beta * e_6 * H_6 * e_6T * T_11_6_new
    ) * u1_filter
)


M_partial_22 = u2_filter' * (
# p2_px2_new
    (
        H_tilde * μ_v * (p2_px2_new)
# p_px_hat_new
    +   T_21_1_new' * (e_1 * H_1 * e_1T)
    +   T_21_2_new' * (e_2 * H_2 * e_2T)
    ) * u2_filter
+
# p2_py2_new
    (
        H_tilde * (K_v + 4/3 * μ_v) * p2_py2_new
# p_py_hat_new
    +   beta * e_3 * H_3 * e_3T * T_22_3_new
    +   beta * e_4 * H_4 * e_4T * T_22_4_new
    )  * u2_filter

# p2_pz2_new
+   (
        H_tilde * μ_v * (p2_pz2_new) 
    +   beta * e_5 * H_5 * e_5T * T_22_5_new
    +   beta * e_6 * H_6 * e_6T * T_22_6_new
    ) * u2_filter
)


M_partial_33 = u3_filter' * (
# p2_px2_new
    (
        H_tilde * μ_v * p2_px2_new 
    +   T_33_1_new' * e_1 * H_1 * e_1T
    +   T_33_2_new' * e_2 * H_2 * e_2T
    ) * u3_filter

+
# p2_py2_new
    (
        H_tilde * μ_v * p2_py2_new
    +   beta * e_3 * H_3 * e_3T * T_33_3_new
    +   beta * e_4 * H_4 * e_4T * T_33_4_new
    ) * u3_filter

# p2_pz2_new
+   (
        H_tilde * (K_v + 4/3 * μ_v) * p2_pz2_new
    +   beta * e_5 * H_5 * e_5T * T_33_5_new
    +   beta * e_6 * H_6 * e_6T * T_33_6_new
    ) * u3_filter
)



temp1 = (u3_filter' * H_tilde * (
        (K_v - 2/3 * μ_v) * (p2_pxpz_new * u1_filter + p2_pypz_new * u2_filter + p2_pz2_new * u3_filter)
    +   2 * μ_v * p2_pz2_new * u3_filter
    )
)

temp2 = (u3_filter' * H_tilde * (
    (K_v + 4/3 * μ_v) * p2_pz2_new
    ) * u3_filter
)

temp3 = (u3_filter' * H_tilde * (
        (K_v - 2/3 * μ_v) * p2_pz2_new
    +   2 * μ_v * p2_pz2_new 
    ) * u3_filter
)


function approx(n)
    if abs(n) ≥ 1e-10
        return 1
    else
        return 0
    end
end


p2_pxpy_2D = kron(D1y, D1x)
u_test = sin.(π.*x .+ π.*y')
plot(x,y,u_test,st=:surface)
plot(x,y,p2_pxpy_2D * u_test[:], st=:surface)
p2_pxpy_2D_alt1 = kron(I_Ny, D1x)' * kron(D1y, I_Ny)
p2_pxpy_2D_alt2 = kron(D1y, I_Nx)'  * kron(I_Ny, D1x)
plot(x,y,p2_pxpy_2D_alt1 * u_test[:], st=:surface)
plot(x,y,p2_pxpy_2D_alt2 * u_test[:], st=:surface)
plot(x,y,(p2_pxpy_2D_alt1 + p2_pxpy_2D_alt2) * u_test[:], st=:surface)


sum(u1' * H_tilde * p2_pxpy_new * u1)
sum(u1' * (H_tilde * p2_pxpy_new * u1))
sum(u1' * (H_tilde * p2_pxpy_alt * u1))
sum(u1' * (H_tilde * p2_pypx_alt * u1))