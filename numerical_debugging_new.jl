s1 = reshape(HI_3 * e_3 * e_3T * (T_11_3 * u1_filter .+ T_12_3 * u2_filter .+ T_13_3 * u3_filter) * analy_sol, Nx, Ny, Nz)
s2 = reshape(HI_3 * e_3 * g₁³[:],Nx,Ny,Nz)

s3 = reshape(HI_4 * e_4 * e_4T * (T_11_4 * u1_filter .+ T_12_4 * u2_filter .+ T_13_4 * u3_filter) * analy_sol, Nx, Ny, Nz)
s4 = reshape(HI_4 * e_4 * g₁⁴[:], Nx, Ny, Nz)

s5 = reshape(HI_5 * e_5 * e_5T * (T_11_5 * u1_filter .+ T_12_5 * u2_filter .+ T_13_5 * u3_filter) * analy_sol, Nx, Ny, Nz)
s6 = reshape(HI_5 * e_5 * g₁⁵[:], Nx, Ny, Nz)

s7 = reshape(HI_6 * e_6 * e_6T * (T_11_6 * u1_filter .+ T_12_6 * u2_filter .+ T_13_6 * u3_filter) * analy_sol, Nx, Ny, Nz)
s8 = reshape(HI_6 * e_6 * g₁⁶[:], Nx, Ny, Nz)

s9 = reshape(HI_3 * e_3 * e_3T * (T_21_3 * u1_filter .+ T_22_3 * u2_filter .+ T_23_3 * u3_filter) * analy_sol, Nx, Ny, Nz)
s10 = reshape(HI_3 * e_3 * g₂³[:], Nx, Ny, Nz)

s11 = reshape(HI_4 * e_4 * e_4T * (T_21_4 * u1_filter .+ T_22_4 * u2_filter .+ T_23_4 * u3_filter) * analy_sol, Nx, Ny, Nz)
s12 = reshape(HI_4 * e_4 * g₂⁴[:], Nx, Ny, Nz)

s13 = reshape(HI_5 * e_5 * e_5T * (T_21_5 * u1_filter .+ T_22_5 * u2_filter .+ T_23_5 * u3_filter) * analy_sol, Nx, Ny, Nz)
s14 = reshape(HI_5 * e_5 * g₂⁵[:], Nx, Ny, Nz)

s15 = reshape(HI_6 * e_5 * e_5T * (T_21_5 * u1_filter .+ T_22_5 * u2_filter .+ T_23_5 * u3_filter) * analy_sol, Nx, Ny, Nz)
s16 = reshape( HI_6 * e_6 * g₂⁶[:], Nx, Ny, Nz)

s17 = reshape(HI_3 * e_3 * e_3T * (T_31_3 * u1_filter .+ T_32_3 * u2_filter .+ T_33_3 * u3_filter) * analy_sol, Nx, Ny, Nz)
s18 = reshape(HI_3 * e_3 * g₃³[:], Nx, Ny, Nz)

s19 = reshape(HI_4 * e_4 * e_4T * (T_31_4 * u1_filter .+ T_32_4 * u2_filter .+ T_33_4 * u3_filter) * analy_sol, Nx, Ny, Nz)
s20 = reshape(HI_4 * e_4 * g₃⁴[:], Nx, Ny, Nz)

s21 = reshape(HI_5 * e_5 * e_5T * (T_31_5 * u1_filter .+ T_32_5 * u2_filter .+ T_33_5 * u3_filter) * analy_sol, Nx, Ny, Nz)
s22 = reshape(HI_5 * e_5 * g₃⁵[:], Nx, Ny, Nz)

s23 = reshape(HI_6 * e_6 * e_6T * (T_31_6 * u1_filter .+ T_32_6 * u2_filter .+ T_33_6 * u3_filter) * analy_sol, Nx, Ny, Nz)
s24 = reshape(HI_6 * e_6 * g₃⁶[:], Nx, Ny, Nz)


s25 = reshape((T_11_1' .- Z_11_1') * (HI_1 * e_1 * (e_1T)) * u1_filter * analy_sol, Nx, Ny, Nz)
s26 = reshape((T_11_1' .- Z_11_1') * (HI_1 * e_1 * g₁¹[:]), Nx, Ny, Nz)

s27 = reshape((T_21_1' .- Z_21_1') * (HI_1 * e_1 * (e_1T)) * u2_filter * analy_sol, Nx, Ny, Nz)
s28 = reshape((T_21_1' .- Z_21_1') * (HI_1 * e_1 * g₂¹[:]), Nx, Ny, Nz)

s29 = reshape((T_31_1' .- Z_31_1') * (HI_1 * e_1 * (e_1T)) * u3_filter * analy_sol, Nx, Ny, Nz)
s30 = reshape((T_31_1' .- Z_31_1') * (HI_1 * e_1 * g₃¹[:]), Nx, Ny, Nz)

s31 = reshape((T_11_2' .- Z_11_2') * (HI_2 * e_2 * (e_2T)) * u1_filter * analy_sol, Nx, Ny, Nz)
s32 = reshape((T_11_2' .- Z_11_2') * (HI_2 * e_2 * g₁²[:]), Nx, Ny, Nz)

s33 = reshape((T_21_2' .- Z_21_2') * (HI_2 * e_2 * (e_2T)) * u2_filter * analy_sol, Nx, Ny, Nz)
s34 = reshape((T_21_2' .- Z_21_2') * (HI_2 * e_2 * g₂²[:]), Nx, Ny, Nz)

s35 = reshape((T_31_2' .- Z_31_2') * (HI_2 * e_2 * (e_2T)) * u3_filter * analy_sol, Nx, Ny, Nz)
s36 = reshape((T_31_2' .- Z_31_2') * (HI_2 * e_2 * g₃²[:]), Nx, Ny, Nz)

s37 = reshape((T_12_1' .- Z_12_1') * (HI_1 * e_1 * (e_1T)) * u1_filter * analy_sol, Nx, Ny, Nz)
s38 = reshape((T_12_1' .- Z_12_1') * (HI_1 * e_1 * g₁¹[:]), Nx, Ny, Nz)

s39 = reshape((T_22_1' .- Z_22_1') * (HI_1 * e_1 * (e_1T)) * u2_filter * analy_sol, Nx, Ny, Nz)
s40 = reshape((T_22_1' .- Z_22_1') * (HI_1 * e_1 * g₂¹[:]), Nx, Ny, Nz)

s41 = reshape((T_32_1' .- Z_32_1') * (HI_1 * e_1 * (e_1T)) * u3_filter * analy_sol, Nx, Ny, Nz)
s42 = reshape((T_32_1' .- Z_32_1') * (HI_1 * e_1 * g₃¹[:]), Nx, Ny, Nz)

s43 = reshape((T_12_2' .- Z_12_2') * (HI_2 * e_2 * (e_2T)) * u1_filter * analy_sol, Nx, Ny, Nz)
s44 = reshape((T_12_2' .- Z_12_2') * (HI_2 * e_2 * g₁²[:]), Nx, Ny, Nz)

s45 = reshape((T_22_2' .- Z_22_2') * (HI_2 * e_2 * (e_2T)) * u2_filter * analy_sol, Nx, Ny, Nz)
s46 = reshape((T_22_2' .- Z_22_2') * (HI_2 * e_2 * g₂²[:]), Nx, Ny, Nz)

s47 = reshape((T_32_2' .- Z_32_2') * (HI_2 * e_2 * (e_2T)) * u3_filter * analy_sol, Nx, Ny, Nz)
s48 = reshape((T_32_2' .- Z_31_2') * (HI_2 * e_2 * g₃²[:]), Nx, Ny, Nz)

s49 = reshape((T_13_1' .- Z_13_1') * (HI_1 * e_1 * (e_1T)) * u1_filter * analy_sol, Nx, Ny, Nz)
s50 = reshape((T_13_1' .- Z_13_1') * (HI_1 * e_1 * g₁¹[:]), Nx, Ny, Nz)

s51 = reshape((T_23_1' .- Z_23_1') * (HI_1 * e_1 * (e_1T)) * u2_filter * analy_sol, Nx, Ny, Nz)
s52 = reshape((T_23_1' .- Z_23_1') * (HI_1 * e_1 * g₂¹[:]), Nx, Ny, Nz)

s53 = reshape((T_33_1' .- Z_33_1') * (HI_1 * e_1 * (e_1T)) * u3_filter * analy_sol, Nx, Ny, Nz)
s54 = reshape((T_33_1' .- Z_33_1') * (HI_1 * e_1 * g₃¹[:]), Nx, Ny, Nz)


s55 = reshape((T_13_2' .- Z_13_2') * (HI_2 * e_2 * (e_2T)) * u1_filter * analy_sol, Nx, Ny, Nz)
s56 = reshape((T_13_2' .- Z_13_2') * (HI_2 * e_2 * g₁²[:]), Nx, Ny, Nz)

s57 = reshape((T_23_2' .- Z_23_2') * (HI_2 * e_2 * (e_2T)) * u2_filter * analy_sol, Nx, Ny, Nz)
s58 = reshape((T_23_2' .- Z_23_2') * (HI_2 * e_2 * g₂²[:]), Nx, Ny, Nz)

s59 = reshape((T_33_2' .- Z_33_2') * (HI_2 * e_2 * (e_2T)) * u3_filter * analy_sol, Nx, Ny, Nz)
s60 = reshape((T_33_2' .- Z_33_2') * (HI_2 * e_2 * g₃²[:]), Nx, Ny, Nz)

s61 = 
s62 = 

s63 = 
s64 = 