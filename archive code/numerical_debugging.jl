
# SAT_1_LHS = -HI_tilde * (
#         e_3 * H_3 * e_3T * (T_11_3 * u1_filter .+ T_12_3 * u2_filter .+ T_13_3 * u3_filter)
#     +   e_4 * H_4 * e_4T * (T_11_4 * u1_filter .+ T_12_4 * u2_filter .+ T_13_4 * u3_filter)
#     +   e_5 * H_5 * e_5T * (T_11_5 * u1_filter .+ T_12_5 * u2_filter .+ T_13_5 * u3_filter)
#     +   e_6 * H_6 * e_6T * (T_11_6 * u1_filter .+ T_12_6 * u2_filter .+ T_13_6 * u3_filter)
# ) 

# SAT_2_LHS = - HI_tilde * (
#         e_3 * H_3 * e_3T * (T_21_3 * u1_filter .+ T_22_3 * u2_filter .+ T_23_3 * u3_filter)
#     +   e_4 * H_4 * e_4T * (T_21_4 * u1_filter .+ T_22_4 * u2_filter .+ T_23_4 * u3_filter)
#     +   e_5 * H_5 * e_5T * (T_21_5 * u1_filter .+ T_22_5 * u2_filter .+ T_23_5 * u3_filter)
#     +   e_6 * H_5 * e_6T * (T_21_6 * u1_filter .+ T_22_6 * u2_filter .+ T_23_6 * u3_filter)
# ) 


# SAT_3_LHS = - HI_tilde * (
#         e_3 * H_3 * e_3T * (T_31_3 * u1_filter .+ T_32_3 * u2_filter .+ T_33_3 * u3_filter)
#     +   e_4 * H_4 * e_4T * (T_31_4 * u1_filter .+ T_32_4 * u2_filter .+ T_33_4 * u3_filter)
#     +   e_5 * H_5 * e_5T * (T_31_5 * u1_filter .+ T_32_5 * u2_filter .+ T_33_5 * u3_filter)
#     +   e_6 * H_6 * e_6T * (T_31_6 * u1_filter .+ T_32_6 * u2_filter .+ T_33_6 * u3_filter)
# )



# SAT_tilde_1_LHS = - HI_tilde * (
#         (T_11_1 .- Z_11_1)' * (e_1 * H_1 * (e_1T)) * u1_filter
#     +   (T_21_1 .- Z_21_1)' * (e_1 * H_1 * (e_1T)) * u2_filter
#     +   (T_31_1 .- Z_31_1)' * (e_1 * H_1 * (e_1T)) * u3_filter
#     +   (T_11_2 .- Z_11_2)' * (e_2 * H_2 * (e_2T)) * u1_filter
#     +   (T_21_2 .- Z_21_2)' * (e_2 * H_2 * (e_2T)) * u2_filter
#     +   (T_31_2 .- Z_31_2)' * (e_2 * H_2 * (e_2T)) * u3_filter
# )

# SAT_tilde_2_LHS = - HI_tilde * (
#         (T_12_1 .- Z_12_1)' * (e_1 * H_1 * (e_1T)) * u1_filter
#     +   (T_22_1 .- Z_22_1)' * (e_1 * H_1 * (e_1T)) * u2_filter
#     +   (T_32_1 .- Z_32_1)' * (e_1 * H_1 * (e_1T)) * u3_filter
#     +   (T_12_2 .- Z_12_2)' * (e_2 * H_2 * (e_2T)) * u1_filter
#     +   (T_22_2 .- Z_22_2)' * (e_2 * H_2 * (e_2T)) * u2_filter
#     +   (T_32_2 .- Z_32_2)' * (e_2 * H_2 * (e_2T)) * u3_filter
# )

# SAT_tilde_3_LHS = - HI_tilde * (
#         (T_13_1 .- Z_13_1)' * (e_1 * H_1 * (e_1T)) * u1_filter
#     +   (T_23_1 .- Z_23_1)' * (e_1 * H_1 * (e_1T)) * u2_filter
#     +   (T_33_1 .- Z_33_1)' * (e_1 * H_1 * (e_1T)) * u3_filter
#     +   (T_13_2 .- Z_13_2)' * (e_2 * H_2 * (e_2T)) * u1_filter
#     +   (T_23_2 .- Z_23_2)' * (e_2 * H_2 * (e_2T)) * u2_filter
#     +   (T_33_2 .- Z_33_2)' * (e_2 * H_2 * (e_2T)) * u3_filter
# )


# ### Assembling SBP terms for right-hand-side (RHS) traction condition
# SAT_1_RHS = - HI_tilde * (
#         e_3 * H_3 * g₁³[:]
#     +   e_4 * H_4 * g₁⁴[:]
#     +   e_5 * H_5 * g₁⁵[:]
#     +   e_6 * H_6 * g₁⁶[:]
# )

# SAT_2_RHS = - HI_tilde * (
#         e_3 * H_3 * g₂³[:]
#     +   e_4 * H_4 * g₂⁴[:]
#     +   e_5 * H_5 * g₂⁵[:]
#     +   e_6 * H_6 * g₂⁶[:]
# )

# SAT_3_RHS = - HI_tilde * (
#         e_3 * H_3 * g₃³[:]
#     +   e_4 * H_4 * g₃⁴[:]
#     +   e_5 * H_5 * g₃⁵[:]
#     +   e_6 * H_6 * g₃⁶[:]
# )



# ### Assembling SBP terms for right-hand-side (RHS) Dirichlet condition
# SAT_tilde_1_RHS = - HI_tilde * (
#         (T_11_1 .- Z_11_1)' * (e_1 * H_1 * g₁¹[:])
#     +   (T_21_1 .- Z_21_1)' * (e_1 * H_1 * g₂¹[:])
#     +   (T_31_1 .- Z_31_1)' * (e_1 * H_1 * g₃¹[:])
#     +   (T_11_2 .- Z_11_2)' * (e_2 * H_2 * g₁²[:])
#     +   (T_21_2 .- Z_21_2)' * (e_2 * H_2 * g₂²[:])
#     +   (T_31_2 .- Z_31_2)' * (e_2 * H_2 * g₃²[:])
# )

# SAT_tilde_2_RHS = - HI_tilde * (
#         (T_12_1 .- Z_12_1)' * (e_1 * H_1 * g₁¹[:])
#     +   (T_22_1 .- Z_22_1)' * (e_1 * H_1 * g₂¹[:])
#     +   (T_32_1 .- Z_32_1)' * (e_1 * H_1 * g₃¹[:])
#     +   (T_12_2 .- Z_12_2)' * (e_2 * H_2 * g₁²[:])
#     +   (T_22_2 .- Z_22_2)' * (e_2 * H_2 * g₂²[:])
#     +   (T_32_2 .- Z_31_2)' * (e_2 * H_2 * g₃²[:])
# )

# SAT_tilde_3_RHS = - HI_tilde * (
#         (T_13_1 .- Z_13_1)' * (e_1 * H_1 * g₁¹[:])
#     +   (T_23_1 .- Z_23_1)' * (e_1 * H_1 * g₂¹[:])
#     +   (T_33_1 .- Z_33_1)' * (e_1 * H_1 * g₃¹[:])
#     +   (T_13_2 .- Z_13_2)' * (e_2 * H_2 * g₁²[:])
#     +   (T_23_2 .- Z_23_2)' * (e_2 * H_2 * g₂²[:])
#     +   (T_33_2 .- Z_33_2)' * (e_2 * H_2 * g₃²[:])
# )


# Manually debugging
# SAT_1 
s1 = reshape(e_3 * H_3 * e_3T * (T_11_3 * u1_filter .+ T_12_3 * u2_filter .+ T_13_3 * u3_filter) * analy_sol, Nx, Ny, Nz)
s2 = reshape(e_3 * H_3 * g₁³[:], Nx,Ny,Nz)

hcat(s1,s2)

s3 = reshape(e_4 * H_4 * e_4T * (T_11_4 * u1_filter .+ T_12_4 * u2_filter .+ T_13_4 * u3_filter) * analy_sol, Nx, Ny, Nz)
s4 = reshape(e_4 * H_4 * g₁⁴[:], Nx, Ny, Nz)
hcat(s3, s4)

s5 = reshape(e_5 * H_5 * e_5T * (T_11_5 * u1_filter .+ T_12_5 * u2_filter .+ T_13_5 * u3_filter) * analy_sol, Nx, Ny, Nz)
s6 = reshape(e_5 * H_5 * g₁⁵[:], Nx, Ny, Nz)
s5[:,:,1]
s6[:,:,1]

s7 = reshape(e_6 * H_6 * e_6T * (T_11_6 * u1_filter .+ T_12_6 * u2_filter .+ T_13_6 * u3_filter) *analy_sol, Nx, Ny, Nz)
s8 = reshape(e_6 * H_6 * g₁⁶[:], Nx, Ny, Nz)
s7[:,:,end]
s8[:,:,end]


# SAT_2
s9 = reshape(e_3 * H_3 * e_3T * (T_21_3 * u1_filter .+ T_22_3 * u2_filter .+ T_23_3 * u3_filter) * analy_sol, Nx, Ny, Nz)
s10 = reshape(e_3 * H_3 * g₂³[:], Nx, Ny, Nz)
hcat(s9,s10)

s11 = reshape(e_4 * H_4 * e_4T * (T_21_4 * u1_filter .+ T_22_4 * u2_filter .+ T_23_4 * u3_filter) * analy_sol, Nx, Ny, Nz)
s12 = reshape( e_4 * H_4 * g₂⁴[:], Nx, Ny, Nz)
hcat(s11, s12)

s13 = reshape(e_5 * H_5 * e_5T * (T_21_5 * u1_filter .+ T_22_5 * u2_filter .+ T_23_5 * u3_filter) * analy_sol, Nx, Ny, Nz)
s14 = reshape(e_5 * H_5 * g₂⁵[:], Nx, Ny, Nz)
# s13 and s14 should be both zero Arrays

s15 = reshape(e_6 * H_5 * e_6T * (T_21_6 * u1_filter .+ T_22_6 * u2_filter .+ T_23_6 * u3_filter) * analy_sol, Nx, Ny, Nz)
s16 = reshape(e_6 * H_6 * g₂⁶[:], Nx, Ny, Nz)
# s15 and s16 should be both zero Arrays 

# SAT_3
s17 = reshape(e_3 * H_3 * e_3T * (T_31_3 * u1_filter .+ T_32_3 * u2_filter .+ T_33_3 * u3_filter) * analy_sol, Nx, Ny, Nz)
s18 = reshape(e_3 * H_3 * g₃³[:], Nx, Ny, Nz)
# s17 and s18 should be both zero Arrays

s19 = reshape(e_4 * H_4 * e_4T * (T_31_4 * u1_filter .+ T_32_4 * u2_filter .+ T_33_4 * u3_filter) * analy_sol, Nx, Ny, Nz)
s20 = reshape(e_4 * H_4 * g₃⁴[:], Nx, Ny, Nz)
# s19 and s20 should be both zero Arrays

s21 = reshape(e_5 * H_5 * e_5T * (T_31_5 * u1_filter .+ T_32_5 * u2_filter .+ T_33_5 * u3_filter) * analy_sol, Nx, Ny, Nz)
s22 = reshape(e_5 * H_5 * g₃⁵[:], Nx, Ny, Nz)
s21[:,:,1]
s22[:,:,1]

s23 = reshape(e_6 * H_6 * e_6T * (T_31_6 * u1_filter .+ T_32_6 * u2_filter .+ T_33_6 * u3_filter) * analy_sol, Nx, Ny, Nz)
s24 = reshape(e_6 * H_6 * g₃⁶[:], Nx, Ny, Nz)
s23[:,:,end]
s24[:,:,end]


# SAT_tilde_1 
s25 = reshape((T_11_1 .- Z_11_1)' * (e_1 * H_1 * (e_1T)) * u1_filter * analy_sol, Nx, Ny, Nz)
s26 = reshape((T_11_1 .- Z_11_1)' * (e_1 * H_1 * g₁¹[:]), Nx, Ny, Nz)

s27 = reshape((T_21_1 .- Z_21_1)' * (e_1 * H_1 * (e_1T)) * u2_filter * analy_sol, Nx, Ny, Nz)
s28 = reshape((T_21_1 .- Z_21_1)' * (e_1 * H_1 * g₂¹[:]), Nx, Ny, Nz)

s29 = reshape((T_31_1 .- Z_31_1)' * (e_1 * H_1 * (e_1T)) * u3_filter * analy_sol, Nx, Ny, Nz)
s30 = reshape((T_31_1 .- Z_31_1)' * (e_1 * H_1 * g₃¹[:]), Nx, Ny, Nz)

s31 = reshape((T_11_2 .- Z_11_2)' * (e_2 * H_2 * (e_2T)) * u1_filter * analy_sol, Nx, Ny, Nz)
s32 = reshape((T_11_2 .- Z_11_2)' * (e_2 * H_2 * g₁²[:]), Nx, Ny, Nz)

s33 = reshape((T_21_2 .- Z_21_2)' * (e_2 * H_2 * (e_2T)) * u2_filter * analy_sol, Nx, Ny, Nz)
s34 = reshape((T_21_2 .- Z_21_2)' * (e_2 * H_2 * g₂²[:]), Nx, Ny, Nz)

s35 = reshape((T_31_2 .- Z_31_2)' * (e_2 * H_2 * (e_2T)) * u3_filter * analy_sol, Nx, Ny, Nz)
s36 = reshape((T_31_2 .- Z_31_2)' * (e_2 * H_2 * g₃²[:]), Nx, Ny, Nz)
# SAT_tilde_2 
s37 = reshape((T_12_1 .- Z_12_1)' * (e_1 * H_1 * (e_1T)) * u1_filter * analy_sol, Nx, Ny, Nz)
s38 = reshape((T_12_1 .- Z_12_1)' * (e_1 * H_1 * g₁¹[:]), Nx, Ny, Nz)

s39 = reshape((T_22_1 .- Z_22_1)' * (e_1 * H_1 * (e_1T)) * u2_filter * analy_sol, Nx, Ny, Nz)
s40 = reshape((T_22_1 .- Z_22_1)' * (e_1 * H_1 * g₂¹[:]), Nx, Ny, Nz)

s41 = reshape((T_32_1 .- Z_32_1)' * (e_1 * H_1 * (e_1T)) * u3_filter * analy_sol, Nx, Ny, Nz)
s42 = reshape((T_32_1 .- Z_32_1)' * (e_1 * H_1 * g₃¹[:]), Nx, Ny, Nz)

s43 = reshape((T_12_2 .- Z_12_2)' * (e_2 * H_2 * (e_2T)) * u1_filter * analy_sol, Nx, Ny, Nz)
s44 = reshape((T_12_2 .- Z_12_2)' * (e_2 * H_2 * g₁²[:]), Nx, Ny, Nz)

s45 = reshape((T_22_2 .- Z_22_2)' * (e_2 * H_2 * (e_2T)) * u2_filter * analy_sol, Nx, Ny, Nz)
s46 = reshape((T_22_2 .- Z_22_2)' * (e_2 * H_2 * g₂²[:]), Nx, Ny, Nz)

s47 = reshape((T_33_2 .- Z_33_2)' * (e_2 * H_2 * (e_2T)) * u3_filter * analy_sol, Nx, Ny, Nz)
s48 = reshape((T_32_2 .- Z_31_2)' * (e_2 * H_2 * g₃²[:]), Nx, Ny, Nz)


# SAT_tilde_3 
s49 = reshape((T_13_1 .- Z_13_1)' * (e_1 * H_1 * (e_1T)) * u1_filter * analy_sol, Nx, Ny, Nz)
s50 = reshape((T_13_1 .- Z_13_1)' * (e_1 * H_1 * g₁¹[:]), Nx, Ny, Nz)

s51 = reshape((T_23_1 .- Z_23_1)' * (e_1 * H_1 * (e_1T)) * u2_filter * analy_sol, Nx, Ny, Nz)
s52 = reshape((T_23_1 .- Z_23_1)' * (e_1 * H_1 * g₂¹[:]), Nx, Ny, Nz)

s53 = reshape((T_33_1 .- Z_33_1)' * (e_1 * H_1 * (e_1T)) * u3_filter * analy_sol, Nx, Ny, Nz)
s54 = reshape((T_33_1 .- Z_33_1)' * (e_1 * H_1 * g₃¹[:]), Nx, Ny, Nz)

s55 = reshape((T_13_2 .- Z_13_2)' * (e_2 * H_2 * (e_2T)) * u1_filter * analy_sol, Nx, Ny, Nz)
s56 = reshape((T_13_2 .- Z_13_2)' * (e_2 * H_2 * g₁²[:]), Nx, Ny, Nz)

s57 = reshape((T_23_2 .- Z_23_2)' * (e_2 * H_2 * (e_2T)) * u2_filter * analy_sol, Nx, Ny, Nz)
s58 = reshape((T_23_2 .- Z_23_2)' * (e_2 * H_2 * g₂²[:]), Nx, Ny, Nz)

s59 = reshape((T_33_2 .- Z_33_2)' * (e_2 * H_2 * (e_2T)) * u3_filter * analy_sol, Nx, Ny, Nz)
s60 = reshape((T_33_2 .- Z_33_2)' * (e_2 * H_2 * g₃²[:]), Nx, Ny, Nz)


reshape(SAT_1_LHS * analy_sol, Nx, Ny, Nz)
reshape(SAT_1_RHS, Nx, Ny, Nz)

reshape(SAT_2_LHS * analy_sol, Nx, Ny, Nz)
reshape(SAT_2_RHS, Nx, Ny, Nz)

reshape(SAT_3_LHS * analy_sol, Nx, Ny, Nz)
reshape(SAT_3_RHS, Nx, Ny, Nz)

reshape(SAT_tilde_1_LHS * analy_sol, Nx, Ny, Nz)
reshape(SAT_tilde_1_RHS, Nx, Ny, Nz)



# Debugging the SPD of the matrix M


p2_px2_result = Matrix(H_tilde * p2_px2)
Matrix(H_tilde * BS_Front)

H_tilde * End_operator' * HI_1 * End_operator * BS_Front

H_tilde * Front_operator' * HI_1 * Front_operator * BS_End


BS_end_result = H_tilde * kron(I_Nz,I_Ny,HIx) * End_operator' * End_operator * BS_End
BS_front_result = H_tilde * kron(I_Nz,I_Ny,HIx) * Front_operator' * Front_operator * BS_Front


BS_end_result'
BS_front_result'


