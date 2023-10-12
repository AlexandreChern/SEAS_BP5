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