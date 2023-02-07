include("elasticity_3D.jl")


SAT_1_LHS = (
        u1_filter' * (- HI_tilde *  (   e_3 * H_3 * e_3T * T_11_3 
                                    .+  e_4 * H_4 * e_4T * T_11_4
                                    .+  e_5 * H_5 * e_5T * T_11_5
                                    .+  e_6 * H_6 * e_6T * T_11_6) * u1_filter)

    +   u2_filter' * (- HI_tilde *  (   e_3 * H_3 * e_3T * T_12_3
                                    .+  e_4 * H_4 * e_4T * T_12_4
                                    .+  e_5 * H_5 * e_5T * T_12_5
                                    .+  e_6 * H_6 * e_6T * T_12_6) * u2_filter)

    +   u3_filter' * (- HI_tilde *  (   e_3 * H_3 * e_3T * T_13_3
                                    .+  e_4 * H_4 * e_4T * T_13_4
                                    .+  e_5 * H_5 * e_5T * T_13_5
                                    .+  e_6 * H_6 * e_6T * T_13_6) * u3_filter)
) 


SAT_2_LHS = (
        u1_filter' * (- HI_tilde *  (   e_3 * H_3 * e_3T * T_21_3 
                                    .+  e_4 * H_4 * e_4T * T_21_4
                                    .+  e_5 * H_5 * e_5T * T_21_5
                                    .+  e_6 * H_6 * e_6T * T_21_6) * u1_filter)

    +   u2_filter' * (- HI_tilde *  (   e_3 * H_3 * e_3T * T_22_3
                                    .+  e_4 * H_4 * e_4T * T_22_4
                                    .+  e_5 * H_5 * e_5T * T_22_5
                                    .+  e_6 * H_6 * e_6T * T_22_6) * u2_filter)

    +   u3_filter' * (- HI_tilde *  (   e_3 * H_3 * e_3T * T_23_3
                                    .+  e_4 * H_4 * e_4T * T_23_4
                                    .+  e_5 * H_5 * e_5T * T_23_5
                                    .+  e_6 * H_6 * e_6T * T_23_6) * u3_filter)
) 


SAT_2_LHS = (
        u1_filter' * (- HI_tilde *  (   e_3 * H_3 * e_3T * T_31_3 
                                    .+  e_4 * H_4 * e_4T * T_31_4
                                    .+  e_5 * H_5 * e_5T * T_31_5
                                    .+  e_6 * H_6 * e_6T * T_31_6) * u1_filter)

    +   u2_filter' * (- HI_tilde *  (   e_3 * H_3 * e_3T * T_32_3
                                    .+  e_4 * H_4 * e_4T * T_32_4
                                    .+  e_5 * H_5 * e_5T * T_32_5
                                    .+  e_6 * H_6 * e_6T * T_32_6) * u2_filter)

    +   u3_filter' * (- HI_tilde *  (   e_3 * H_3 * e_3T * T_33_3
                                    .+  e_4 * H_4 * e_4T * T_33_4
                                    .+  e_5 * H_5 * e_5T * T_33_5
                                    .+  e_6 * H_6 * e_6T * T_33_6) * u3_filter)
) 




### Assembling SBP terms for Dirichlet 

SAT_tilde_1_LHS = - HI_tilde * (
        (T_11_1 .- Z_11_1)' * (e_1 * H_1 * (e_1T)) * u1_filter
    +   (T_21_1 .- Z_21_1)' * (e_1 * H_1 * (e_1T)) * u2_filter
    +   (T_31_1 .- Z_31_1)' * (e_1 * H_1 * (e_1T)) * u3_filter
    +   (T_11_2 .- Z_11_2)' * (e_2 * H_2 * (e_2T)) * u1_filter
    +   (T_21_2 .- Z_21_2)' * (e_2 * H_2 * (e_2T)) * u2_filter
    +   (T_31_2 .- Z_31_2)' * (e_2 * H_2 * (e_2T)) * u3_filter
)


SAT_tilde_1_LHS = (
        u1_filter' * (- HI_tilde * (    
                (T_11_1 .- Z_11_1)' * (e_1 * H_1 * e_1T)
            .+  (T_11_2 .- Z_11_2)' *(e_2 * H_2 * e_2T)
            ) ) * u1_filter
    +   u2_filter' * (- HI_tilde * (
                (T_21_1 .- Z_21_1)' * (e_1 * H_1 * e_1T)
            .+  (T_21_1 .- Z_21_1)' * (e_2 * H_2 * e_2T)
            )) * u2_filter
    +   u3_filter' * (- HI_tilde * (
                (T_31_1 .- Z_31_1)' * (e_1 * H_1 * e_1T)
            .+  (T_31_1 .- Z_31_1)'
            )) * u3_filter
        # TO DO
)

SAT_tilde_2_LHS = (
        u1_filter' * (- HI_tilde * (
                (T_12_1 .- Z_12_1)' * (e_1 * H_1 * e_1T)
            .+  (T_12_2 .- Z_12_2)' * (e_2 * H_2 * e_2T)
        )) * u1_filter
    +   u2_filter' * (- HI_tilde * (
                (T_22_2 .- Z_22_2)' * (e_1 * H_1 * e_1T)
            .+  (T_22_2 .- Z_22_2)' * (e_2 * H_2 * e_2T)
        )) * u2_filter
    +   u3_filter' * (- HI_tilde * (
                (T_32_1 .- Z_32_1)' * (e_1 * H_1 * e_1T)
            .+  (T_32_2 .- Z_32_2)' * (e_2 * H_2 * e_2T)
        )) * u3_filter
)


SAT_tilde_2_LHS = (
        u1_filter' * (- HI_tilde * (
                (T_13_1 .- Z_13_1)' * (e_1 * H_1 * e_1T)
            .+  (T_13_2 .- Z_13_2)' * (e_2 * H_2 * e_2T)
        )) * u1_filter
    +   u2_filter' * (- HI_tilde * (
                (T_23_2 .- Z_23_2)' * (e_1 * H_1 * e_1T)
            .+  (T_23_2 .- Z_23_2)' * (e_2 * H_2 * e_2T)
        )) * u2_filter
    +   u3_filter' * (- HI_tilde * (
                (T_33_1 .- Z_33_1)' * (e_1 * H_1 * e_1T)
            .+  (T_33_2 .- Z_33_2)' * (e_2 * H_2 * e_2T)
        )) * u3_filter
)