# Rate and State friction file from Brittany

# function rateandstate(V2, V3, psi, σn, τ2, τ3, η, a, V0)

#     V = sqrt(V2^2 + V3^2)
#     dV_dV2 = 0.5*(V2^2 + V3^2)^(-0.5) .* 2 .* V2
#     dV_dV3 = 0.5*(V2^2 + V3^2)^(-0.5) .* 2 .* V3
     
#     Y = (1 ./ (2 .* V0)) .* exp.(psi ./ a)
#     f = a .* asinh.(V .* Y)  # compute friction coefficient 
#     df_dV2  = a .* (1 ./ sqrt.(1 + (V .* Y).^2)) .* (dV_dV2 .* Y)  # derivative wrt V_2
#     df_dV3  = a .* (1 ./ sqrt.(1 + (V .* Y).^2)) .* (dV_dV3 .* Y)  # derivative wrt V_2

#     g1 = σn .* f .* V2 / V   + η .* V2 - τ2
#     g2 = σn .* f .* V3 / V   + η .* V3 - τ3
   

#     A2 = V2/V
#     A3 = V3/V

#     dA2_dV2 = (V - V2*dV_dV2)/V^2
#     dA2_dV3 = (-V2*dV_dV3)/V^2
#     dA3_dV2 = (-V3*dV_dV2)/V^2
#     dA3_dV3 = (V - V3*dV_dV3)/V^2

#     dg1_dV2 = σn .* (df_dV2 .* V2 / V + f .* dA2_dV2) + η
#     dg1_dV3 = σn .* (df_dV3 .* V2 / V + f .* dA2_dV3)

#     dg2_dV2 = σn .* (df_dV2 .* V3 / V + f .* dA3_dV2) 
#     dg2_dV3 = σn .* (df_dV3 .* V3 / V + f .* dA3_dV3) + η 

#     return (g1, g2, dg1_dV2, dg1_dV3, dg2_dV2, dg2_dV3)
# end
  

# function newtbndv(func, x, y; ftol = 1e-12, maxiter = 500, 
#                     atolx = 1e-4, rtolx = 1e-4)

#     (f, g, dfx, dfy, dgx, dgy) = func(x, y)
#     for iter = 1:maxiter

#         z = [x; y] 
#         (f, g, dfx, dfy, dgx, dgy) = func(x, y)
#         if (abs(f) < ftol && abs(g) < ftol)
#             return (x, y, f, g, iter)
#         end
       
      
#         J = [dfx dfy; dgx dgy] 
#         dx, dy = -J\[f; g]
      
#         x = x + dx
#         y = y + dy

#         @show dx, atolx + rtolx * (abs(dx) + abs(x))
#         @show dy, atolx + rtolx * (abs(dy) + abs(y))
  
#        if abs(dx) < atolx + rtolx * (abs(dx) + abs(x)) && abs(dy) < atolx + rtolx * (abs(dy) + abs(y))
            
#             return (x, y, f, g, iter)
#        end
#     end
#     return (x, y, f, g, -maxiter)
# end



function rateandstate(V2, V3, psi, σn, τ2, τ3, η, a, V0)

    V = sqrt(V2^2 + V3^2)
    dV_dV2 = 0.5*(V2^2 + V3^2)^(-0.5) .* 2 .* V2
    dV_dV3 = 0.5*(V2^2 + V3^2)^(-0.5) .* 2 .* V3
    
    Y = (1 ./ (2 .* V0)) .* exp.(psi ./ a)
    f = a .* asinh.(V .* Y)  # compute friction coefficient 
    df_dV2  = a .* (1 ./ sqrt.(1 + (V .* Y).^2)) .* (dV_dV2 .* Y)  # derivative wrt V_2
    df_dV3  = a .* (1 ./ sqrt.(1 + (V .* Y).^2)) .* (dV_dV3 .* Y)  # derivative wrt V_2

    g1 = σn .* f .* V2 / V   + η .* V2 - τ2
    g2 = σn .* f .* V3 / V   + η .* V3 - τ3
   

    A2 = V2/V
    A3 = V3/V

    dA2_dV2 = (V - V2*dV_dV2)/V^2
    dA2_dV3 = (-V2*dV_dV3)/V^2
    dA3_dV2 = (-V3*dV_dV2)/V^2
    dA3_dV3 = (V - V3*dV_dV3)/V^2

    dg1_dV2 = σn .* (df_dV2 .* V2 / V + f .* dA2_dV2) + η
    dg1_dV3 = σn .* (df_dV3 .* V2 / V + f .* dA2_dV3)

    dg2_dV2 = σn .* (df_dV2 .* V3 / V + f .* dA3_dV2) 
    dg2_dV3 = σn .* (df_dV3 .* V3 / V + f .* dA3_dV3) + η 

    return (g1, g2, dg1_dV2, dg1_dV3, dg2_dV2, dg2_dV3)
end
  

function newtbndv(func, x, y; ftol = 1e-12, maxiter = 500, 
                    atolx = 1e-4, rtolx = 1e-4)

    (f, g, dfx, dfy, dgx, dgy) = func(x, y)
    for iter = 1:maxiter

        z = [x; y] 
        (f, g, dfx, dfy, dgx, dgy) = func(x, y)
      
        J = [dfx dfy; dgx dgy] 
        dx, dy = -J\[f; g]
      
        x = x + dx
        y = y + dy
  
       if abs(f) < ftol && abs(dx) < atolx + rtolx * (abs(dx) + abs(x)) && abs(g) < ftol && abs(dy) < atolx + rtolx * (abs(dy) + abs(y))
            
            return (x, y, f, g, iter)
       end
    end
    return (x, y, f, g, -maxiter)
end

function rateandstate_vectorized(V2, V3, psi_v, σn, τ2_v, τ3_v, η, a_v, V0)
    V  = sqrt.(V2.^2 .+ V3.^2) 
    dV_dV2 = 0.5 .* sqrt.(V2.^2 .+ V3.^2) .* 2 .* V2
    dV_dV3 = 0.5 .* sqrt.(V2.^2 .+ V3.^2) .* 2 .* V3

    Y = (1 ./ (2 .* V0)) .* exp.(psi_v ./ 2)
    f = a_v .* asinh.(V .* Y)
    df_dV2 = a_v .* (1 ./ sqrt.(1 .+ (V .* Y).^2)) .* (dV_dV2 .* Y)
    df_dV3 = a_v .* (1 ./ sqrt.(1 .+ (V .* Y).^2)) .* (dV_dV3 .* Y)

    g1 = σn .* f .* V2 ./ V .+ η .* V2 .- τ2_v
    g2 = σn .* f .* V3 ./ V .+ η .* V3 .- τ3_v

    dA2_dV2 = (V .- V2 .* dV_dV2) ./ V.^2
    dA2_dV3 = (-V2 .* dV_dV3) ./ V.^2
    dA3_dV2 = (-V3 .* dV_dV2) ./ V.^2
    dA3_dV3 = (V .- V3 .* dV_dV3) ./ V.^2

    dg1_dV2 = σn .* (df_dV2 .* V2 ./ V .+ f .* dA2_dV2) .+ η
    dg1_dV3 = σn .* (df_dV3 .* V2 ./ V .+ f .* dA2_dV3)

    dg2_dV2 = σn .* (df_dV3 .* V3 ./ V .+ f .* dA3_dV2)
    dg2_dV3 = σn .* (df_dV3 .* V3 ./ V .+ f .* dA3_dV3) .+ η

    return (g1, g2, dg1_dV2, dg1_dV3, dg2_dV2, dg2_dV3)
end

function newtbndv_vectorized(rateandstate_vectorized, V2, V3; 
                            ftol=1e-12, maxiter=500, atolx = 1e-4, rtolx=1e-4)
    for iter = 1:maxiter
        (f_v, g_v, dfx_v, dfy_v, dgx_v, dgy_v) = rateandstate_vectorized(V2, V3,  
                                        psi_v, σn, τ2_v, τ3_v, η, a_v, V0)
        inv_J = map_jacobian_inv.(dfx_v, dfy_v, dgx_v, dgy_v)
        dV2, dV3 = -inv_J .* map_cat.(f_v, g_v)

        V2 .+= df_v
        V3 .+= dg_v
        
        # TODO vectorized control flow
        if all(abs(f_v) .< ftol) && all(abs(df_v) .< atolx .+ rtolx .* (abs.(df_v) + abs(x))) && abs(g) < ftol && abs(dy) < atolx + rtolx * (abs(dy) + abs(y))
    end
end
# TESTING:
ψn = 0.6
an = 0.015
η = 32/6
σn = 50
RSV0 = 1e-6
V1_actual = 1
V2_actual = 2
Vactual = sqrt(V1_actual^2 + V2_actual^2)
τ = σn * an * asinh(Vactual/(2*RSV0) * exp(ψn/an)) * V1_actual/Vactual + η*V1_actual
τz =  σn * an * asinh(Vactual/(2*RSV0) * exp(ψn/an)) * V2_actual/Vactual + η*V2_actual
Vn1 = 2              #initial guess
Vn2 = 1      #initial guess

obj_rs(V2, V3) = rateandstate(V2, V3, ψn, σn, τ, τz, η, an, RSV0)
(Vn2, Vn3, f, g, iter) = newtbndv(obj_rs, Vn1, Vn2; ftol = 1e-12,
                            atolx = 1e-12, rtolx = 1e-12)


# Testing vectorized 
V2=[1,2,3]
V3=[1,3,4]
psi_v = [0.6, 0.6, 0.6]
a_v = [0.015, 0.015, 0.6]
τ2_v = [τ, τ, τ]
τ3_v = [τz, τz, τz]

rateandstate_vectorized(V2, V3, psi_v, σn, τ2_v, τ3_v, η, a_v, V0)


# # TESTING:
# ψn = 0.6
# an = 0.015
# η = 32/6
# σn = 50
# RSV0 = 1e-6
# V1_actual = 1
# V2_actual = 2
# Vactual = sqrt(V1_actual^2 + V2_actual^2)
# τ = σn * an * asinh(Vactual/(2*RSV0) * exp(ψn/an)) * V1_actual/Vactual + η*V1_actual
# τz =  σn * an * asinh(Vactual/(2*RSV0) * exp(ψn/an)) * V2_actual/Vactual + η*V2_actual
# Vn1 = V1_actual               #initial guess
# Vn2 = V2_actual + .00001      #initial guess

# obj_rs(V2, V3) = rateandstate(V2, V3, ψn, σn, τ, τz, η, an, RSV0)
# (Vn2, Vn3, f, g, iter) = newtbndv(obj_rs, Vn1, Vn2; ftol = 1e-12,
#                             atolx = 1e-12, rtolx = 1e-12)


function map_jacobian(x, y, a, b)
    return [x y; a b]
end

function map_jacobian_inv(x, y, a, b)
    return inv([x y; a b])
end

function map_cat(x,y)
    return [x;y]
end