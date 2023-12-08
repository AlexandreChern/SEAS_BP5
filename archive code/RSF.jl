# Rate and State friction file from Brittany

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
        if (abs(f) < ftol && abs(g) < ftol)
            return (x, y, f, g, iter)
        end
       
      
        J = [dfx dfy; dgx dgy] 
        dx, dy = -J\[f; g]
      
        x = x + dx
        y = y + dy

        @show dx, atolx + rtolx * (abs(dx) + abs(x))
        @show dy, atolx + rtolx * (abs(dy) + abs(y))
  
       if abs(dx) < atolx + rtolx * (abs(dx) + abs(x)) && abs(dy) < atolx + rtolx * (abs(dy) + abs(y))
            
            return (x, y, f, g, iter)
       end
    end
    return (x, y, f, g, -maxiter)
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
Vn1 = V1_actual               #initial guess
Vn2 = V2_actual + .00001      #initial guess

obj_rs(V2, V3) = rateandstate(V2, V3, ψn, σn, τ, τz, η, an, RSV0)
(Vn2, Vn3, f, g, iter) = newtbndv(obj_rs, Vn1, Vn2; ftol = 1e-12,
                            atolx = 1e-12, rtolx = 1e-12)

