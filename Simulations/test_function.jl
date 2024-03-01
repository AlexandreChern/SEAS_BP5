# This is a test file to test the Newton solver
# Will be archived later

function tau(V;a=0.04)
    ψ = 0.8072326581844863
    Y = exp(ψ / a) / (2 * BP5_coeff.V0)
    return a * asinh(V * Y)
end


xs = range(-1e-10, 1e-10, length=1000)
ys = tau.(xs)
plot(xs, ys)


xs2 = range(-1e-10,1e-10,length=1000)
ys2 = tau.(xs2; a=0.004)
plot(xs2, ys2)


using Optim

RS_filter_2D_nzind
VW_favorable_filter_RS_nzind

τ2_f = τ2[VW_favorable_filter_RS_nzind]
τ3_f = τ3[VW_favorable_filter_RS_nzind]
V2_f = V2[VW_favorable_filter_RS_nzind]
V3_f = V3[VW_favorable_filter_RS_nzind]
RSas_f = RSas[VW_favorable_filter_RS_nzind]
ψ_f = ψ[VW_favorable_filter_RS_nzind]
f_f, g_f, dfx_f, dfy_f, dgx_f, dgy_f = rateandstate_vectorized(V2_f, V3_f, ψ_f, σn, τ2_f, τ3_f, η, RSas_f, RSV0)
inv_F = map_jacobian_inv.(dfx_f, dfy_f, dgx_f, dgy_f)
V2_f_tmp, V3_f_tmp, f_f_tmp, g_f_tmp, maxiter = newtbndv_vectorized(
    rateandstate_vectorized, V2_f, V3_f, ψ_f, σn, τ2_f, τ3_f, η, RSas_f, RSV0; maxiter=1)