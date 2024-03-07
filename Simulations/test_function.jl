# This is a test file to test the Newton solver
# Will be archived later

function F(V_f;a=0.04)
    ψ_f = 0.8072326581844863
    Y_f = exp(ψ_f / a) / (2 * BP5_coeff.V0)
    return σn * a * asinh(V_f * Y_f)
end


function tau(V_f, τ)
    return τ - η * V_f
end

function root(V_f, τ; a = 0.004)
    return F.(V_f;a = a) .- tau.(V_f, τ)
end

xs = range(-1e-9, 1e-9, length=1000)
ys = F.(xs)
plot(xs, ys)


xs2 = range(-1e-3,1e-3,length=1000)
ys2 = F.(xs2; a=0.004)
plot(xs2, ys2)

ys2_v = tau.(xs2, τ2_f[1])
plot!(xs2, ys2_v)

plot(xs2, root(xs2, τ2_f[1]; a=0.004))

root(xs2, τ2_f[1]; a=0.004)

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
    rateandstate_vectorized, V2_f, V3_f, ψ_f, σn, τ2_f, τ3_f, η, RSas_f, RSV0; maxiter=100, α=1)


function rateandstate_vectorized(V2_f, V3_f, ψ_f, σn, τ2_f, τ3_f, η, RSas_f, RSV0)
    V_f  = sqrt.(V2_f.^2 .+ V3_f.^2) 
    dV_dV2_f = 0.5 ./ sqrt.(V2_f.^2 .+ V3_f.^2) .* 2 .* V2_f
    dV_dV3_f = 0.5 ./ sqrt.(V2_f.^2 .+ V3_f.^2) .* 2 .* V3_f

    Y_f = (1 ./ (2 .* RSV0)) .* exp.(ψ_f ./ RSas_f)
    f = RSas_f .* asinh.(V_f .* Y_f)
    df_dV2_f = RSas_f .* (1 ./ sqrt.(1 .+ (V_f .* Y_f).^2)) .* (dV_dV2_f .* Y_f)
    df_dV3_f = RSas_f .* (1 ./ sqrt.(1 .+ (V_f .* Y_f).^2)) .* (dV_dV3_f .* Y_f)

    g1 = σn .* f .* V2_f ./ V_f .+ η .* V2_f .- τ2_f
    g2 = σn .* f .* V3_f ./ V_f .+ η .* V3_f .- τ3_f

    dA2_dV2_f = (V_f .- V2_f .* dV_dV2_f) ./ V_f.^2
    dA2_dV3_f = (-V2_f .* dV_dV3_f) ./ V_f.^2
    dA3_dV2_f = (-V3_f .* dV_dV2_f) ./ V_f.^2
    dA3_dV3_f = (V_f .- V3_f .* dV_dV3_f) ./ V_f.^2

    dg1_dV2_f = σn .* (df_dV2_f .* V2_f ./ V_f .+ f .* dA2_dV2_f) .+ η
    dg1_dV3_f = σn .* (df_dV3_f .* V2_f ./ V_f .+ f .* dA2_dV3_f)

    dg2_dV2_f = σn .* (df_dV2_f .* V3_f ./ V_f .+ f .* dA3_dV2_f)
    dg2_dV3_f = σn .* (df_dV3_f .* V3_f ./ V_f .+ f .* dA3_dV3_f) .+ η

    return (g1, g2, dg1_dV2_f, dg1_dV3_f, dg2_dV2_f, dg2_dV3_f)
end

function newtbndv_vectorized(rateandstate_vectorized, V2_f, V3_f, ψ_f, σn, τ2_f, τ3_f, η, RSas_f, RSV0;  
                            ftol=1e-12, maxiter=100, atolx = 1e-4, rtolx=1e-4, α=1.0) # change atolx to 1e-8 and rtolx to 1e-8 for better stability
    (f_f, g_f, dfx_f, dfy_f, dgx_f, dgy_f) = rateandstate_vectorized(V2_f, V3_f,  
                        ψ_f, σn, τ2_f, τ3_f, η, RSas_f, RSV0)
    # @show V2_f[1], V3_f[1], ψ_f[1], τ2_f[1], τ3_f[1], RSV0
    for iter = 1:maxiter
        (f_f, g_f, dfx_f, dfy_f, dgx_f, dgy_f) = rateandstate_vectorized(V2_f, V3_f,  
                                        ψ_f, σn, τ2_f, τ3_f, η, RSas_f, RSV0)
        inv_J = map_jacobian_inv.(dfx_f, dfy_f, dgx_f, dgy_f)
        # @show dfx_f[1], dfy_f[1], dgx_f[1], dgy_f[1]
        # @show inv_J[1]
        # println()
        dV2_fV3_f =  -inv_J .* map_cat.(f_f, g_f)
        dV2_f = get_first.(dV2_fV3_f)
        dV3_f = get_second.(dV2_fV3_f)
        V2_f = V2_f .+ α * dV2_f
        V3_f = V3_f .+ α * dV3_f
        # @show dV2_f[1], dV3_f[1]
        
        # TODO vectorized control flow
        if all(abs.(f_f) .< ftol) && all(abs.(dV2_f) .< atolx .+ rtolx .* (abs.(dV2_f) .+ abs.(V2_f))) && all(abs.(g_f) .< ftol) && all(abs.(dV2_f) .< atolx .+ rtolx .* (abs.(dV3_f) .+ abs.(V3_f)))
            return (V2_f, V3_f, f_f, g_f, iter)
        end
    end
    return (V2_f, V3_f, f_f, g_f, -maxiter)
end
    

index = 1:4
f_i, g_i, dfx_i, dfy_i, dgx_i, dgy_i  = rateandstate_vectorized(V2_f[index], V3_f[index], ψ_f[index], σn, τ2[index], τ3[index], η, RSas_f[index], V0)
V2_i_tmp, V3_i_tmp, f_i_tmp, g_i_tmp, maxiter = newtbndv_vectorized(
    rateandstate_vectorized, V2_f[index], V3_f[index], ψ_f[index], σn, τ2_f[index], τ3_f[index], η, RSas_f[index], RSV0; maxiter=100, α=0.2)



# Storing struct in JSON
using JSON
# Define a dictionary
data = Dict(
    "name" => "Alice",
    "age" => 30,
    "city" => "New York",
    "vector" => [1, 2, 3]
)

# write to json file
open("data.json", "w") do file 
    JSON.print(file, data)
end

# load json file
file_contents = read("data.json", String)
json_data = JSON.parse(file_contents)