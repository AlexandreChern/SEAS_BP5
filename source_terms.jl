using Symbolics

@syms x_s y_s z_s
@variables K_s μ_s
@variables π_s # = Symbolics.pi

u1_s = cos(2 * π_s * x_s) * cos(π_s * y_s) * cos(π_s * z_s)
u2_s = cos(π_s * x_s) * cos(3 * π_s * y_s) * cos(π_s * z_s)
u3_s = cos(π_s * x_s) * cos(π_s * y_s) * cos(π_s * 4 * z_s)

# First derivatives
Dx = Differential(x_s)
Dy = Differential(y_s)
Dz = Differential(z_s)

# Second derivatives
Dxx = Differential(x_s)^2
Dyy = Differential(y_s)^2
Dzz = Differential(z_s)^2

Dxy = Dx * Dy
Dxz = Dx * Dz
Dyz = Dy * Dz

Dzx = Dxz
Dzy = Dyz
Dyx = Dxy

Dx_u1 = expand_derivatives(Dx(u1_s)) 
Dx_u2 = expand_derivatives(Dy(u1_s)) 

Dxy_u1 = expand_derivatives(Dxy(u1_s)) 
Dyx_u1 = expand_derivatives(Dyx(u1_s)) 
Dyz_u1 = expand_derivatives(Dyz(u1_s))
Dzx_u1 = expand_derivatives(Dzx(u2_s))


expand_derivatives(Dxx(u1_s)) 
expand_derivatives(Dyy(u1_s))

# substitute(Dxy_u1,Dict(π_s=>π))

# Dx_u1


# Equation for u1

u1_source = expand_derivatives(
        (K_s - 2//3 * μ_s) * (Dxx(u1_s) + Dxy(u2_s) + Dxz(u3_s)) 
    +   2 * μ_s * Dxx(u1_s)
    +   μ * (Dyy(u1) + Dyx(u2))
    +   μ * (Dzz(u3) + Dzx(u3))
)

substitute(u1_source, Dict(K_s=>1,μ_s=>1)) # how to substitute symbolics with real values

u2_source = expand_derivatives(
        μ * (Dxx(u2_s) + Dxy(u1_s))
    +   (K - 2//3 * μ) * (Dyx(u1_s) + Dyy(u2_s) + Dyz(u3_s))
    +   2 * μ * (Dyy(u2_s))
    +   μ * (Dzz(u2_s) + Dzy(u3_s))
)

substitute(u2_source, Dict(K=>1,μ=>1))

u3_source = expand_derivatives(
        μ * (Dxx(u3_s) + Dxz(u1_s))
    +   μ * (Dyy(u3_s) + Dyz(u2_s))
    +   (K - 2//3 * μ) * (Dzx(u1_s) + Dzy(u2_s) + Dzz(u3_s))
    +   2 * μ * (Dzz(u3_s))
)


# substitute((K_s - 2//3 * μ_s) * (Dxx(u1_s) + Dxy(u2_s) + Dxz(u3_s)),Dict(K_s=>1, μ_s=>1,π_s => π, x_s=>x))


# Testing symbolics vector
@syms x_test[]
substitute(cos.(x_test), Dict(x_test=>[0,π/2,π])) # this is working


@syms x_s_new[] 
Dx_broadcast = Differential.(x_s_new)
Dxx_expand = expand_derivatives(Dx_broadcast(cos.(π_s .* x_s_new)))

substitute(Dxx_expand, Dict(π_s=>π,x_s_new=>[0,1]))