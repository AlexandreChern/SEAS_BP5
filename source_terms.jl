using Symbolics

@variables x y z
@variables K μ
@variables π_s # = Symbolics.pi

u1 = cos(2*π_s*x) * cos(π_s*y) * cos(π_s*z)
u2 = cos(π_s*x) * cos(3*π_s*y) * cos(π_s*z)
u3 = cos(π_s*x) * cos(π_s*y) * cos(π_s*4*z)

# First derivatives
Dx = Differential(x)
Dy = Differential(y)
Dz = Differential(z)

# Second derivatives
Dxx = Differential(x)^2
Dyy = Differential(y)^2
Dzz = Differential(z)^2

Dxy = Dx * Dy
Dxz = Dx * Dz
Dyz = Dy * Dz

Dzx = Dxz
Dzy = Dyz
Dyx = Dxy

Dx_u1 = expand_derivatives(Dx(u1)) 
Dx_u2 = expand_derivatives(Dy(u1)) 

Dxy_u1 = expand_derivatives(Dxy(u1)) 
Dyx_u1 = expand_derivatives(Dyx(u1)) 
Dyz_u1 = expand_derivatives(Dyz(u1))
Dzx_u1 = expand_derivatives(Dzx(u2))


expand_derivatives(Dxx(u1)) 
expand_derivatives(Dyy(u1))

# substitute(Dxy_u1,Dict(π_s=>π))

# Dx_u1


# Equation for u1

u1_source = expand_derivatives(
        (K - 2//3 * μ) * (Dxx(u1) + Dxy(u2) + Dxz(u3)) 
    +   2 * μ * Dxx(u1)
    +   μ * (Dyy(u1) + Dyx(u2))
    +   μ * (Dzz(u3) + Dzx(u3))
)

substitute(u1_source, Dict(K=>1,μ=>1))

u2_source = expand_derivatives(
        μ * (Dxx(u2) + Dxy(u1))
    +   (K - 2//3 * μ) * (Dyx(u1) + Dyy(u2) + Dyz(u3))
    +   2μ * (Dyy(u2))
    +   μ * (Dzz(u2) + Dzy(u3))
)

substitute(u2_source, Dict(K=>1,μ=>1))

u3_source = expand_derivatives(
        μ * (Dxx(u3) + Dxz(u1))
    +   μ * (Dyy(u3) + Dyz(u2))
    +   (K - 2//3 * μ) 
)