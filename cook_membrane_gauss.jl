
using  ApproxOperator, LinearAlgebra, Printf
ndiv=5
include("input.jl")
elements, nodes ,elms= import_gauss_quadratic("./msh/cook_membrance_"*string(ndiv)*".msh",:TriGI3)

κ = 400942
μ = 80.1938
E = 9*κ*μ/(3*κ+μ)
# ν = (3*κ-2*μ)/2/(3*κ+μ)
ν =0.499999999
# E = 70.0
# ν = 0.3333
Cᵢᵢᵢᵢ = E*(1-ν)/(1+ν)/(1-2*ν)
Cᵢᵢⱼⱼ = E*ν/(1+ν)/(1-2*ν)
Cᵢⱼᵢⱼ = E/(1+ν)/2

nₚ = length(nodes)
nₑ = length(elements["Ω"])
s = 2.5*44/ndiv*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γᵗ"])
set𝝭!(elements["Γᵍ"])
set∇𝝭!(elements["Ωᶜ"])

prescribe!(elements["Γᵗ"],:t₁=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)

ops = [
    Operator{:Δ∫∫EᵢⱼSᵢⱼdxdy_NeoHookean}(:E=>E,:ν=>ν),
    Operator{:∫∫EᵢⱼSᵢⱼdxdy_NeoHookean}(:E=>E,:ν=>ν),
    Operator{:∫vᵢtᵢds}(),
    Operator{:∫vᵢuᵢds}(:α=>1e7*E),
]

k = zeros(2*nₚ,2*nₚ)
kα = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)
fα = zeros(2*nₚ)
fint = zeros(2*nₚ)
fext = zeros(2*nₚ)
d = zeros(2*nₚ)
Δd= zeros(2*nₚ)
d₁ = zeros(nₚ)
d₂ = zeros(nₚ)

push!(nodes,:d₁=>d₁,:d₂=>d₂)

nmax = 20
P = 0:6.25/nmax:6.25
tolerance=1.0e-10;maxiters=1000;
for (n,p) in enumerate(P)
    if n == 1
        continue
    end
    err_Δd = 1.0
    dnorm = 0.0
    # err_Δf = 1.0
    # fnorm = 0.0
    @printf "Load step=%i,p=%e \n" n p
    fill!(fext,0.0)
    prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->p)
    ops[3](elements["Γᵗ"],fext)
    # fill!(k,0.0)
    # fill!(kα,0.0)
    # fill!(fα,0.0)
    # ops[1](elements["Ω"],k)
    # ops[4](elements["Γᵍ"],kα,fα)
    # k⁻¹ .= inv(k+kα)
    iter = 0
    while err_Δd>tolerance && iter<maxiters
        iter += 1
        fill!(fint,0.0)
        ops[2](elements["Ω"],fint)
        f .= fext-fint

        fill!(k,0.0)
        fill!(kα,0.0)
        fill!(fα,0.0)
        ops[1](elements["Ω"],k)
        ops[4](elements["Γᵍ"],kα,fα)

        # if iter == 1
        #     Δd .= k⁻¹*(f+fα)
        # else
        #     Δd .= k⁻¹*f
        # end

        Δd .= (k+kα)\(f+fα)

        # fnorm = norm(f)
        # fᵗnorm = fnorm+1.0
        # fᵗ .= f
        # λ = 2.0
        # while fᵗnorm ≥ fnorm && λ > tolerance
        #     # println(λ)
        #     fill!(fint,0.0)
        #     λ *= 0.5
        #     d₁ .= d[1:2:2*nₚ]+λ*Δd[1:2:2*nₚ]
        #     d₂ .= d[2:2:2*nₚ]+λ*Δd[2:2:2*nₚ]
        #     ops[2](elements["Ω"],fint)
        #     fᵗ = fext-fint
        #     fᵗnorm = norm(fᵗ)
        #     # println(fnorm)
        #     # println(fᵗnorm)
        # end
        # d .+= λ*Δd 

        d .+= Δd 
        d₁ .= d[1:2:2*nₚ]
        d₂ .= d[2:2:2*nₚ]

        Δdnorm = LinearAlgebra.norm(Δd)
        # Δdnorm = LinearAlgebra.norm(λ*Δd)
        dnorm += Δdnorm
        err_Δd = Δdnorm/dnorm
        # err_Δd = Δdnorm
        # Δfnorm = LinearAlgebra.norm(f+fα)
        # fnorm += Δfnorm
        # err_Δf = Δfnorm/fnorm

        # @printf "iter = %i, err_Δf = %e, err_Δd = %e \n" iter err_Δf err_Δd
        @printf "iter = %i, err_Δd = %e \n" iter err_Δd
    end
end 

u=d₂[3]
println(u)

# fo = open("./vtk/cook_membrance_rkgsi_mix_"*string(ndiv_𝑢)*".vtk","w")
# fo = open("./vtk/cook_membrance_rkgsi_"*string(ndiv_𝑢)*".vtk","w")
fo = open("./vtk/cook_membrance_guass3_"*string(ndiv)*".vtk","w")
@printf fo "# vtk DataFile Version 2.0\n"
# @printf fo "cook_membrance_rkgsi_mix\n"
@printf fo "cook_membrance_guass3\n"
@printf fo "ASCII\n"
@printf fo "DATASET POLYDATA\n"
@printf fo "POINTS %i float\n" nₚ
for p in nodes
    @printf fo "%f %f %f\n" p.x p.y p.z
end
@printf fo "POLYGONS %i %i\n" nₑ 4*nₑ
for ap in elms["Ω"]
    𝓒 = ap.vertices
    @printf fo "%i %i %i %i\n" 3 (x.i-1 for x in 𝓒)...
end
@printf fo "POINT_DATA %i\n" nₚ
@printf fo "VECTORS U float\n"
for p in elements["Ωᶜ"]
    ξ = collect(p.𝓖)[1]
    N = ξ[:𝝭]
    u₁ = 0.0
    u₂ = 0.0
    for (i,x) in enumerate(p.𝓒)
        u₁ += N[i]*x.d₁
        u₂ += N[i]*x.d₂
    end
    @printf fo "%f %f %f\n" u₁ u₂ 0.0
end

@printf fo "TENSORS STRESS float\n"
for p in elements["Ωᶜ"]
    𝓒 = p.𝓒
    𝓖 = p.𝓖
    ε₁₁ = 0.0
    ε₂₂ = 0.0
    ε₁₂ = 0.0

    for (i,ξ) in enumerate(𝓖)
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        for (j,xⱼ) in enumerate(𝓒)
            ε₁₁ += B₁[j]*xⱼ.d₁
            ε₂₂ += B₂[j]*xⱼ.d₂
            ε₁₂ += B₁[j]*xⱼ.d₂ + B₂[j]*xⱼ.d₁
        end
    end
    σ₁₁ = Cᵢᵢᵢᵢ*ε₁₁+Cᵢᵢⱼⱼ*ε₂₂
    σ₂₂ = Cᵢᵢⱼⱼ*ε₁₁+Cᵢᵢᵢᵢ*ε₂₂
    σ₁₂ = Cᵢⱼᵢⱼ*ε₁₂
    @printf fo "%f %f %f\n" σ₁₁ σ₁₂ 0.0
    @printf fo "%f %f %f\n" σ₁₂ σ₂₂ 0.0
    @printf fo "%f %f %f\n" 0.0 0.0 0.0
end
close(fo)
