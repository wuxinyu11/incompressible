using Revise, ApproxOperator, TOML, LinearAlgebra

config = TOML.parsefile("./toml/rk6.toml")
elements, nodes = importmsh("./msh/square_rk6_1.msh",config)

dgt = 10
dofs = Tuple{Float64,Float64}[]
for elm in elements["Ω"]
    𝓖 = elm.𝓖
    for ξ in 𝓖
        x = trunc(ξ.x; digits=dgt)
        y = trunc(ξ.y; digits=dgt)
        push!(dofs,(x,y))
    end
end
unique!(dofs)


nᵢ = length(dofs)
k = zeros(nᵢ,nᵢ)
g₁ᵢ = zeros(6,3)
g₂ᵢ = zeros(6,3)
g₁ⱼ = zeros(6,3)
g₂ⱼ = zeros(6,3)
G⁻¹ = [9. -12. -12.;-12. 24. 12.;-12. 12. 24.]
get𝒒(ξ,η) = (1.0,ξ,η)
get∂𝒒∂ξ(ξ,η) = (0.0,1.0,0.0)
get∂𝒒∂η(ξ,η) = (0.0,0.0,1.0)

for elm in elements["Ω"][1:2]
    𝓖 = elm.𝓖
    𝐴 = ApproxOperator.get𝐴(elm)
    x₁ = elm.𝓒[1].x
    x₂ = elm.𝓒[2].x
    x₃ = elm.𝓒[3].x
    y₁ = elm.𝓒[1].y
    y₂ = elm.𝓒[2].y
    y₃ = elm.𝓒[3].y
    ∂ξ∂x = y₃-y₂
    ∂ξ∂y = x₂-x₃
    ∂η∂x = y₁-y₃
    ∂η∂y = x₃-x₁
    for (i,ξᵢ) in enumerate(𝓖)
        ξ = ξᵢ.ξ
        η = ξᵢ.η
        q = get𝒒(ξ,η)
        ∂q∂ξ = get∂𝒒∂ξ(ξ,η)
        ∂q∂η = get∂𝒒∂η(ξ,η)
        wb = ξᵢ.wᵇ
        wi = ξᵢ.w
        D₁ = ξᵢ.D₁
        D₂ = ξᵢ.D₂
        for (k,qₖ) in enumerate(q)
            ∂qₖ∂ξ = ∂q∂ξ[k]
            ∂qₖ∂η = ∂q∂η[k]
            g₁ᵢ[i,k] = wi/2.0*(∂qₖ∂ξ*∂ξ∂x+∂qₖ∂η*∂η∂x)+wb*qₖ*D₁
            g₂ᵢ[i,k] = wi/2.0*(∂qₖ∂ξ*∂ξ∂y+∂qₖ∂η*∂η∂y)+wb*qₖ*D₂
        end
    end
    for (j,ξⱼ) in enumerate(𝓖)
        ξ = ξⱼ.ξ
        η = ξⱼ.η
        q = get𝒒(ξ,η)
        ∂q∂ξ = get∂𝒒∂ξ(ξ,η)
        ∂q∂η = get∂𝒒∂η(ξ,η)
        wb = ξⱼ.wᵇ
        wi = ξⱼ.w
        D₁ = ξⱼ.D₁
        D₂ = ξⱼ.D₂
        for (k,qₖ) in enumerate(q)
            ∂qₖ∂ξ = ∂q∂ξ[k]
            ∂qₖ∂η = ∂q∂η[k]
            g₁ⱼ[j,k] = wi/2.0*(∂qₖ∂ξ*∂ξ∂x+∂qₖ∂η*∂η∂x)+wb*qₖ*D₁
            g₂ⱼ[j,k] = wi/2.0*(∂qₖ∂ξ*∂ξ∂y+∂qₖ∂η*∂η∂y)+wb*qₖ*D₂
        end
    end
    for (i,ξᵢ) in enumerate(𝓖)
        x = trunc(ξᵢ.x; digits=dgt)
        y = trunc(ξᵢ.y; digits=dgt)
        I = findfirst(x_->x_==(x,y),dofs)
        for (j,ξⱼ) in enumerate(𝓖)
            x = trunc(ξⱼ.x; digits=dgt)
            y = trunc(ξⱼ.y; digits=dgt)
            J = findfirst(x_->x_==(x,y),dofs)
            for ii in 1:3
                for jj in 1:3
                    k[I,J] += g₁ᵢ[i,ii]*G⁻¹[ii,jj]*g₁ⱼ[j,jj]/𝐴 + g₂ᵢ[i,ii]*G⁻¹[ii,jj]*g₂ⱼ[j,jj]/𝐴
                    # k[I,J] += g₁ᵢ[i,ii]*G⁻¹[ii,jj]*g₁ⱼ[j,jj]/𝐴
                end
            end
        end
    end
end

rk = rank(k)
