using ApproxOperator, GLMakie, CairoMakie, XLSX, YAML

include("input.jl")
elements, nodes = import_gauss_quadratic("./msh/cook_membrance_25.msh",:TriGI3)

nκ = 400942
μ = 80.1938
E = 9*κ*μ/(3*κ+μ)
ν = (3*κ-2*μ)/2/(3*κ+μ)
# E = 70.0
# ν = 0.3333

nₚ = length(nodes)
nₑ = length(elements["Ω"])
s = 2.5*44/25*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)


data = getfield(nodes[1],:data)
sp = ApproxOperator.RegularGrid(data[:x][2],data[:y][2],data[:z][2],n=3,γ=5)

sheet = XLSX.readtable("./xlsx/triangular_contour.xlsx", "Sheet1")

inte = 300
x = [-10/3+10/inte*i for i in 0:inte]
y = [-10/3^0.5+20/3^0.5/inte*i for i in 0:inte]


D = 1.0
ν = 0.3
# w₁₁(x,y) = -1/4320*(10+3*x)*(800-420*x+45*x^2-27*y^2)
# w₂₂(x,y) = -1/4320*(10+3*x)*(800-60*x-9*x^2-81*y^2)
# w₁₂(x,y) = -1/4320*18*y*(200-60*x-9*x^2-9*y^2)

w₁₁(x,y) = 1/640*(6*x-20)*(4/9*100-x^2-y^2)+1/640*(3*x^2-3*y^2-20*x)*(-2*x)*2-2/640*(x^3-3*y^2*x-10(x^2+y^2)+4/27*1000)
w₂₂(x,y) = 1/640*(-6*x-20)*(4/9*100-x^2-y^2)+1/640*(0-6*y*x-20*y)*(-2*y)*2-2/640*(x^3-3*y^2*x-10(x^2+y^2)+4/27*1000)
w₁₂(x,y) = 1/640*(-6*y)*(4/9*100-x^2-y^2)+1/640*(3*x^2-3*y^2-20*x)*(-2*y)+1/640*(0-6*y*x-20*y)*(-2*x)
M₁₁(x,y) = - D*(w₁₁(x,y)+ν*w₂₂(x,y))
M₂₂(x,y) = - D*(ν*w₁₁(x,y)+w₂₂(x,y))
M₁₂(x,y) = - D*(1-ν)*w₁₂(x,y)

Mᵣᵣ = Array{Union{Missing,Float64}}(missing,inte+1,inte+1)
for j in 1:inte+1
    for i in 1:inte+1
        xᵢ = x[i]
        yᵢ = y[j]
        if yᵢ<=-3^0.5/3*xᵢ+20*3^0.5/9 && yᵢ>=3^0.5/3*xᵢ-20*3^0.5/9
        # if yᵢ<=-3^0.5/3*xᵢ+20*3^0.5/9
            # Mᵣᵣ[i,j] = M₁₁(xᵢ,yᵢ)
            Mᵣᵣ[i,j] = M₁₂(xᵢ,yᵢ)
            # Mᵣᵣ[i,j] = M₂₂(xᵢ,yᵢ)
        end
    end
end

f = Figure(fontsize=24,fonts = (; regular = "Times New Roman"))
ax = Axis(f[1, 1])
ax.aspect = 1
hidespines!(ax)
hidedecorations!(ax)
surface!(x,y,Mᵣᵣ,
    # colormap = :balance,
    # colormap = :bluesreds,
    colormap = :haline,
    shading = false,
    colorrange=(-1.,1.)
    )
contour!(x,y,Mᵣᵣ,color=:whitesmoke,levels=-0.9:0.2:1.1)
# contour!(x,y,Mᵣᵣ,color=:whitesmoke)
lines!([-10/3.,20/3,-10/3,-10/3],[-10/3^0.5,0.0,10/3^0.5,-10/3^0.5],color=:black)
Colorbar(f[1,2], colormap=:haline,limits = (-1, 1),ticks = -1.0:0.5:1.0, size=50)
f
# save("./figure/triangular_exact.png",f)

# Mᵣᵣ = Array{Union{Missing,Float64}}(missing,inte+1,inte+1)
# n = 5
# d = sheet[1][n]
# for j in 1:inte+1
#     for i in 1:inte+1
#         xᵢ = x[i]
#         yᵢ = y[j]
#         if yᵢ<=-3^0.5/3*xᵢ+20*3^0.5/9 && yᵢ>=3^0.5/3*xᵢ-20*3^0.5/9
#             Mᵣᵣ[i,j] = 0.0
#             indices = sp(xᵢ,yᵢ,0.0)
#             ξ = ApproxOperator.SNode((1,1,0),Dict([:x=>(2,[xᵢ]),:y=>(2,[yᵢ]),:z=>(2,[0.])]))
#             ap = ApproxOperator.ReproducingKernel{:Cubic2D,:□,:QuinticSpline,:Tri3}([nodes[ind] for ind in indices],[ξ],Dict{Symbol,ApproxOperator.SymMat}())
#             set_memory_𝗠!(ap,:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂²𝝭∂x²,:∂²𝝭∂x∂y,:∂²𝝭∂y²)
#             set_memory_𝝭!([ap],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂²𝝭∂x²,:∂²𝝭∂x∂y,:∂²𝝭∂y²)
#             set∇²₂𝝭!(ap,ξ)
#             B₁₁ = ξ[:∂²𝝭∂x²]
#             B₁₂ = ξ[:∂²𝝭∂x∂y]
#             B₂₂ = ξ[:∂²𝝭∂y²]
#             for (ind,I) in enumerate(indices)
#                 M₁₁ = -D*(B₁₁[ind]+ν*B₂₂[ind])*d[I]
#                 M₂₂ = -D*(ν*B₁₁[ind]+B₂₂[ind])*d[I]
#                 M₁₂ = -D*(1-ν)*B₁₂[ind]*d[I]
#                 Mᵣᵣ[i,j] += M₁₂
#             end
#         end
#     end
# end
# f = Figure()
# ax = Axis(f[1, 1])
# ax.aspect = 1
# hidespines!(ax)
# hidedecorations!(ax)
# surface!(x,y,Mᵣᵣ,colormap = :haline, shading = false, colorrange=(-1.,1.))
# contour!(x,y,Mᵣᵣ,color=:whitesmoke,levels=-0.9:0.2:1.1)
# lines!([-10/3.,20/3,-10/3,-10/3],[-10/3^0.5,0.0,10/3^0.5,-10/3^0.5],color=:black)
# save("./figure/triangular_"*string(n)*".png",f)