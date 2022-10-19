
using Revise,ApproxOperator,DataFrames,XLSX,TimerOutputs

to = TimerOutput()

ndiv = 32
elements, nodes = importmsh("./msh/cantilever_"*string(ndiv)*".msh")
nₚ = length(nodes[:x])
nₑ = length(elements["Ω"])

type = (SNode,:Quadratic2D,:□,:CubicSpline)
s = 2.5*12.0/ndiv*ones(nₚ)
# type = (SNode,:Cubic2D,:□,:CubicSpline)
# s = 3.5*12.0/ndiv*ones(nₚ)
sp = RegularGrid(nodes[:x],nodes[:y],nodes[:z],n = 2,γ = 5)
@timeit to "Total Time" begin
@timeit to "searching" begin
elements["Ω"] = ReproducingKernel{type...,:Tri3}(elements["Ω"],sp)
elements["Ω̃"] = ReproducingKernel{type...,:Tri3}(elements["Ω"])
elements["Γᵗ"] = ReproducingKernel{type...,:Seg2}(elements["Γᵗ"],sp)
elements["Γᵍ"] = ReproducingKernel{type...,:Seg2}(elements["Γᵍ"])
elements["Ω∩Γᵍ"] = elements["Ω"]∩elements["Γᵍ"]
end
@timeit to "prescribling" begin
set𝓖!(elements["Ω"],:TriRK6,:∂1,:∂x,:∂y,:∂z)
set𝓖!(elements["Ω̃"],:TriGI3,:∂1,:∂x,:∂y,:∂z)
set𝓖!(elements["Γᵍ"],:SegRK3,:∂1,:∂x,:∂y,:∂z,:∂̄x,:∂̄y)
set𝓖!(elements["Γᵗ"],:SegRK3,:∂1)
# set𝓖!(elements["Ω"],:TriRK13,:∂1,:∂x,:∂y,:∂z)
# set𝓖!(elements["Ω̃"],:TriGI6,:∂1,:∂x,:∂y,:∂z)
# set𝓖!(elements["Γᵍ"],:SegRK5,:∂1,:∂x,:∂y,:∂z)
# set𝓖!(elements["Γᵗ"],:SegRK5,:∂1)
elements["Γᵍ"] = ReproducingKernel{type...,:Tri3}(elements["Ω"],elements["Γᵍ"])

E = 3E6;ν = 0.3;P = 1000;L = 48;D = 12;
I = D^3/12
EI = E*I

prescribe!(elements["Γᵗ"],:t₂,(x,y,z)->P/2/I*(D^2/4-y^2))
prescribe!(elements["Γᵍ"],:g₁,(x,y,z)->-P*y/6/EI*((6*L-3x)*x + (2+ν)*(y^2-D^2/4)))
prescribe!(elements["Γᵍ"],:g₂,(x,y,z)->P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2))
prescribe!(elements["Γᵍ"],:n₁₁,(x,y,z)->1.0)
prescribe!(elements["Γᵍ"],:n₂₂,(x,y,z)->1.0)
end

push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)
@timeit to "shape functions Ω" set𝝭!(elements["Ω"])
@timeit to "shape functions Ω" set∇̃𝝭!(elements["Ω̃"],elements["Ω"])
@timeit to "shape functions Γᵍ" set∇̃𝝭!(elements["Γᵍ"],elements["Ω∩Γᵍ"])
@timeit to "shape functions Γᵍ" set∇̄𝝭!(elements["Γᵍ"])
@timeit to "shape functions Γᵗ" set𝝭!(elements["Γᵗ"])

coefficient = (:E=>E,:ν=>ν)
ops = [Operator(:∫∫εᵢⱼσᵢⱼdxdy,coefficient...),
       Operator(:∫vᵢtᵢds,coefficient...),
       Operator(:∫σᵢⱼnⱼgᵢds,coefficient...),
       Operator(:∫σ̄ᵢⱼnⱼgᵢds,coefficient...),
       Operator(:Hₑ_PlaneStress,coefficient...)]

k = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)

@timeit to "assembly in Ω" ops[1](elements["Ω̃"],k)
@timeit to "assembly in Γᵗ" ops[2](elements["Γᵗ"],f)
@timeit to "assembly in Γᵍ" ops[3](elements["Γᵍ"],k,f)
@timeit to "assembly in Γᵍ" ops[4](elements["Γᵍ"],k,f)

@timeit to "solve" d = k\f
end
d₁ = d[1:2:2*nₚ]
d₂ = d[2:2:2*nₚ]
push!(nodes,:d₁=>d₁,:d₂=>d₂)
set𝓖!(elements["Ω"],:TriGI16,:∂1,:∂x,:∂y,:∂z)
set∇𝝭!(elements["Ω"])
prescribe!(elements["Ω"],:u,(x,y,z)->-P*y/6/EI*((6*L-3x)*x + (2+ν)*(y^2-D^2/4)))
prescribe!(elements["Ω"],:v,(x,y,z)->P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2))
prescribe!(elements["Ω"],:∂u∂x,(x,y,z)->-P/EI*(L-x)*y)
prescribe!(elements["Ω"],:∂u∂y,(x,y,z)->-P/6/EI*((6*L-3*x)*x + (2+ν)*(3*y^2-D^2/4)))
prescribe!(elements["Ω"],:∂v∂x,(x,y,z)->P/6/EI*((6*L-3*x)*x - 3*ν*y^2 + (4+5*ν)*D^2/4))
prescribe!(elements["Ω"],:∂v∂y,(x,y,z)->P/EI*(L-x)*y*ν)
h1,l2 = ops[5](elements["Ω"])
l2 = log10(l2)
h1 = log10(h1)
h = log10(12.0/ndiv)

if ndiv == 4
Cᵢᵢᵢᵢ = E/(1-ν^2)
Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
Cᵢⱼᵢⱼ = E/2/(1+ν)
inte = 100
n̄ₚ = (inte+1)^2
x = zeros(n̄ₚ)
y = zeros(n̄ₚ)
σ₁₁_ = zeros(n̄ₚ)
σ₂₂_ = zeros(n̄ₚ)
σ₁₂_ = zeros(n̄ₚ)
𝗠 = elements["Ω"][1].𝗠
𝝭 = elements["Ω"][1].𝝭
ap = ReproducingKernel{type...,:Node}([Node(i,nodes) for i in 1:nₚ],Node[],𝗠,𝝭)
for i in 0:inte
    for j in 0:inte
        xᵢ = 48.0/inte*i
        yᵢ = -6.0+12.0/inte*j
        x[(inte+1)*j+i+1] = xᵢ
        y[(inte+1)*j+i+1] = yᵢ
        𝒙 = (xᵢ,yᵢ,0.0)
        uᵢ,ε₁₁,ε₂₂,ε₁₂ = get𝝐(ap,𝒙,sp)
        σ₁₁_[(inte+1)*j+i+1] = Cᵢᵢᵢᵢ*ε₁₁ + Cᵢᵢⱼⱼ*ε₂₂
        σ₂₂_[(inte+1)*j+i+1] = Cᵢᵢⱼⱼ*ε₁₁ + Cᵢᵢᵢᵢ*ε₂₂
        σ₁₂_[(inte+1)*j+i+1] = Cᵢⱼᵢⱼ*ε₁₂
    end
end

df = DataFrame(x=x,y=y,σ₁₁=σ₁₁_,σ₂₂=σ₂₂_,σ₁₂=σ₁₂_)
XLSX.openxlsx("./xlsx/cantilever.xlsx", mode="rw") do xf
    name = "rigsi_hr"
    name∉XLSX.sheetnames(xf) ? XLSX.addsheet!(xf,name) : nothing
    XLSX.writetable!(xf[name],df)
end
end

show(to)
