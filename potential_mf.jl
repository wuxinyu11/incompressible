using Revise, ApproxOperator, BenchmarkTools
include("input.jl")

elements, nodes = import_rkgsi("./msh/patch_test.msh")

nₚ = length(nodes)
s = 2.5/10*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)
set𝝭!.(elements["Ω"])
set∇𝝭!.(elements["Ω̃"])
set∇𝝭!.(elements["Γᵍ"])

r = 1
u(x,y) = (x+y)^r
∂u∂x(x,y) = r*(x+y)^abs(r-1)
∂u∂y(x,y) = r*(x+y)^abs(r-1)
∂²u∂x²(x,y) = r*(r-1)*(x+y)^abs(r-2)
∂²u∂x∂y(x,y) = r*(r-1)*(x+y)^abs(r-2)
∂²u∂y²(x,y) = r*(r-1)*(x+y)^abs(r-2)
prescribe!(elements["Γᵍ"],:g=>(x,y,z)->u(x,y))
prescribe!(elements["Ω"],:b=>(x,y,z)->-∂²u∂x²(x,y)-∂²u∂y²(x,y))

ops = [
    Operator{:∫∫∇v∇udxdy}(:k=>1.0),
    Operator{:∫vbdΩ}(),
    Operator{:∫vtdΓ}(),
    Operator{:∫∇𝑛vgds}(:k=>1.0,:α=>1e3),
    Operator{:H₁}()
]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)
ops[1](elements["Ω̃"],k)
ops[2](elements["Ω"],f)
# ops[3].(elements["Γᵗ"],f=f)
ops[4](elements["Γᵍ"],k,f)

d = k\f