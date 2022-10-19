using YAML, ApproxOperator, XLSX, TimerOutputs

ndiv = 8
config = YAML.load_file("./yml/cantilever_rkgsi_hr.yml")
elements, nodes = importmsh("./msh/cantilever_"*string(ndiv)*".msh", config)
nₚ = length(nodes)
nₑ = length(elements["Ω"])
 s = 3.5*12 / ndiv * ones(nₚ)
 push!(nodes, :s₁ => s, :s₂ => s, :s₃ => s)


 set_memory_𝗠!(elements["Ω̃"],:∇̃)
 set_memory_𝗠!(elements["Γᵍ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∇̃)


 elements["Ω∩Γᵍ"] = elements["Ω"]∩elements["Γᵍ"]

set∇₂𝝭!(elements["Ω"])
set∇̃𝝭!(elements["Ω̃"],elements["Ω"])
set∇̃𝝭!(elements["Γᵍ"],elements["Ω∩Γᵍ"])
set∇₂𝝭!(elements["Γᵍ"])
set𝝭!(elements["Γᵗ"])

E = 3E6;ν = 0.3;P = 1000;L = 48;D = 12;
I = D^3/12
EI = E*I

prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->P/2/I*(D^2/4-y^2))
prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->-P*y/6/EI*((6*L-3x)*x + (2+ν)*(y^2-D^2/4)))
prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2))
prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)

coefficient = (:E=>E,:ν=>ν)
ops = [Operator{:∫∫εᵢⱼσᵢⱼdxdy}(coefficient...),
       Operator{:∫vᵢtᵢds}(coefficient...),
       Operator{:∫σᵢⱼnⱼgᵢds}(coefficient...),
       Operator{:∫σ̄ᵢⱼnⱼgᵢds}(coefficient...),
       Operator{:Hₑ_PlaneStress}(coefficient...)]

k = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)

 ops[1](elements["Ω̃"],k)
 ops[2](elements["Γᵗ"],f)
 ops[3](elements["Γᵍ"],k,f)
 ops[4](elements["Γᵍ"],k,f)

 d = k\f

d₁ = d[1:2:2*nₚ]
d₂ = d[2:2:2*nₚ]
push!(nodes,:d₁=>d₁,:d₂=>d₂)
set𝓖!(elements["Ω"],:TriGI16,:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂𝝭∂z)

# set𝓖!(elements["Ω"],:TriGI16,:∂1,:∂x,:∂y,:∂z)
# set𝝭!(elements["Ω"])
set∇₂𝝭!(elements["Ω"])
prescribe!(elements["Ω"],:u=>(x,y,z)->-P*y/6/EI*((6*L-3x)*x + (2+ν)*(y^2-D^2/4)))
prescribe!(elements["Ω"],:v=>(x,y,z)->P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2))
prescribe!(elements["Ω"],:∂u∂x=>(x,y,z)->-P/EI*(L-x)*y)
prescribe!(elements["Ω"],:∂u∂y=>(x,y,z)->-P/6/EI*((6*L-3*x)*x + (2+ν)*(3*y^2-D^2/4)))
prescribe!(elements["Ω"],:∂v∂x=>(x,y,z)->P/6/EI*((6*L-3*x)*x - 3*ν*y^2 + (4+5*ν)*D^2/4))
prescribe!(elements["Ω"],:∂v∂y=>(x,y,z)->P/EI*(L-x)*y*ν)
h1,l2 = ops[5](elements["Ω"])
# L2 = log10(l2)
# h1 = log10(h1)
# h = log10(12.0/ndiv)
# L2 = log10(l2)
