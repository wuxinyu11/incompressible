
using YAML, ApproxOperator,Revise,LinearAlgebra,TOML
file_nodes = "./msh/cantilever_8.msh"
file_elements = "./msh/cantilever_8.msh"

# config = YAML.load_file("./yml/cantilever_rkgsi_nitsche.yml")

config = TOML.parsefile("./toml/cantilever_2D_rkgsi_nitsche_cubic.toml")
elements,nodes = importmsh(file_nodes,config)
nₚ = length(nodes)
 nₑ = length(elements["Ω"])
 s = 3.5*12.0/8*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

set_memory_𝗠!(elements["Ω̃"],:∇̃)
set∇₂𝝭!(elements["Ω"])
set∇̃𝝭!(elements["Ω̃"],elements["Ω"])
set𝝭!(elements["Γᵗ"])
set∇₂𝝭!(elements["Γᵍ"])

E = 3E6;ν = 0.3;P = 1000;L = 48;D = 12;
I = D^3/12
EI = E*I

prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->P/2/I*(D^2/4-y^2))
prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->-P*y/6/EI*((6*L-3x)*x + (2+ν)*(y^2-D^2/4)))
prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2))
prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)

coefficient = (:E=>E,:ν=>ν,:α=>1E3*E)
ops = [Operator{:∫∫εᵢⱼσᵢⱼdxdy}(coefficient...),
              Operator{:∫vᵢtᵢds}(coefficient...),
              Operator{:∫σᵢⱼnⱼgᵢds}(coefficient...),
              Operator{:∫vᵢgᵢds}(coefficient...),
              Operator{:Hₑ_PlaneStress}(coefficient...)]

k = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)

ops[1](elements["Ω̃"],k)
ops[2](elements["Γᵗ"],f)
  ops[3](elements["Γᵍ"],k,f)
ops[4](elements["Γᵍ"],k,f) 

d = k\f
