
using YAML, ApproxOperator,LinearAlgebra 

ndiv = 10
config = YAML.load_file("./yml/bar_z.yml")
path = "./msh/bar_10.msh"
elements, nodes = importmsh(path,config)

nₚ = length(nodes)
s = 2.5/ndiv*ones(nₚ)


push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)
set∇𝝭!(elements["Ω"])
set∇𝝭!(elements["Γᵍ"])


prescribe!(elements["Ω"],:b=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:g=>(x,y,z)->x)
prescribe!([elements["Γᵍ"][1]],:n₁=>(x,y,z)->-1.0)
prescribe!(elements["Γᵍ"][2],:n₁=>(x,y,z)->1.0)


ops = [
    Operator{:∫∇v∇udΩ}(:k=>1.0),
    Operator{:∫vbdΩ}(),
    Operator{:∫vtdΓ}(),
    Operator{:∫vgdΓ}(:α=>1e3),
    Operator{:L₂}()
]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)


ops[1](elements["Ω"],k)
ops[2](elements["Ω"],f)
#  ops[3](elements["Γᵗ"],f)
ops[4](elements["Γᵍ"],k,f)

d = k\f

# push!(nodes,:d=>d)
#     set𝝭!(elements["Ωᴳ"])
#     prescribe!(elements["Ωᴳ"],:u=>(x,y,z)->x)
#     l2 = ops[5](elements["Ωᴳ"])
#     L2=log10(l2)