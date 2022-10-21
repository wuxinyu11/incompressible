
using YAML, ApproxOperator

ndiv = 10
config1 = YAML.load_file("./yml/bar.yml")
path = "./msh/bar_10.msh"

elements1, nodes1 = importmsh(path,config1)
config2 = YAML.load_file("./yml/bar1.yml")
elements2, nodes2 = importmsh(path,path,config2)

nₚ₁ = length(nodes1)
nₚ₂ = length(nodes2)

set_memory_𝗠!(elements1["Ω̃"],:∇̃)
set_memory_𝗠!(elements2["Ω̃"],:∇̃)

s = 2.5/ndiv*ones(nₚ₁)
s = 2.5/ndiv*ones(nₚ₂)

push!(nodes1,:s₁=>s,:s₂=>s,:s₃=>s)
push!(nodes2,:s₁=>s,:s₂=>s,:s₃=>s)

set𝝭!(elements1["Ω"])
set∇̃𝝭!(elements1["Ω̃"],elements1["Ω"])
set∇𝝭!(elements1["Γᵍ"])

set𝝭!(elements2["Ω"])
set∇̃𝝭!(elements2["Ω̃"],elements2["Ω"])
set∇𝝭!(elements2["Γᵍ"])

prescribe!(elements1["Ω"],:b=>(x,y,z)->1.0)
prescribe!(elements1["Γᵍ"],:g=>(x,y,z)->x)
prescribe!([elements1["Γᵍ"][1]],:n₁=>(x,y,z)->-1.0)
prescribe!(elements1["Γᵍ"][2],:n₁=>(x,y,z)->1.0)

prescribe!(elements2["Ω"],:b=>(x,y,z)->1.0)
prescribe!(elements2["Γᵍ"],:g=>(x,y,z)->x)
prescribe!([elements2["Γᵍ"][1]],:n₁=>(x,y,z)->-1.0)
prescribe!(elements2["Γᵍ"][2],:n₁=>(x,y,z)->1.0)

ops = [
    Operator{:∫∇v∇udΩ}(:k=>1.0),
    Operator{:∫vbdΩ}(),
    Operator{:∫vtdΓ}(),
    Operator{:∫∇𝑛vgdΓ}(:k=>1.0),
    Operator{:∫vgdΓ}(:α=>1e3)
]

k1 = zeros(nₚ₁,nₚ₁)
f1 = zeros(nₚ₁)
k2 = zeros(nₚ₂,nₚ₂)
f2= zeros(nₚ₂)
ops[1](elements1["Ω̃"],k1)
# ops[2](elements1["Ω"],f1)
# ops[3](elements["Γᵗ"],f)
# ops[4](elements1["Γᵍ"],k1,f1)
# ops[5](elements1["Γᵍ"],k1,f1)

ops[1](elements2["Ω̃"],k2)
# ops[2](elements2["Ω"],f2)
# ops[3](elements["Γᵗ"],f)
# ops[4](elements2["Γᵍ"],k2,f2)
# ops[5](elements2["Γᵍ"],k2,f2)

# d1 = k1\f1
# d2 = k2\f2
F=f1-f2
K=k1-k2
# D=d1-d2
