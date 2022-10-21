
using YAML, ApproxOperator

ndiv2 = 12

index = [12,24,48,96]

for n in 1:12

ndiv1 = ndiv2/n

path1 = "./msh/bar_"*string(ndiv1)*".msh"
path2 = "./msh/bar_"*string(ndiv2)*".msh"

config2 = YAML.load_file("./yml/bar1.yml")
elements, nodes = importmsh(path1,path2,config2)

nₚ = length(nodes)

set_memory_𝗠!(elements["Ω̃"],:∇̃)

s = 2.5/ndiv2*ones(nₚ)

push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

set𝝭!(elements["Ω"])
set∇̃𝝭!(elements["Ω̃"],elements["Ω"])
set∇𝝭!(elements["Γᵍ"])

prescribe!(elements["Ω"],:b=>(x,y,z)->-6*x)
prescribe!(elements["Γᵍ"],:g=>(x,y,z)->x^3)
prescribe!([elements["Γᵍ"][1]],:n₁=>(x,y,z)->-1.0)
prescribe!(elements["Γᵍ"][2],:n₁=>(x,y,z)->1.0)

ops = [
    Operator{:∫∇v∇udΩ}(:k=>1.0),
    Operator{:∫vbdΩ}(),
    Operator{:∫vtdΓ}(),
    Operator{:∫∇𝑛vgdΓ}(:k=>1.0),
    Operator{:∫vgdΓ}(:α=>1e3),
    Operator{:L₂}()
]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

ops[1](elements["Ω̃"],k)
ops[2](elements["Ω"],f)
# ops[3](elements["Γᵗ"],f)
ops[4](elements["Γᵍ"],k,f)
ops[5](elements["Γᵍ"],k,f)

 d = k\f

push!(nodes,:d=>d)
set𝓖!(elements["Ω"],:SegGI6,:𝝭)
set𝝭!(elements["Ω"])
prescribe!(elements["Ω"],:u=>(x,y,z)->x^3)
l2 = ops[6](elements["Ω"])
L2 = log10(l2)
logs = log10(ndiv1)

XLSX.openxlsx("./xlsx/bar.xlsx", mode="rw") do xf
    row = "A"
    𝐿₂ = xf[2]
    𝐻₁ = xf[3]
    ind = findfirst(n->n==ndiv2,index)+1
    row = row*string(ind)
    𝐿₂[row] = log10(l2)
    𝐻₁[row] = log10(h1)
end
end