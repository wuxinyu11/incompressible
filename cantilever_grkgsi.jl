
using Revise,ApproxOperator
include("input.jl")

fid_𝑢 = "./msh/cantilever_8.msh"
fid_𝑝 = "./msh/cantilever_8.msh"
elements, nodes = import_rkgsi_mix(fid_𝑢,fid_𝑝)

ApproxOperator.cal𝗠!(elements["Ω̄"][1])

elements["Ω̄"][1].𝗚