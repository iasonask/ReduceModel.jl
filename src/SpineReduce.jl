# __precompile__()
module SpineReduce

# Load packages
using Revise
using JuMP
using PowerModels
using Ipopt
using SparseArrays
using LinearAlgebra
using PyCall

# Export utility
export call_rei
export REIOptions

include("util/types.jl")
include("util/util.jl")
include("util/partition.jl")

include("rei/construct_rei.jl")
include("rei/call_rei.jl")
include("rei/reduce_areas.jl")

const net_layout = PyNULL()

function __init__()
    push!(pyimport("sys")["path"], joinpath(dirname(@__FILE__), "util/"))
    copy!(net_layout, pyimport("network_layout"))
end


end
