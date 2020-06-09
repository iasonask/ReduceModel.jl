__precompile__()
module SpineReduce

# Load packages
using JuMP
using PowerModels
using Ipopt
using SparseArrays
using LinearAlgebra
using PyCall

# Export utility
export call_rei

include("util/types.jl")
include("util/util.jl")
include("util/partition.jl")

include("rei/construct_rei.jl")
include("rei/call_rei.jl")
include("rei/reduce_areas.jl")

const nl = PyNULL()

function __init__()
    # push!(pyimport("sys")["path"], joinpath(dirname(@__FILE__), "/util"))
    copy!(nl, pyimport("util/network_layout"))
end


end
