module ReduceModel

# Load packages
using PowerModels
using SparseArrays
using LinearAlgebra
using PyCall
using Plots
using SpineInterface
using URIParser

# Export utility
export call_rei
export REIOptions
export makePlots
export plot_grid

include("util/types.jl")
include("util/util.jl")
include("util/partition.jl")
include("util/reduced_plots.jl")

include("rei/construct_rei.jl")
include("rei/call_rei.jl")
include("rei/reduce_areas.jl")

const net_layout = PyNULL()

function __init__()
    push!(pyimport("sys")["path"], joinpath(dirname(@__FILE__), "util/"))
    copy!(net_layout, pyimport("network_layout"))
end


end
