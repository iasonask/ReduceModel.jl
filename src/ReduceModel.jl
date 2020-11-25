module ReduceModel

# Load packages
using PowerModels
using SparseArrays
using LinearAlgebra
using PyCall
using Requires

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
    # only load ploting functions if Plots is imported
    @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
        println("Loading Ploting functionality...")
        using Plots
        include("util/reduced_plots.jl")
        export makePlots
        export plot_grid
    end
end

end
