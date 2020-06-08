module SpineReduce

# Load packages
using JuMP
using PowerModels
using Ipopt
using SparseArrays
using LinearAlgebra

# Export utility
export call_rei

include("util/types.jl"))
include("util/util.jl"))
include("util/partition.jl"))

include("rei/construct_rei.jl"))
include("rei/call_rei.jl"))
include("rei/reduce_areas.jl"))

end
