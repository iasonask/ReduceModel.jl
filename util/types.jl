# constants definition for referencing PowerModels data
# buses
const BUS_ID = 1
const BUS_TYPE = 2
const VM = 3
const VA = 4
const PD = 5
const QD = 6

# branches
const L_ID = 1
const F_BUS = 2
const T_BUS = 3

# generators
const GEN_ID = 1
const GEN_BUS = 2


struct PMAreas
    name::Symbol
    no_of_areas::Int64
    clusters::Dict{Int64,Array{Int64,1}}
    data::Dict{Int64,Dict{String,Any}}
    function PMAreas(name, no_of_areas, clusters)
    length(clusters) == no_of_areas || error("Ill-defined clustering of network")
    new(name, no_of_areas, clusters, Dict())
    end
end

struct REIOptions
    pf_model::DataType # ?? AbstractACPModel
    pf_method::Function
end
