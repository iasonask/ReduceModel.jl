# constants definition for referencing PowerModels data
# buses
const BUS_ID = 1
const BUS_TYPE = 2
const PD = 3
const QD = 4
const GS = 5
const BS = 6
const BUS_AREA = 7
const VM = 8
const VA = 9
const BASE_KV = 10 # check numbering
const ZONE = 11 # check numbering
const VMAX = 12
const VMIN = 13
const BUS_ARRAY_SIZE = 13

# branches
const F_BUS = 1
const T_BUS = 2
const BR_R = 3
const BR_X = 4
const BR_B = 5
const RATE_A = 6
const RATE_B = 7
const RATE_C = 8
const TAP = 9
const SHIFT = 10
const BR_STATUS = 11
const ANGMIN = 12
const ANGMAX = 13
const BR_ARRAY_SIZE = 13

# generators
const GEN_BUS = 1
const PG = 2
const QG = 3
const QMAX = 4
const QMIN = 5
const VG = 6
const MBASE = 7
const GEN_STATUS = 8
const PMAX = 9
const PMIN = 10
const GEN_ARRAY_SIZE = 10


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
    pf_model::DataType
    pf_method::Function
end
