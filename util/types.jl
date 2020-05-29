

struct PMAreas
    name::Symbol
    no_of_areas::Int64
    clusters::Dict{String,Array{Int64,1}}
    function PMAreas(name, no_of_areas, clusters)
    length(clusters) == no_of_areas || error("Ill-defined clustering of power model")
    new(name, no_of_areas, clusters)
    end
end

struct REIOptions
    pf_model::DataType # ?? AbstractACPModel
    pf_method::Function
end
