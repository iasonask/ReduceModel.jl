using SparseArrays


## admittance matrix calculation
# from: https://kersulis.github.io/2015/11/05/admittance-matrix/

"""
    createYbus(f,t,x [,r,b]) -> Y
Create an admittance matrix for AC power flow.
All inputs are real. The output matrix is real if no line
resistances are provided (DC case), and complex otherwise.
* `f`,`t`: vectors encoding all lines (fi,ti)
* `x`: per-unit reactance xi for all lines
* `r`: per-unit resistance ri for all lines
* `b`: per-unit susceptance bi for all lines
"""
function createYbus(
    f::Vector{Int64},
    t::Vector{Int64},
    x::Vector{Float64},
    r = 0.0::Union{Vector{Float64},Float64},
    b = 0.0::Union{Vector{Float64},Float64},
)
    z = r + x * 1im
    y = 1 ./ z
    b = b * 1im
    Y = sparse([f; t; t; f], [t; f; t; f], [-y; -y; y + b ./ 2; y + b ./ 2])

    # for DC power flow, we typically want a matrix with real entries:
    if r == 0
        return imag(Y)
    else
        return Y
    end
end


function ismember(a, b)
    Lia = in.(a, [b])
    Lib = indexin(a, b)
    Lib = [i == nothing ? 0 : i for i in Lib]
    return Lia, Lib
end

function readBrDict(br_i)
    rate_a = get(br_i, "rate_a", 0)
    rate_b = get(br_i, "rate_b", 0)
    rate_c = get(br_i, "rate_c", 0)
    b_ = (br_i["b_fr"] + br_i["b_to"])
    [br_i["f_bus"] br_i["t_bus"] br_i["br_r"] br_i["br_x"] b_ rate_a rate_b rate_c br_i["tap"] br_i["shift"] br_i["br_status"] br_i["angmin"] br_i["angmax"]]
end

function busArrayToDict(bus::Array{Float64,1})
    Dict(
        "zone" => Int64(bus[BUS_AREA]),
        "bus_i" => Int64(bus[BUS_ID]),
        "bus_type" => Int64(bus[BUS_TYPE]),
        "vmax" => bus[VMAX],
        "source_id" => Any["bus", Int64(bus[BUS_ID])],
        "area" => Int64(bus[BUS_AREA]),
        "vmin" => bus[VMIN],
        "index" => Int64(bus[BUS_ID]),
        "va" => bus[VA],
        "vm" => bus[VM],
        "base_kv" => bus[BASE_KV],
    )
end

function genArrayToDict(id::Int64, gen::Array{Float64,1})
    # populate with default values for cost parameters
    Dict(
        "ncost" => 3,
        "qc1max" => 0.0,
        "pg" => gen[PG],
        "model" => 2,
        "shutdown" => 0.0,
        "startup" => 0.0,
        "qc2max" => 0.0,
        "ramp_agc" => 0.0,
        "qg" => gen[QG],
        "gen_bus" => Int64(gen[GEN_BUS]),
        "pmax" => gen[PMAX],
        "ramp_10" => 0.0,
        "vg" => gen[VG],
        "mbase" => gen[MBASE],
        "source_id" => Any["gen", id],
        "pc2" => 0.0,
        "index" => id,
        "cost" => [100.0, 30.0, 0.2],
        "qmax" => gen[QMAX],
        "gen_status" => Int64(gen[GEN_STATUS]),
        "qmin" => gen[QMIN],
        "qc1min" => 0.0,
        "qc2min" => 0.0,
        "pc1" => 0.0,
        "ramp_q" => 0.0,
        "ramp_30" => 0.0,
        "pmin" => gen[PMIN],
        "apf" => 0.0,
    )
end

function genArrayToDict(
    id::Int64,
    gen::Array{Float64,1},
    gen_cost::Dict{String,Any},
)
    # Firt copy all dict and then make necessary changes
    new_gen = deepcopy(gen_cost)
    new_gen["pg"] = gen[PG]
    new_gen["qg"] = gen[QG]
    new_gen["gen_bus"] = Int64(gen[GEN_BUS])
    new_gen["pmax"] = gen[PMAX]
    new_gen["vg"] = gen[VG]
    new_gen["mbase"] = gen[MBASE]
    new_gen["source_id"] = Any["gen", id]
    new_gen["index"] = id
    new_gen["qmax"] = gen[QMAX]
    new_gen["gen_status"] = Int64(gen[GEN_STATUS])
    new_gen["qmin"] = gen[QMIN]
    new_gen["pmin"] = gen[PMIN]
    return new_gen
end

function branchArrayToDict(id::Int64, branch::Array{Float64,1})
    # From PowerModels documenation:
    # All non-transformer branches are given nominal transformer values
    #  (i.e. a tap of 1.0 and a shift of 0.0).
    # All branches have a transformer field indicating if they are a transformer or not.
    # In Matpower a zero TAP indicates that there is no transformer
    br_dict = Dict(
        "br_r" => branch[BR_R],
        "shift" => branch[SHIFT],
        "br_x" => branch[BR_X],
        "g_to" => 0.0, # :TODO check those values
        "g_fr" => 0.0, # :TODO check those values
        "source_id" => Any["branch", id],
        "b_fr" => branch[BR_B] / 2,
        "f_bus" => Int64(branch[F_BUS]),
        "br_status" => Int64(branch[BR_STATUS]),
        "t_bus" => Int64(branch[T_BUS]),
        "b_to" => branch[BR_B] / 2,
        "index" => id,
        "angmin" => branch[ANGMIN],
        "angmax" => branch[ANGMAX],
        "transformer" => (branch[TAP] != 0),
        "tap" => (branch[TAP] == 0 ? 1 : branch[TAP]),
    )

    if branch[RATE_A] > 0
        branch["rate_a"] = branch[RATE_A]
    end
    if branch[RATE_B] > 0
        branch["rate_b"] = branch[RATE_B]
    end
    if branch[RATE_C] > 0
        branch["rate_c"] = branch[RATE_C]
    end
    return br_dict
end

function shuntsToDict(id::Int64, bus_id::Int64, bus::Array{Float64,1})
    Dict(
        "source_id" => Any["bus", bus_id],
        "shunt_bus" => bus_id,
        "status" => 1,
        "gs" => bus[GS],
        "bs" => bus[BS],
        "index" => id,
    )
end

function loadToDict(id::Int64, bus_id::Int64, bus::Array{Float64,1})
    Dict(
        "source_id" => Any["bus", bus_id],
        "load_bus" => bus_id,
        "status" => 1,
        "qd" => bus[QD],
        "pd" => bus[PD],
        "index" => id,
    )
end
