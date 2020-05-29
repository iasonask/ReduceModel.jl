## Imports

using Revise
using PowerModels
using Ipopt
using Gurobi

## Initializations

# REI monolithic script
println("Calculating REI...")

# REI options
# opf will change them for using functions as arguments :TODO
# Define whether we want to run a power flow (1) or an OPF (2)
pfMethod = 2
# genGroup: 1: group all generators into a big one
# 0: keep existing generators but moved to the common bus
genGroup = 0
# selectPV: 1: select the PV buses based on the net power injections, i.e.
# some PV buses can become PQ
# 0: keep the original PV buses
selectPV = 0

# load network_data
# network_data = PowerModels.parse_file("../examples/case5.m")
# network_data = PowerModels.parse_file("../data/Matpower/case118.m")
# network_data = PowerModels.parse_file("/Users/iasonas/Documents/MATLAB/matpower7.0/data/case118_mod.m")
network_data = PowerModels.parse_file("/Users/iasonas/Documents/MATLAB/matpower7.0/data/case39.m")


# areas
# hardcoded areas for now
# area1 = [1:23..., 25:32..., 113:115..., 117]
# area2 = [33:67...]
# area3 = [24, 68:112..., 116, 118]
area1 = [1:14..., 25, 30, 31, 32, 37, 39]
area2 = [15:24..., 26:29..., 33:36..., 38]
areas = Dict(1 => area1, 2 => area2)
no_areas = length(areas)

# number of buses
no_buses = length(network_data["bus"])
# mpc = ext2int(mpcOri) internal mapping when bus ids is not in order :TODO

baseMVA = network_data["baseMVA"]
bus_full = network_data["bus"]
gen_full = network_data["gen"]
branch_full = network_data["branch"]

# Get the indices of the Ref, PV and PQ buses
indRef = [bus.second["bus_i"] for bus in bus_full if bus.second["bus_type"] == 3]
indPV = [bus.second["bus_i"] for bus in bus_full if bus.second["bus_type"] == 2]
indPQ = [bus.second["bus_i"] for bus in bus_full if bus.second["bus_type"] == 1]

# calculate admittance matrix of full network
Ybus = calc_admittance_matrix(network_data)
conj!(Ybus.matrix) # PowerModels returns the conjugate for some reason
# sum each column to get the diagonal elements
Ydiag = sum(Ybus.matrix, dims=1)

# Save area information in PowerModel object
for area in areas
    for bus in area.second
        bus_full["$bus"]["area"] = area.first
    end
end

## Get voltages and apparent powers from the load flow case

_pf = (pfMethod == 2 ? PowerModels.build_opf : PowerModels.build_pf)
# create PF model
pm = instantiate_model(network_data, ACPPowerModel, _pf)
# run the power flow
results = optimize_model!(pm, optimizer=Ipopt.Optimizer)

# results = ext2int(results) internal mapping when bus ids is not in order :TODO
solution = results["solution"]
(baseMVA, bus_sol, gen_sol) = (solution["baseMVA"], solution["bus"], solution["gen"])
(bus_data, gen_data, branch_data, load_data) = (pm.data["bus"], pm.data["gen"], pm.data["branch"], pm.data["load"])

# pass the results from the dictionary to an array for the next calculations
bus = zeros(Float64, (no_buses, 5))
# assuming ordered buses in 1:no_buses
for bs in solution["bus"]
    id = parse(Int64, bs.first)
    bus[id, 1:VA] = [id, bs.second["vm"], bs.second["va"]]
end

# Get the apparent power injections and voltages at each bus
voltagesRefCase = bus[:, VM] .* exp.(1im * bus[:, VA])
powerRefCase = zeros(Complex, no_buses, 1)
powerRefCaseG = zeros(Complex, no_buses, 1)
powerRefCaseL = zeros(Complex, no_buses, 1)

# Fill in the apparent power (sum of all injections at each bus)
for gn in gen_data
    if gn.second["gen_status"] > 0
        # check if gen is load :TODO
        # g =find(gen(:, GEN_STATUS) > 0 & gen(:, GEN_BUS) == bus(i, BUS_I) & ~isload(gen))
        id = parse(Int64, gn.first)
        powerRefCaseG[id] += gen_sol[gn.first]["pg"] + 1im * gen_sol[gn.first]["qg"]
    end
end

for ld in load_data
    id = parse(Int64, ld.first)
    powerRefCaseL[id] += load_data[ld.first]["pd"] + 1im * load_data[ld.first]["qd"]
    bus[id, PD:QD] = [load_data[ld.first]["pd"], load_data[ld.first]["qd"]]
end

powerRefCase = powerRefCaseG .+ powerRefCaseL

# keep the indices for buses with non-zero injections
indNZ = findall(x -> x > 0, real.(powerRefCase))
indNZ2 = findall(x -> x < 0, real.(powerRefCase))
indRefNZ = intersect(indNZ,indRef)
indPVNZ = setdiff(indNZ, indRefNZ)
indPQNZ = indNZ2 # why is that? :TODO
indPQori = findall(x -> x != 0, bus[:,PD])


print(pm.model)
print(result["solution"])
