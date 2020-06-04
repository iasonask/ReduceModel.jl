## Imports

using Revise
using PowerModels
using Ipopt, Gurobi
using SparseArrays, LinearAlgebra

include(joinpath(dirname(@__FILE__), "..", "util/types.jl"))
include(joinpath(dirname(@__FILE__), "..", "util/util.jl"))

## Initializations

# REI monolithic script
println("Calculating REI...")

# REI options
# opf will change them for using functions as arguments :TODO
# Define whether we want to run a power flow (1) or an OPF (2)
pfMethod = 2
# genGroup: true: group all generators into a big one
# false: keep existing generators but moved to the common bus
genGroup = false
# selectPV: true: select the PV buses based on the net power injections, i.e.
# some PV buses can become PQ
# false: keep the original PV buses
selectPV = false

# load network_data
# network_data = PowerModels.parse_file("../examples/case5.m")
# network_data = PowerModels.parse_file("../data/Matpower/case118.m")
# network_data = PowerModels.parse_file("/Users/iasonas/Documents/MATLAB/matpower7.0/data/case118_mod.m")
network_data = PowerModels.parse_file("/Users/iasonas/Documents/MATLAB/matpower7.0/data/case39.m")
caseName = "case39.m"

# Set power flow model
PFModel = ACPPowerModel

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
no_branches = length(network_data["branch"])
no_gen = length(network_data["gen"])
# mpc = ext2int(mpcOri) internal mapping when bus ids is not in order :TODO

baseMVA = network_data["baseMVA"]
bus_full = network_data["bus"]
gen_full = network_data["gen"]
branch_full = network_data["branch"]

# Get the indices of the Ref, PV and PQ buses
indRef = [bus.second["bus_i"] for bus in bus_full if bus.second["bus_type"] == 3]
indPV = sort([bus.second["bus_i"] for bus in bus_full if bus.second["bus_type"] == 2])
indPQ = sort([bus.second["bus_i"] for bus in bus_full if bus.second["bus_type"] == 1])

# calculate admittance matrix of full network
Ybus = calc_admittance_matrix(network_data)
Yadm = Ybus.matrix
conj!(Yadm) # PowerModels returns the conjugate for some reason
# sum each column to get the diagonal elements
Ydiag = sum(Yadm, dims=1)

# Save area information in PowerModel object
for area in areas
    for bus in area.second
        bus_full["$bus"]["area"] = area.first
    end
end

## Get voltages and apparent powers from the load flow case

_pf = (pfMethod == 2 ? PowerModels.build_opf : PowerModels.build_pf)
# create PF model
pm = instantiate_model(network_data, PFModel, _pf)
# run the power flow
results = optimize_model!(pm, optimizer=Ipopt.Optimizer)

# results = ext2int(results) internal mapping when bus ids are not in order :TODO
solution = results["solution"]
(baseMVA, bus_sol, gen_sol) = (solution["baseMVA"], solution["bus"], solution["gen"])
(bus_data, gen_data, branch_data, load_data) = (pm.data["bus"], pm.data["gen"], pm.data["branch"], pm.data["load"])

# pass the results from the dictionary to arrays for easier calculations
bus = zeros(Float64, (no_buses, 6))
# sorted buses in 1:no_buses
for bs in bus_sol
    id = parse(Int64, bs.first)
    bus[id, BUS_ID:VA] = [id, bus_data[bs.first]["bus_type"], bs.second["vm"], bs.second["va"]]
end

# create branch table for simplifying operations
branch = zeros(Int64, (no_branches, 2))
for ln in branch_data
    id = parse(Int64, ln.first)
    branch[id, F_BUS:T_BUS] = [ln.second["f_bus"], ln.second["t_bus"]]
end

# generators
gen = zeros(Float64, (no_gen, 10))
# sorted buses in 1:no_buses
for gn in gen_data
    id = parse(Int64, gn.first)
    gen[id, GEN_BUS:QG] = [gn.second["gen_bus"], gn.second["pg"], gn.second["qg"]]
    gen[id, QMAX:VG] = [gn.second["qmax"], gn.second["qmin"], gn.second["vg"]]
    gen[id, MBASE:PMIN] = [gn.second["mbase"], gn.second["gen_status"], gn.second["pmax"], gn.second["pmin"]]
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
        id = gn.second["gen_bus"]
        powerRefCaseG[id] += gen_sol[gn.first]["pg"] + 1im * gen_sol[gn.first]["qg"]
    end
end

# Fill in load
for ld in load_data
    id = ld.second["load_bus"]
    powerRefCaseL[id] = -(ld.second["pd"] + 1im * ld.second["qd"])  # :TODO should it -= here instead?
    bus[id, PD:QD] = [ld.second["pd"], ld.second["qd"]]
end

powerRefCase = powerRefCaseG .+ powerRefCaseL

# keep the indices for buses with non-zero injections
indNZ = findall(x -> x > 0, real.(powerRefCase))
indPQNZ = findall(x -> x < 0, real.(powerRefCase))
indRefNZ = intersect(indNZ, indRef)
indPVNZ = setdiff(indNZ, indRefNZ)
indPQori = findall(x -> x != 0, bus[:,PD])

## Construct the new area matrices
# We then use the information stored in areaNames, to identify the
# generator and load buses in each area.
# areaInfo is a structure containing information for the REI of each area
areaInfo = PMAreas(:cluster, 2, areas)

# We will also need to keep information about interarea lines.
nl = no_branches
lineInfos = zeros(nl,10)
let
    # The size of the new system that we update during reducing
    nNew = 0

    # The total number of generators
    nGenNew = 0

    for i in 1:areaInfo.no_of_areas
        # Extract the names of the buses in area i
        indAreaI = sort(areaInfo.clusters[i])

        # Identify the frontier buses
        linesFr = in.(branch[:, F_BUS], [indAreaI])
        linesTo = in.(branch[:, T_BUS], [indAreaI])
        bordLines = xor.(linesFr,linesTo) # logical indexing of the border lines
        bordLineFr = bordLines .& linesFr # indices of the border lines from area i
        bordLineTo = bordLines .& linesTo # indices of the border lines to area i
        indBusesFr = sort(unique(branch[bordLineFr, F_BUS])) # indices of the border buses (direction from)
        indBusesTo = sort(unique(branch[bordLineTo, T_BUS])) # indices of the border buses (direction to)
        indBusBord = sort(union(indBusesFr, indBusesTo))

        # # Making sure that indBusBord is a column of indices :TODO check if necessary
        # if size(indBusBord,1) == 1
        #     indBusBord = indBusBord';
        # end

        # Essential buses: border buses and the two new REI buses
        (_, bordInArea) = ismember(indBusBord, indAreaI)
        sort!(bordInArea)

        # log_ind = in.(a, [b])
        # inds = [1:length(log_ind)...][log_ind]

        # the dict areaInfo is used to store information about the area.
        areaInfoI = Dict()
        areaInfoI["buses"] = indAreaI
        areaInfoI["indE"] = indAreaI[bordInArea] # storing the essential border indices (original indices)
        areaInfoI["indNE"] = setdiff(indAreaI, indAreaI[bordInArea]) # storing non essential buses (original indices)

        # area information: area from which or to which a inter area line goes
        lineInfos[bordLineFr, 1] .= i
        lineInfos[bordLineTo, 2] .= i

        # Identify the PV and PQ buses
        if selectPV # select based on net power injection
            refYN = in.(indAreaI, [indRefNZ]) .& .~in.(indAreaI, [indBusBord])
            pvYN =  in.(indAreaI, [indPVNZ]) .& .~in.(indAreaI, [indBusBord])
            pqYN = in.(indAreaI, [indPQNZ]) .& .~in.(indAreaI, [indBusBord])
            allPV_log = (in.(indAreaI, [indPVNZ]) .| in.(indAreaI, [indRefNZ]))
        else # select the original PV buses
            refYN = in.(indAreaI, [indRef]) .& .~in.(indAreaI, [indBusBord])
            pvYN =  in.(indAreaI, [indPV]) .& .~in.(indAreaI, [indBusBord])
            pqYN = in.(indAreaI, [indPQori]) .& .~in.(indAreaI, [indBusBord])
            allPV_log = (in.(indAreaI, [indPV]) .| in.(indAreaI, [indRef]))
        end

        # original PV indices (not based on the net real power injections)
        oriPV = in.(indAreaI, [indPV]) .| in.(indAreaI, [indRef])


        areaIPVbuses = sort(union(indAreaI[pvYN], indAreaI[refYN]))
        areaIRefbuses = indAreaI[refYN]
        areaIPQbuses = indAreaI[pqYN]
        areaIoriPV = indAreaI[oriPV]
        allPV = indAreaI[allPV_log]

        """
        The buses that were PV in the original system but are not based on
        the net real power injection criteria will be counted as PQ. We keep
        a list of such buses to remove the generators from such buses in the
        end when reducing the model.
        """
        PVtoPQbuses = setdiff(areaIoriPV, allPV)

        # Storing the EXTERNAL = ORIGINAL indices
        areaInfoI["PVtoPQbuses"] = PVtoPQbuses
        areaInfoI["PVbuses"] = areaIPVbuses
        areaInfoI["Refbuses"] = areaIRefbuses
        areaInfoI["PQbuses"] = areaIPQbuses

        # check whether one of the aggregated generators was the slack
        # not done!!! :TODO what does he mean here?
        isRefAgg = in.(indRefNZ, [areaIPVbuses])
        if ~isempty(isRefAgg)
            warn("One of the aggregated generators was the slack bus!")
        end

        (areaIBordRef_log, areaIBordRef_int) = ismember(indRef, indAreaI[bordInArea])
        areaIBordRef_int = filter(x -> x != 0, areaIBordRef_int)
        (areaIBordPV_log, areaIBordPV_int) = ismember(indPV, indAreaI[bordInArea])
        areaIBordPV_int = filter(x -> x != 0, areaIBordPV_int)
        areaIBordRef = indRef[areaIBordRef_log]
        areaIBordPV = indPV[areaIBordPV_log]

        # Definition of the size of the full REI matrix
        sizeAreaI = length(indAreaI)
        totalSizeREI = sizeAreaI + 2 * (~isempty(areaIPVbuses)) + 2 * (~isempty(areaIPQbuses))

        # Warning: in the next lines, the numbering is intra area, i.e., the buses
        # are identified by a number between 1 and the number of buses in this area.
        newPD = []
        newQD = []
        newPG = []
        newQG = []
        # variable to store whether the new REI buses exist in this area
        existPV = 1
        existPQ = 1
        if ((~isempty(areaIPVbuses)) && (~isempty(areaIPQbuses)))
            indE = [bordInArea; sizeAreaI + 3; sizeAreaI + 4]
            newTypes = [maximum(bus[areaIPVbuses, BUS_TYPE]); 1] # types of the new REI buses

            if selectPV
                newPD = [0; -sum(real(powerRefCase[areaIPQbuses]))] * baseMVA
                # the minus sign is because we want the value of the consumption
                newQD = [0; -sum(imag(powerRefCase[areaIPQbuses]))] * baseMVA
                newPG = sum(real(powerRefCase[areaIPVbuses])) * baseMVA
                newQG = sum(imag(powerRefCase[areaIPVbuses])) * baseMVA
            else
                newPD = [0; -sum(real(powerRefCaseL[areaIPQbuses]))] * baseMVA
                # the minus sign is because we want the value of the consumption
                newQD = [0; -sum(imag(powerRefCaseL[areaIPQbuses]))] * baseMVA
                newPG = sum(real(powerRefCaseG[areaIPVbuses])) * baseMVA
                newQG = sum(imag(powerRefCaseG[areaIPVbuses])) * baseMVA
            end

            areaInfoI["shiftGen"] = -2 # indicating that the new REI gen bus is the one before the last
            areaInfoI["shiftLoad"] = -1

        elseif (~isempty(areaIPVbuses))
            indE = [bordInArea; sizeAreaI+2]
            newTypes = maximum(bus[areaIPVbuses, BUS_TYPE])
            if selectPV
                newPG = sum(real(powerRefCase[areaIPVbuses])) * baseMVA
                newQG = sum(imag(powerRefCase[areaIPVbuses])) * baseMVA
            else
                newPG = sum(real(powerRefCaseG[areaIPVbuses])) * baseMVA
                newQG = sum(imag(powerRefCaseG[areaIPVbuses])) * baseMVA
            end
            areaInfoI["shiftGen"] = -1 # indicating that the new REI gen bus is the one before the last
            areaInfoI["shiftLoad"] = 0 # indicating that the new REI load bus does not exist
            existPQ = 0

        elseif (~isempty(areaIPQbuses))
            indE = [bordInArea; sizeAreaI+2]
            newTypes = 1
            if selectPV
                newPD = -sum(real(powerRefCase[areaIPQbuses])) * baseMVA
                newQD = -sum(imag(powerRefCase[areaIPQbuses])) * baseMVA
            else
                newPD = -sum(real(powerRefCaseL[areaIPQbuses])) * baseMVA
                newQD = -sum(imag(powerRefCaseL[areaIPQbuses])) * baseMVA
            end
            areaInfoI["shiftGen"] = 0
            areaInfoI["shiftLoad"] = -1
            existPV = 0
        else
            indE = bordInArea
            newTypes = []
            areaInfoI["shiftGen"] = 0
            areaInfoI["shiftLoad"] = 0
            existPV = 0
            existPQ = 0
        end
        # storing the essential border indices (new indices)
        areaInfoI["indEnew"] = indE

        # Non essential buses: all the others
        areaBuses = [1:totalSizeREI...]
        indNE = setdiff(areaBuses, indE)
        areaInfoI["indNEnew"] = indNE # storing non essential buses (new indices)

        # We now have to update lineInfos: getting the internal indices of the
        # border buses
        (_, indBordFr) = ismember(branch[bordLineFr, F_BUS], indAreaI[bordInArea])
        (_, indBordTo) = ismember(branch[bordLineTo, T_BUS], indAreaI[bordInArea])
        lineInfos[bordLineFr, 3] = indBordFr
        lineInfos[bordLineTo, 4] = indBordTo
        lineInfos[bordLineFr, 5] = areaInfoI["indE"][indBordFr]
        lineInfos[bordLineTo, 6] = areaInfoI["indE"][indBordTo]
        lineInfos[bordLineTo, 10] = findall(x -> x, bordLineTo)

        # Voltage at the REI PV bus (if it exists)
        Vtot_pv = []
        # Construct the zero-loss network for PV and PQ buses
        if (~isempty(areaIPVbuses))
            if selectPV
                # Stot
                Stot_pv = sum(powerRefCase[areaIPVbuses])
                # Y0i
                Y0i_pv = -conj(powerRefCase[areaIPVbuses]) ./ abs.(voltagesRefCase[areaIPVbuses]).^2
                # Ii
                ii_pv = conj(powerRefCase[areaIPVbuses] ./ voltagesRefCase[areaIPVbuses])
            else
                # Stot
                Stot_pv = sum(powerRefCaseG[areaIPVbuses])
                # Y0i
                Y0i_pv = -conj(powerRefCaseG[areaIPVbuses]) ./ abs.(voltagesRefCase[areaIPVbuses]).^2
                # Ii
                ii_pv = conj(powerRefCaseG[areaIPVbuses] ./ voltagesRefCase[areaIPVbuses])
            end
            # Itot
            itot_pv = sum(ii_pv)
            # Vtot
            Vtot_pv = Stot_pv / conj(itot_pv)
            # Ytot
            Ytot_pv = conj(Stot_pv) / abs(Vtot_pv)^2
        end

        if (~isempty(areaIPQbuses))
            if selectPV
                Stot_pq = sum(powerRefCase[areaIPQbuses])
                Y0i_pq = -conj(powerRefCase[areaIPQbuses]) ./ abs.(voltagesRefCase[areaIPQbuses]).^2
                ii_pq = conj(powerRefCase[areaIPQbuses] ./ voltagesRefCase[areaIPQbuses])
            else
                Stot_pq = sum(powerRefCaseL[areaIPQbuses])
                Y0i_pq = -conj(powerRefCaseL[areaIPQbuses]) ./ abs.(voltagesRefCase[areaIPQbuses]).^2
                ii_pq = conj(powerRefCaseL[areaIPQbuses] ./ voltagesRefCase[areaIPQbuses])
            end
            itot_pq = sum(ii_pq)
            Vtot_pq = Stot_pq / conj(itot_pq)
            Ytot_pq = conj(Stot_pq) / abs(Vtot_pq)^2
        end

        # Saving the values of the voltages of the PV buses
        areaInfoI["Vtot"] = [voltagesRefCase[sort([areaIBordPV; areaIBordRef])]; Vtot_pv]
        areaInfoI["Vtotpq"] = Vtot_pq

        # Filling in the info about the buses
        # column 1: bus type
        # column 2: PD
        # column 3: QD
        # column 4: cell array containing the lists of the generators connected
        # to the new REI bus.
        busTypes = [bus[indBusBord, BUS_TYPE]; newTypes]
        busPD = [bus[indBusBord, PD] .* baseMVA; newPD]
        busQD = [bus[indBusBord, QD] .* baseMVA; newQD]
        (busGen, _) = ismember(gen[:, GEN_BUS], areaIPVbuses)

        (bT_2, bT_3) = (busTypes .== 2, busTypes .== 3)
        nGenNew = nGenNew + length(busTypes[bT_2 .| bT_3])

        # Storing bus information
        areaInfoI["busTypes"] = busTypes
        areaInfoI["busPD"] = busPD
        areaInfoI["busQD"] = busQD
        areaInfoI["REIPG"] = newPG
        areaInfoI["REIQG"] = newQG
        areaInfoI["busGen"] = gen[busGen, GEN_BUS]

        # Building the new admittance matrix
        YadmAreaI = sparse(zeros(Complex, totalSizeREI, totalSizeREI))

        # Copying the original admittance matrix in the upper right-hand corner
        # of the matrix
        YadmAreaI[1:sizeAreaI, 1:sizeAreaI] = Yadm[indAreaI, indAreaI]

        # Finding the indices of the PV and PQ buses in the area
        if ((~isempty(areaIPVbuses)) && (~isempty(areaIPQbuses)))
            (_, pvInd) = ismember(areaIPVbuses, indAreaI)
            # Filling the admittance matrix for the 0 buses (symmetric matrix)
            YadmAreaI[pvInd, sizeAreaI + 1] = -Y0i_pv
            YadmAreaI[sizeAreaI + 1, pvInd] = -Y0i_pv
            # Filling for the total equivalent buses
            YadmAreaI[sizeAreaI + 1, sizeAreaI + 3] = -Ytot_pv
            YadmAreaI[sizeAreaI + 3, sizeAreaI + 1] = -Ytot_pv
            (_, pqInd) = ismember(areaIPQbuses, indAreaI)
            YadmAreaI[pqInd, sizeAreaI + 2] = -Y0i_pq
            YadmAreaI[sizeAreaI + 2, pqInd] = -Y0i_pq
            YadmAreaI[sizeAreaI + 2, sizeAreaI + 4] = -Ytot_pq
            YadmAreaI[sizeAreaI + 4, sizeAreaI + 2] = -Ytot_pq

        elseif (~isempty(areaIPQbuses))
            (_, pqInd) = ismember(areaIPQbuses, indAreaI)
            YadmAreaI[pqInd, sizeAreaI + 1] = -Y0i_pq
            YadmAreaI[sizeAreaI + 1, pqInd] = -Y0i_pq
            YadmAreaI[sizeAreaI + 1, sizeAreaI + 2] = -Ytot_pq
            YadmAreaI[sizeAreaI + 2, sizeAreaI + 1] = -Ytot_pq

        elseif (~isempty(areaIPVbuses))
            (_, pvInd) = ismember(areaIPVbuses, indAreaI)
            # Filling the admittance matrix for the 0 buses (symmetric matrix)
            YadmAreaI[pvInd, sizeAreaI + 1] = -Y0i_pv
            YadmAreaI[sizeAreaI + 1, pvInd] = -Y0i_pv
            # Filling for the total equivalent buses
            YadmAreaI[sizeAreaI + 1, sizeAreaI + 2] = -Ytot_pv
            YadmAreaI[sizeAreaI + 2, sizeAreaI + 1] = -Ytot_pv
        end
        # areaInfoI["busTypes"] = busTypes :TODO check if that is required

        # Modifying the diagonal, not forgetting the shunts and B from the
        # lines
        YdiagArea = [Ydiag[indAreaI]; zeros(totalSizeREI - sizeAreaI, 1)]
        YadmAreaI[(I(totalSizeREI) .== 1)] = YdiagArea

        for d in 1:totalSizeREI
            YadmAreaI[d, d] = YadmAreaI[d, d] - sum(YadmAreaI[d, 1:(d-1)]) -
            sum(YadmAreaI[d, (d+1):totalSizeREI])
        end

        # Gaussian reduction - we use pinv to take care of cases when the
        # admittance matrix is singular
        # real(sum(inv(YadmAreaI(indNE,indNE)))) :TODO check this comment
        YadmAreaRed = YadmAreaI[indE, indE] - YadmAreaI[indE, indNE] *
                      inv(Matrix(YadmAreaI[indNE, indNE])) * YadmAreaI[indNE, indE]

        # Here we store the reduced REI matrix, the full REI admittance matrix
        # and the original admittance matrix
        areaInfoI["YadmRed"] = YadmAreaRed
        areaInfoI["YadmREI"] = YadmAreaI
        areaInfoI["YadmOri"] = Yadm[indAreaI, indAreaI]

        # Updating the reduced system's size
        nNew += size(YadmAreaRed, 1)

        # storing the indices of the PV buses (for the new REI pv buses and the
        # non aggregated PV buses, i.e. those lying on the border).
        if existPV == 1
            indREIpv = size(YadmAreaRed, 1) - existPQ
        else
            indREIpv = []
        end

        areaInfoI["areaIPV"] = sort([areaIBordRef_int; areaIBordPV_int; indREIpv])

        # Storing the area data
        areaInfo.data[i] = areaInfoI
    end

    # keep data regarding the reduction in .data[0]
    areaInfo.data[0] = Dict("nNew" => nNew, "nGenNew" => nGenNew)
end

## We can now build the new admittance matrix and create the new test case
let
    println("Building admittance matrix...")
    # We need the base quantities
    Vbase = bus_data["1"]["base_kv"]

    # The size is
    nNew = areaInfo.data[0]["nNew"]
    YadmNew = sparse(zeros(Complex, nNew, nNew))

    # in indAreas, we will keep the first and last index of area i in the
    # admittance matrix YadmNew
    indAreas = zeros(2, areaInfo.no_of_areas)

    # Create the bus, gen and branch arrays
    busNew = zeros(nNew, BUS_ARRAY_SIZE)
    branchNew = zeros(Int(nNew * (nNew-1)/2), BR_ARRAY_SIZE)

    if genGroup
        genNew = zeros(nGenNew, GEN_ARRAY_SIZE)
    else
        genNew = gen
    end

    # gencostNew = gencost :TODO check on generators cost

    # initializing the bus numbers
    busNew[:, BUS_ID] = 1:nNew

    # Some common values for the buses
    busNew[:, VM] .= 1
    busNew[:, VA] .= 1
    busNew[:, VMAX] .= 1.05
    busNew[:, VMIN] .= 0.95
    busNew[:, BASE_KV] .= Vbase # :TODO check on that
    busNew[:, ZONE] .= 1

    indBegin = 1
    # indices to fill in the genNew array
    indGenBegin = 1
    indGenEnd = 0

    for a in 1:areaInfo.no_of_areas

        ai = areaInfo.data[a]

        sizeAi = size(ai["YadmRed"], 1)
        indEnd = indBegin + sizeAi - 1

        # We fill the diagonal blocks corresponding to the admittance matrices
        # of the area.
        YadmNew[indBegin:indEnd, indBegin:indEnd] = ai["YadmRed"]

        # Filling in the infos about the buses
        busNew[indBegin:indEnd, BUS_AREA] .= a
        busNew[indBegin:indEnd, BUS_TYPE] .= ai["busTypes"]
        busNew[indBegin:indEnd, PD] .= ai["busPD"]
        busNew[indBegin:indEnd, QD] = ai["busQD"]

        # indices of the PV bus to set voltages
        indPV_int = ai["areaIPV"][] + indBegin - 1
        busNew[indPV_int, VM] = abs(ai["Vtot"][])
        busNew[indPV_int, VA] = angle(ai["Vtot"][]) # :TODO check degs or rads -> 180/pi *

        # Copying the BS and GS values of the border buses from the original
        # buses
        indBorderEnd = indBegin + size(ai["indE"], 1) - 1
        indBorderOri = ai["indE"]
        indNonBorderOri = ai["indNE"]
        shunts = pm.data["shunt"]
        if ~isempty(shunts)
            # a temporary mapping of the shunts to buses
            bus_shunts = zeros(no_buses, 2)
            for sh in shunts
                id = parse(Int64, sh.second["shunt_bus"]) # :TODO check if is p.u. or not
                bus_shunts[id, :] = [sh.second["gs"], sh.second["bs"]]
            end
            for i in indBegin:indBorderEnd
                busNew[i, GS] = bus_shunts[indBorderOri[i], 1]
                busNew[i, BS] = bus_shunts[indBorderOri[i], 2]
            end
        end

        # The column loadMap stores the bus number to which each original load is
        # moved in the reduced system
        loadMap = zeros(no_buses, 1)

        # We can now fill in the column mapping the original buses to the new
        # ones (for the load):
        loadMap[indBorderOri] .= indBegin:indBorderEnd
        if (ai["shiftLoad"] == -1) # otherwise there is no REI load bus
            loadMap[indNonBorderOri] .= indEnd + ai["shiftLoad"] + 1
        end

        # Filling info about the generators lying on the border buses
        (cond1, cond2) = (ai["busTypes"][1:length(ai["indE"])] .== 2, ai["busTypes"][1:length(ai["indE"])] .== 3)
        condition = cond1 .| cond2
        indGenBord_log = findall(x -> x, condition)
        (indGenBord, _) = ismember(gen[:, GEN_BUS], ai["indE"][indGenBord_log])

        if genGroup
            if sum(indGenBord) > 0
                indGenEnd = indGenBegin + sum(indGenBord) - 1
                genNew[indGenBegin:indGenEnd, :] = gen[indGenBord, :]
                genNew[indGenBegin:indGenEnd, GEN_BUS] = indBegin .+ indGenBord_log .- 1 # :TODO this doesn't seem right
                # Update the indices if any change has been made
                indGenBegin = indGenEnd + 1
                indGenEnd = indGenEnd + 1
            end
        else
            # Connecting them to the right buses in the reduced model
            idxBusGen_areai = ai["indE"][indGenBord_log] # getting the PV buses
            for i in 1:length(idxBusGen_areai)
                # Getting all generators connected to bus idxBusGen_areai(i)
                (indGenBord_i, _) = ismember(gen[:, GEN_BUS], idxBusGen_areai[i])
                # Changing the bus number of all these generators
                genNew[indGenBord_i, GEN_BUS] = indBegin + indGenBord_log[i] - 1 # :TODO check this as well
            end
        end

        # Filling info about the generators aggregated to the new REI bus (if any)
        if ~isempty(ai["REIPG"])
            (indGen, _) = ismember(gen[:, GEN_BUS], ai["busGen"])
            if genGroup
                genNew[indGenBegin, GEN_BUS] = indEnd + ai["shiftGen"] + 1
                genNew[indGenBegin, PG] = ai["REIPG"]
                genNew[indGenBegin, QG] = ai["REIQG"]
                genNew[indGenBegin, QMAX] = sum(gen[indGen, QMAX]) # :TODO check p.u.!!
                genNew[indGenBegin, QMIN] = -sum(abs.(gen[indGen, QMIN]))
                genNew[indGenBegin, VG] = busNew[indEnd + ai["shiftGen"] + 1, VM]
                genNew[indGenBegin, MBASE] = baseMVA
                genNew[indGenBegin, GEN_STATUS] = 1
                genNew[indGenBegin, PMAX] = sum(gen[indGen, PMAX])
                genNew[indGenBegin, PMIN] = sum(gen[indGen, PMIN])
                indGenBegin = indGenBegin + 1
                indGenEnd = indGenEnd + 1
            else
                # Getting the original buses to which the generators are connected
                indGenTo = genNew[indGen, GEN_BUS]
                # Updating their production by taking away the local load
                if selectPV
                    genNew[indGen, PG] = real(powerRefCase[Int64.(indGenTo)]) * baseMVA
                    genNew[indGen, QG] = imag(powerRefCase[Int64.(indGenTo)]) * baseMVA
                end
                # Connecting them to the righ buses in the reduced model
                genNew[indGen, GEN_BUS] .= indEnd + ai["shiftGen"] + 1
                # Moving the generators connected to buses that have become
                # PQ buses (because the net real power injections were
                # negative, i.e. the generators generate less than the local
                # load) to the common PV REI bus BUT deactivating them
                (indGenPV2PQ, _) = ismember(gen[:, GEN_BUS], ai["PVtoPQbuses"])
                genNew[indGenPV2PQ, GEN_BUS] .= indEnd + ai["shiftGen"] + 1
                genNew[indGenPV2PQ, GEN_STATUS] .= -1
            end
        end

        # Saving the indices and updating indBegin
        indAreas[1, a] = indBegin
        indAreas[2, a] = indEnd
        indBegin = indBegin + sizeAi
    end

    # Change the initial voltages of the generators to be equal to those to
    # the bus to which they are connected
    genBus = Int64.(genNew[:, GEN_BUS])
    genNew[:, VG] = busNew[genBus, VM]

    # We now have to fill the admittances of the inter area lines.
    borderLines = findall(x -> x != 0, lineInfos[:, 1]) # indices of the inter area lines
    # Getting the info about the border lines and adding two columns to store
    # the indices of the end buses in the numbering of YadmNew and two columns
    # to store the number of the border lines among all lines
    borderLinesInfo = lineInfos[borderLines, :]
    # At the end we may need to swap some of the lines
    swapFrTo = ones(size(borderLinesInfo, 1), 1)

    # To keep track of parallel lines and avoiding updating the corresponding
    # diagonal element each time, I keep the already processed buses in one
    # array and see whether a line has already been processed between these two
    # buses.
    dealtWith = zeros(length(borderLines), 2)

    for l in 1:length(borderLines)
        # Getting the right indices
        aFr = Int64.(borderLinesInfo[l, 1]) # area from
        aTo = Int64.(borderLinesInfo[l, 2]) # area to
        bFr = Int64.(borderLinesInfo[l, 3]) # internal numbering of bus from
        bTo = Int64.(borderLinesInfo[l, 4]) # internal numbering of bus to
        indFr = Int64.(indAreas[1, aFr] + bFr - 1) # index bus from in YadmNew
        indTo = Int64.(indAreas[1, aTo] + bTo - 1) # index bus to in YadmNew
        indSysFr = Int64.(borderLinesInfo[l, 5]) # index bus from in the original system
        indSysTo = Int64.(borderLinesInfo[l, 6]) # index bus to in the original system

        # Checking whether another line between the two same buses has already
        # been processed in which case the equivalent admittance of all
        # parallel lines between these buses has already been considered.
        condition = ((dealtWith[:, 1] .== indFr) .& (dealtWith[:, 2] .== indTo)) .|
                    ((dealtWith[:, 2] .== indFr) .& (dealtWith[:, 1] .== indTo))
        critDealtWith = findall(x -> x > 0, condition)
        if isempty(critDealtWith) # :TODO check if nonzeros() is required
            # Updating the new admittance matrix
            YadmNew[indFr, indTo] = Yadm[indSysFr, indSysTo]
            YadmNew[indTo, indFr] = Yadm[indSysFr, indSysTo]
            # Updating the diagonal
            YadmNew[indFr, indFr] = YadmNew[indFr, indFr] - Yadm[indSysFr, indSysTo] # :TODO check those values
            YadmNew[indTo, indTo] = YadmNew[indTo, indTo] - Yadm[indSysFr, indSysTo]
            # Keeping the processed buses
            dealtWith[l, 1] = indFr
            dealtWith[l, 2] = indTo
        else
            print("The buses $(dealtWith[critDealtWith, 1]) and ")
            print("$(dealtWith[critDealtWith, 2]) were already encountered ")
            print("for line(s) $(critDealtWith)\n")
        end
        # Storing the indices of the buses from and to according to the
        # numbering in YadmNew
        borderLinesInfo[l, 7] = indFr
        borderLinesInfo[l, 8] = indTo
    end

    # Filling in the shunt values
    busNew[:, GS] = real(sum(YadmNew, dims=1)) * baseMVA #Shunt susceptance
    busNew[:, BS] = imag(sum(YadmNew, dims=1)) * baseMVA # shunt values must NOT be in pu

    # Taking care of the branches: for each bus i, we have n-i possible
    # branches. Some of them will be offline, status=0, meaning that there is
    # no such branch in reality.
    indBrBegin = 1
    for i in 1:nNew
        indBrEnd = indBrBegin + (nNew-i) - 1
        branchNew[indBrBegin:indBrEnd, F_BUS] .= i # all branches that come from bus i
        branchNew[indBrBegin:indBrEnd, T_BUS] = (i+1):nNew
        for j in (i+1):nNew
            if YadmNew[i, j] != 0
                # First I check whether this is an inter-area branch and, if it
                # is the case, I check whether the two buses are connected by
                # parallel branches
                cond = (((borderLinesInfo[:, 7] .== i) .& (borderLinesInfo[:, 8] .== j)) .|
                       ((borderLinesInfo[:, 7] .== j) .& (borderLinesInfo[:, 8] .== i)))
                indBorderLine2 = findall(x -> x > 0, cond)

                # If one of the border line goes from bus i to bus j and from
                # bus j to bus i
                if ~isempty(indBorderLine2)
                    println("Line non considere parce que elle est entre deux zones:($(i),$(j))")
                else
                    branchNew[indBrBegin + j - (i+1), BR_R] = -real(1/YadmNew[i, j])
                    branchNew[indBrBegin + j - (i+1), BR_X] = -imag(1/YadmNew[i, j])
                    branchNew[indBrBegin + j - (i+1), BR_STATUS] = 1
                end
            end
        end
        #indNonZero = ~ismember(YadmNew(i,(i+1):nNew),0);
        indBrBegin = indBrBegin + (nNew-i)
    end
    branchNew[:, ANGMIN] .= -60 / 180 * π # PowerModels accepts [-60, 60] degs
    branchNew[:, ANGMAX] .= 60 / 180 * π
    ln_non_def = findall(x -> x > 0, branchNew[:, BR_STATUS])
    branchNew = branchNew[ln_non_def, :] # removing the lines non defined

    # Finding the indices of the border branches
    branchNewBorder = zeros(size(borderLinesInfo, 1), BR_ARRAY_SIZE)

    for i in 1:size(borderLinesInfo, 1)
        # Copying the original inter-area branches
        br_i = branch_data["$(Int(borderLinesInfo[i,10]))"]
        branchNewBorder[i, :] = readBrDict(br_i)
        branchNewBorder[i, BR_B] = 0
        branchNewBorder[i, F_BUS] = borderLinesInfo[i, 7]
        branchNewBorder[i, T_BUS] = borderLinesInfo[i, 8]

        # Saving the indices of the reduced border lines in the branch array
        borderLinesInfo[i, 9] = size(branchNew, 1) + i # :TODO why?
    end

    branchNew = [branchNew; branchNewBorder]

    if (pfMethod == 2)
        busNew[:, VMAX] .= 1.06
        busNew[:, VMIN] .= 0.95
    end

    # Some of the lines must be swapped to be consistent when it comes to the
    # from and to buses
    indSwap = findall(x -> x == -1, swapFrTo)
    if ~isempty(indSwap)
        borderLinesInfo[indSwap, [1, 2, 3, 4, 5, 6]] = borderLinesInfo[indSwap, [2, 1, 4, 3, 6, 5]]
    end

    ## Save the reduced case :TODO
    mpcNew = Dict()
    mpcNew["baseMVA"] = baseMVA
    mpcNew["bus"] = busNew
    mpcNew["branch"] = branchNew
    mpcNew["gen"] = genNew
    mpcNew["gencost"] = [] # :TODO must implement that
    loadmap_name = string("$(caseName)_reduced.m")
    areaInfo.data[0]["pm_reduced"] = mpcNew
end

save = false

if save
    pmodel = export_matpower(mpcNew)
    # save model in file
    s = sprint(print, pmodel)
    write(loadmap_name, s)
end
