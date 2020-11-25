"""
     aggregateAreas!(areaInfo::PMAreas, network_data::Dict{String,Any}, options::REIOptions) -> areaInfo::PMAreas

Create a Reduced Equivalent Independent power model for each area


    # Arguments


"""

function aggregateAreas!(areaInfo::PMAreas, pm::ACPPowerModel, options::REIOptions)

    network_data = pm.data
    selectPV = options.selectPV
    genGroup = options.genGroup

    # Similar to Matpower's functions
    ext2int = Dict(bus["index"] => i for (i, (k, bus)) in enumerate(sort(network_data["bus"], by=x->parse(Int, x))))
    int2ext = Dict(i => bus["index"] for (i, (k, bus)) in enumerate(sort(network_data["bus"], by=x->parse(Int, x))))
    # Initializations
    # number of buses
    no_buses = length(network_data["bus"])
    no_branches = length(network_data["branch"])
    # remove decommitioned generators? :TODO check also unused branches
    for gn in network_data["gen"]
        if gn.second["gen_status"] == 0
            delete!(network_data["gen"], gn.first)
        end
    end
    genext2int = Dict(gen["index"] => i for (i, (k, gen)) in enumerate(sort(network_data["gen"], by=x->parse(Int, x))))
    genint2ext = Dict(i => gen["index"] for (i, (k, gen)) in enumerate(sort(network_data["gen"], by=x->parse(Int, x))))
    no_gen = length(network_data["gen"])

    baseMVA = network_data["baseMVA"]
    bus_full = network_data["bus"]
    gen_full = network_data["gen"]
    branch_full = network_data["branch"]

    # Get the indices of the Ref, PV and PQ buses
    indRef = [ext2int[bus.second["bus_i"]] for bus in bus_full if bus.second["bus_type"] == 3]
    indPV = sort([ext2int[bus.second["bus_i"]] for bus in bus_full if bus.second["bus_type"] == 2])
    indPQ = sort([ext2int[bus.second["bus_i"]] for bus in bus_full if bus.second["bus_type"] == 1])

    # calculate admittance matrix of full network
    Ybus = calc_admittance_matrix(network_data)
    Yadm = Ybus.matrix
    conj!(Yadm) # PowerModels returns the conjugate for some reason
    # sum each column to get the diagonal elements
    Ydiag = sum(Yadm, dims=1)

    # Save area information in PowerModel object
    for area in areaInfo.clusters
        for bus in area.second
            bus_full["$(int2ext[bus])"]["area"] = area.first
        end
    end

    # Get voltages and apparent powers from the load flow case
    solution = pm.solution
    (baseMVA, bus_sol, gen_sol) = (solution["baseMVA"], solution["bus"], solution["gen"])
    (bus_data, gen_data, branch_data, load_data) = (pm.data["bus"], pm.data["gen"], pm.data["branch"], pm.data["load"])

    # pass the results from the dictionary to arrays for easier calculations
    bus = zeros(Float64, (no_buses, 9))
    # sorted buses in 1:no_buses
    for bs in bus_sol
        id = ext2int[parse(Int64, bs.first)]
        bus[id, BUS_ID:BUS_TYPE] = [id, bus_data[bs.first]["bus_type"]]
        bus[id, VM:VA] = [bs.second["vm"], bs.second["va"]]
    end

    # create branch table for simplifying operations
    branch = zeros(Int64, (no_branches, 2))
    for ln in branch_data
        id = parse(Int64, ln.first)
        branch[id, F_BUS:T_BUS] = [ext2int[ln.second["f_bus"]], ext2int[ln.second["t_bus"]]]
    end

    # generators
    gen = zeros(Float64, (no_gen, GEN_ARRAY_SIZE))
    # sorted buses in 1:no_buses
    for gn in gen_data
        id = genext2int[parse(Int64, gn.first)]
        gen[id, GEN_BUS:QG] = [ext2int[gn.second["gen_bus"]], gn.second["pg"], gn.second["qg"]]
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
            id = ext2int[gn.second["gen_bus"]]
            powerRefCaseG[id] += gen_sol[gn.first]["pg"] + 1im * gen_sol[gn.first]["qg"]
        end
    end

    # Fill in load
    for ld in load_data
        id = ext2int[ld.second["load_bus"]]
        powerRefCaseL[id] = -(ld.second["pd"] + 1im * ld.second["qd"])  # :TODO should it -= here instead?
        bus[id, PD:QD] = [ld.second["pd"], ld.second["qd"]]
    end

    powerRefCase = powerRefCaseG .+ powerRefCaseL

    # keep the indices for buses with non-zero injections
    # Take only x coordinates of CartesianIndex
    indNZ = findall(x -> x > 0, real.(powerRefCase))
    indNZ = [i[1] for i in indNZ]
    indPQNZ = findall(x -> x < 0, real.(powerRefCase))
    indPQNZ = [i[1] for i in indPQNZ]
    indRefNZ = intersect(indNZ, indRef)
    indPVNZ = setdiff(indNZ, indRefNZ)
    indPQori = findall(x -> x != 0, bus[:,PD])

    ## Construct the new area matrices
    # We then use the information stored in areaNames, to identify the
    # generator and load buses in each area.

    # We will also need to keep information about interarea lines.
    nl = no_branches
    lineInfos = zeros(nl,10)

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
            @warn "One of the aggregated generators was the slack bus!"
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
                newPD = [0; -sum(real(powerRefCase[areaIPQbuses]))] # for matpower -> * baseMVA
                # the minus sign is because we want the value of the consumption
                newQD = [0; -sum(imag(powerRefCase[areaIPQbuses]))] # for matpower -> * baseMVA
                newPG = sum(real(powerRefCase[areaIPVbuses])) # for matpower -> * baseMVA
                newQG = sum(imag(powerRefCase[areaIPVbuses])) # for matpower -> * baseMVA
            else
                newPD = [0; -sum(real(powerRefCaseL[areaIPQbuses]))] # for matpower -> * baseMVA
                # the minus sign is because we want the value of the consumption
                newQD = [0; -sum(imag(powerRefCaseL[areaIPQbuses]))] # for matpower -> * baseMVA
                newPG = sum(real(powerRefCaseG[areaIPVbuses])) # for matpower -> * baseMVA
                newQG = sum(imag(powerRefCaseG[areaIPVbuses])) # for matpower -> * baseMVA
            end

            areaInfoI["shiftGen"] = -2 # indicating that the new REI gen bus is the one before the last
            areaInfoI["shiftLoad"] = -1

        elseif (~isempty(areaIPVbuses))
            indE = [bordInArea; sizeAreaI+2]
            newTypes = maximum(bus[areaIPVbuses, BUS_TYPE])
            if selectPV
                newPG = sum(real(powerRefCase[areaIPVbuses])) # for matpower -> * baseMVA
                newQG = sum(imag(powerRefCase[areaIPVbuses])) # for matpower -> * baseMVA

                newPG = sum(real(powerRefCaseG[areaIPVbuses])) # for matpower -> * baseMVA
                newQG = sum(imag(powerRefCaseG[areaIPVbuses])) # for matpower -> * baseMVA
            end
            areaInfoI["shiftGen"] = -1 # indicating that the new REI gen bus is the one before the last
            areaInfoI["shiftLoad"] = 0 # indicating that the new REI load bus does not exist
            existPQ = 0

        elseif (~isempty(areaIPQbuses))
            indE = [bordInArea; sizeAreaI+2]
            newTypes = 1
            if selectPV
                newPD = -sum(real(powerRefCase[areaIPQbuses])) # for matpower -> * baseMVA
                newQD = -sum(imag(powerRefCase[areaIPQbuses])) # for matpower -> * baseMVA
            else
                newPD = -sum(real(powerRefCaseL[areaIPQbuses])) # for matpower -> * baseMVA
                newQD = -sum(imag(powerRefCaseL[areaIPQbuses])) # for matpower -> * baseMVA
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
        areaInfoI["Vtotpq"] = ~(@isdefined Vtot_pq) ? Vtot_pq = 1 + 0im : Vtot_pq

        # Filling in the info about the buses
        # column 1: bus type
        # column 2: PD
        # column 3: QD
        # column 4: cell array containing the lists of the generators connected
        # to the new REI bus.
        busTypes = [bus[indBusBord, BUS_TYPE]; newTypes]
        busPD = [bus[indBusBord, PD]; newPD]
        busQD = [bus[indBusBord, QD]; newQD]
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
        rtol = sqrt(eps(real(float(one(eltype(Matrix(YadmAreaI[indNE, indNE])))))))
        YadmAreaRed = YadmAreaI[indE, indE] - YadmAreaI[indE, indNE] *
                      pinv(Matrix(YadmAreaI[indNE, indNE]); rtol=rtol) * YadmAreaI[indNE, indE]

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
    areaInfo.data[0]["base_kv"] = bus_data["$(int2ext[indRef[]])"]["base_kv"]
    areaInfo.data[0]["gen"] = gen
    areaInfo.data[0]["genint2ext"] = genint2ext
    areaInfo.data[0]["pm"] = pm
    areaInfo.data[0]["no_buses"] = no_buses
    areaInfo.data[0]["lineInfos"] = lineInfos
    areaInfo.data[0]["Yadm"] = Yadm
    areaInfo.data[0]["original"] = network_data
    areaInfo.data[0]["powerRefCase"] = powerRefCase

    areaInfo
end
