"""
     reduce_network(areaInfo::PMAreas, options::REIOptions) -> case::Dict{String, Any}

Create a Reduced Equivalent Independent power model for each area


    # Build the new admittance matrix and create the new test case


"""

function reduce_network(areaInfo::PMAreas, options::REIOptions)

    pm = areaInfo.data[0]["pm"]
    no_buses = areaInfo.data[0]["no_buses"]
    gen = areaInfo.data[0]["gen"]
    genint2ext = areaInfo.data[0]["genint2ext"]
    lineInfos = areaInfo.data[0]["lineInfos"]
    Yadm = areaInfo.data[0]["Yadm"]
    selectPV = options.selectPV
    genGroup = options.genGroup
    network_data = areaInfo.data[0]["original"]
    powerRefCase = areaInfo.data[0]["powerRefCase"]

    ext2int = Dict(bus["index"] => i for (i, (k, bus)) in enumerate(sort(network_data["bus"], by=x->parse(Int, x))))
    println("Building admittance matrix...")
    # We need the base quantities
    Vbase = areaInfo.data[0]["base_kv"]

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
        genNew = zeros(areaInfo.data[0]["nGenNew"], GEN_ARRAY_SIZE)
        # fill with ones the GEN_BUS column for avoiding errors when the
        # clustering solution is problematic
        genNew[:, GEN_BUS] .= 1
    else
        genNew = areaInfo.data[0]["gen"]
    end

    # gencostNew = gencost :TODO check on generators cost

    # initializing the bus numbers
    busNew[:, BUS_ID] = 1:nNew

    # Some common values for the buses
    busNew[:, VM] .= 1
    busNew[:, VA] .= 0
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
        try
            busNew[indBegin:indEnd, PD] = ai["busPD"]
            busNew[indBegin:indEnd, QD] = ai["busQD"]
        catch e
           println("Problematic clustering solution, reducing the number of clusters could help.")
           busNew[indBegin:indEnd, PD] .= 0
           busNew[indBegin:indEnd, QD] .= 0
        end

        # indices of the PV bus to set voltages
        # display(ai["areaIPV"])
        indPV_int = ai["areaIPV"] .+ indBegin .- 1
        busNew[indPV_int, VM] = abs.(ai["Vtot"])
        busNew[indPV_int, VA] = angle.(ai["Vtot"]) # in rad -> / π * 180 for degs

        # Copying the BS and GS values of the border buses from the original
        # buses
        indBorderEnd = indBegin + size(ai["indE"], 1) - 1
        indBorderOri = ai["indE"]
        indNonBorderOri = ai["indNE"]
        shunts = pm.data["shunt"]
        if ~isempty(shunts)
            # a temporary mapping of the shunts to buses
            # should be implement more efficient
            bus_shunts = zeros(no_buses, 2)
            for sh in shunts
                id = ext2int[Int64(sh.second["shunt_bus"])] # :TODO check if is p.u. or not
                bus_shunts[id, :] = [sh.second["gs"], sh.second["bs"]]
            end
            for (i, _) in enumerate(indBegin:indBorderEnd)
                busNew[i, GS] = bus_shunts[indBorderOri[i], 1]
                busNew[i, BS] = bus_shunts[indBorderOri[i], 2]
            end
        end

        # The column loadMap stores the bus number to which each original load is
        # moved in the reduced system
        loadMap = zeros(Int64, no_buses)

        # We can now fill in the column mapping the original buses to the new
        # ones (for the load):

        loadMap[indBorderOri] = collect(indBegin:indBorderEnd)
        if (ai["shiftLoad"] == -1) # otherwise there is no REI load bus
            loadMap[indNonBorderOri] .= indEnd + ai["shiftLoad"] + 1
        end

        # keep buses Mapping in areaInfo as well
        if ~haskey(areaInfo.data[0], "mapping")
            areaInfo.data[0]["mapping"] = loadMap
        else
            I = findall(x -> x > 0, loadMap)
            areaInfo.data[0]["mapping"][I] = loadMap[I]
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
                (indGenBord_i, _) = ismember(gen[:, GEN_BUS], [idxBusGen_areai[i]])
                # Changing the bus number of all these generators
                genNew[indGenBord_i, GEN_BUS] .= indBegin + indGenBord_log[i] - 1 # :TODO check this as well
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
                genNew[indGenBegin, MBASE] = network_data["baseMVA"]
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
                    genNew[indGen, PG] = real(powerRefCase[Int64.(indGenTo)])
                    genNew[indGen, QG] = imag(powerRefCase[Int64.(indGenTo)])
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
        indBegin += sizeAi
    end

    # Change the initial voltages of the generators to be equal to those to
    # the bus to which they are connected
    genBus = Int64.(genNew[:, GEN_BUS])
    # genBus = [gen for gen in genBus]
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
    busNew[:, GS] = real(sum(YadmNew, dims=1)) # Shunt susceptance -> * baseMVA for matpower
    busNew[:, BS] = imag(sum(YadmNew, dims=1)) # :OLD "shunt values must NOT be in pu"

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
                    #println("Line not considered:($(i),$(j))")
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
        br_i = pm.data["branch"]["$(Int(borderLinesInfo[i,10]))"]
        branchNewBorder[i, :] = readBrDict(br_i)
        branchNewBorder[i, BR_B] = 0
        branchNewBorder[i, F_BUS] = borderLinesInfo[i, 7]
        branchNewBorder[i, T_BUS] = borderLinesInfo[i, 8]

        # Saving the indices of the reduced border lines in the branch array
        borderLinesInfo[i, 9] = size(branchNew, 1) + i # :TODO why?
    end

    branchNew = [branchNew; branchNewBorder]
    busNew[:, VMAX] .= 1.06
    busNew[:, VMIN] .= 0.95

    # Some of the lines must be swapped to be consistent when it comes to the
    # from and to buses
    indSwap = findall(x -> x == -1, swapFrTo)
    if ~isempty(indSwap)
        borderLinesInfo[indSwap, [1, 2, 3, 4, 5, 6]] = borderLinesInfo[indSwap, [2, 1, 4, 3, 6, 5]]
    end

    ## create a PowerModels dict to store the reduced network
    caseName = string("$(pm.data["name"])_reduced")
    # prepare a powermodels network dict and populate it with the calculated values
    case = Dict{String, Any}()
    case["bus"] = Dict{String,Any}()
    for id in 1:size(busNew, 1)
        case["bus"]["$(id)"] = busArrayToDict(busNew[id, :])
    end

    case["source_type"] = "matpower"
    case["name"] = caseName
    case["dcline"] = Dict{String, Any}() # will be left empty
    case["source_version"] = "2"

    case["gen"] = Dict{String, Any}()
    if genGroup
        @warn "Aggregating the generators, no cost parameters are introduced!"
        for id in 1:size(genNew, 1)
            case["gen"]["$(id)"] = genArrayToDict(id, genNew[id, :])
        end
    else
        for id in 1:size(genNew, 1)
            case["gen"]["$(id)"] = genArrayToDict(id, genNew[id, :], pm.data["gen"]["$(genint2ext[id])"])
        end
    end

    case["branch"] = Dict{String, Any}()
    for id in 1:size(branchNew, 1)
        case["branch"]["$(id)"] = branchArrayToDict(id, branchNew[id, :])
    end

    if ~isempty(pm.data["storage"])
        @warn "\"storage\" dict has storage units which are not considered here."
    end
    case["storage"] = Dict{String, Any}() # will be left empty

    if ~isempty(pm.data["switch"])
        @warn "\"switch\" dict has storage units which are not considered here."
    end
    case["switch"] = Dict{String, Any}() # will be left empty

    case["baseMVA"] = pm.data["baseMVA"]
    case["per_unit"] = true

    case["shunt"] = Dict{String, Any}()
    sh_id = 1
    for (bus_id, (gs, bs)) in enumerate(zip(busNew[:, GS], busNew[:, BS]))
        if (abs(gs) + abs(bs) > 0)
            case["shunt"]["$(sh_id)"] = shuntsToDict(sh_id, bus_id, busNew[bus_id, :])
            sh_id += 1
        end
    end

    case["load"] = Dict{String, Any}()
    ld_id = 1
    for (bus_id, (pd, qd)) in enumerate(zip(busNew[:, PD], busNew[:, QD]))
        if (abs(pd) + abs(qd) > 0)
            case["load"]["$(ld_id)"] = loadToDict(ld_id, bus_id, busNew[bus_id, :])
            ld_id += 1
        end
    end

    PowerModels.correct_network_data!(case)
    # append also area information to case dict
    case["areas"] = areaInfo
    case
end
