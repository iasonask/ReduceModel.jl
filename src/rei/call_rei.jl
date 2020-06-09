"""
    call_rei(pm::PowerModel, areas::PMAreas, options::REIOptions ) -> case::PowerModel

Create a Reduced Equivalent Independent power model


    # Arguments


"""

function call_rei(
    file::String,
    no_areas::Int64;
    options::REIOptions=REIOptions(),
    optimizer=Ipopt.Optimizer,
    export_file=true,
    path=@__DIR__,
    no_tries=NO_TRIES,
    )

    println("Starting REI calculation...")
    network_data = PowerModels.parse_file(file)
    caseName = split(file, "/")[end]
    # Set power flow model
    PFModel = options.pf_model
    _pf = options.pf_method

    # reduced model
    case = Dict{String, any}

    # areas
    ext2int = Dict(bus["index"] => i for (i, (k, bus)) in enumerate(sort(network_data["bus"], by=x->parse(Int, x))))

    # the non-deterministic clustering solution might cause admittance singularity
    # issues, the procedure is repeated no_tries times for improving the chances
    # of calculating a convergent power flow model
    for tr in 1:no_tries
        println("Trying to calculate REI, iteration: $(tr).")
        areas = partition(no_areas, ext2int, network_data["bus"], network_data["branch"])
        # a PMAreas data structure holds all relevant parameters and values for
        # the reduction procedure
        areaInfo = PMAreas(:cluster, no_areas, areas)

        # calculate the reduced network of each area
        aggregateAreas!(areaInfo, network_data, options, optimizer)

        # combine reduced areas back together to a single power flow model
        case = reduce_network(areaInfo, options)

        # test that the calculated model converges
        pm_red = instantiate_model(case, PFModel, _pf)
        results_red = optimize_model!(pm_red, optimizer=optimizer)

        if termination_status(pm_red.model) == MOI.LOCALLY_SOLVED ||
            termination_status(pm_red.model) == MOI.OPTIMAL
            println("Calulation of REI completed successfully after $(tr) try(ies)!")
            break
        else
            println("Calulation of REI failed: $(termination_status(pm_red.model))," * (tr < NO_TRIES ?  " trying again." : " exiting."))
        end
    end

    # export matpower file
    if export_file
        pmodel_string = export_matpower(case)
        # save model in file
        s = sprint(print, pmodel_string)
        println("Writing results case...")
        write(joinpath(path, case["name"]*".m"), s)
    end
    case
end
