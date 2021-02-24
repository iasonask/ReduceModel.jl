"""
    call_rei(pm::PowerModel, areas::PMAreas, options::REIOptions ) -> case::PowerModel

Create a Reduced Equivalent Independent power model


    # Arguments


"""

function call_rei(
    file::String,
    no_areas::Int64;
    optimizer,
    options::REIOptions = REIOptions(),
    export_file = true,
    path = @__DIR__,
    no_tries = NO_TRIES,
)

    println("Starting REI calculation...")
    network_data = PowerModels.parse_file(file)
    caseName = split(file, "/")[end]
    # Set power flow model
    PFModel = options.pf_model
    _pf = options.pf_method

    # reduced model
    case = Dict{String,any}

    # areas
    ext2int = Dict(
        bus["index"] => i
        for
        (i, (k, bus)) in
        enumerate(sort(network_data["bus"], by = x -> parse(Int, x)))
    )
    # create PF model
    pm = instantiate_model(network_data, PFModel, _pf)
    # run the power flow
    optimize_model!(pm, optimizer = optimizer)

    # the non-deterministic clustering solution might cause admittance singularity
    # issues, the procedure is repeated no_tries times for improving the chances
    # of calculating a convergent power flow model
    for tr = 1:no_tries
        println("Trying to calculate REI, iteration: $(tr).")
        areas = partition(
            no_areas,
            ext2int,
            network_data["bus"],
            network_data["branch"],
        )
        # a PMAreas data structure holds all relevant parameters and values for
        # the reduction procedure
        areaInfo = PMAreas(:cluster, no_areas, areas)

        try
            # calculate the reduced network of each area
            aggregateAreas!(areaInfo, pm, options)

            # combine reduced areas back together to a single power flow model
            case = reduce_network(areaInfo, options)

            # test that the calculated model converges
            pm_red = instantiate_model(case, PFModel, _pf)
            results_red = optimize_model!(pm_red, optimizer = optimizer)

            if results_red["termination_status"] == LOCALLY_SOLVED ||
               results_red["termination_status"] == OPTIMAL
                println(
                    "Calulation of REI completed successfully after $(tr) ",
                    tr > 1 ? "tries!" : "try!",
                )
                break
            else
                println(
                    "Calulation of REI PowerFlow failed: $(results_red["termination_status"])," *
                    (
                        tr < NO_TRIES ? " trying again." :
                        " exiting, consider choosing different options for the REI calculation."
                    ),
                )
            end

        catch e
            println(
                "Calulation of REI failed: possibly problematic REI configuration," *
                (
                    tr < NO_TRIES ? " trying again." :
                    " exiting, consider choosing different options for the REI calculation."
                ),
            )
            sprint(showerror, e)
        end
    end

    # export matpower file
    if export_file
        pmodel_string = export_matpower(case)
        # save model in file
        s = sprint(print, pmodel_string)
        println("Writing results case...")
        write(joinpath(path, case["name"] * ".m"), s)
    end
    case
end
