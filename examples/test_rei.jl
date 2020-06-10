## Test the REI functionality
using ReduceModel

# load network_data
# file = joinpath(dirname(@__FILE__), "cases/case_ACTIVSg2000.m")
# file = joinpath(dirname(@__FILE__), "cases/case_ACTIVSg500.m")
# file = joinpath(dirname(@__FILE__), "cases/case39.m")
# file = joinpath(dirname(@__FILE__), "cases/case300.m")

# calculate rei and return a PowerModel dict,
# choosing number of areas and default values
# case = call_rei(file, 4)

## use custom preferences

using PowerModels
using Ipopt
rei_opt = REIOptions(ACPPowerModel, build_opf, true, false)


file = joinpath(dirname(@__FILE__), "cases/case118.m")
reduced_case = call_rei(
    file,
    3;
    options = rei_opt,
    optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0),
    export_file = false,
    path = "./examples",
    no_tries = ReduceModel.NO_TRIES,
)


## Plots

using ReduceModel
using PowerModels, Plots
using Ipopt
using Pkg; Pkg.build("GR")

rei_opt = REIOptions(ACPPowerModel, build_pf, false, false)


file = joinpath(dirname(@__FILE__), "cases/case118.m")
# file = joinpath(dirname(@__FILE__), "cases/case_ACTIVSg2000.m")
# file = joinpath(dirname(@__FILE__), "cases/case300.m")
# file = joinpath(dirname(@__FILE__), "cases/case39.m")


original_net = parse_file(file)

reduced_net = call_rei(
    file,
    5;
    options = rei_opt,
    optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0),
    export_file = false,
    path = "./examples",
)

pl = makePlots(original_net, reduced_net)
# display(pl)
