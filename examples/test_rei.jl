## Test the REI functionality
using ReduceModel

# load network_data
# file = "examples/case_ACTIVSg500.m"
# file = "examples/case_ACTIVSg2000.m"
# file = "examples/case39.m"
file = joinpath(dirname(@__FILE__), "cases/case300.m")

# calculate rei and return a PowerModel dict,
# choosing number of areas and default values
case = call_rei(file, 4)

## use custom preferences

using PowerModels
using Ipopt
rei_opt = REIOptions(ACPPowerModel, build_pf, true, false)


file = joinpath(dirname(@__FILE__), "cases/case118.m")
case2 = call_rei(file, 3;
                 options=rei_opt,
                 optimizer=optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0),
                 export_file=true,
                 path="./examples",
                 no_tries=ReduceModel.NO_TRIES,
                 )
