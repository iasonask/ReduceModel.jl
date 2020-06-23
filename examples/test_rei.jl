## Test the REI functionality
using ReduceModel
using Ipopt

# load network_data
file = joinpath(dirname(@__FILE__), "cases/case118.m")

# calculate rei and return a PowerModel dict,
# choosing number of areas and default values
case = call_rei(file, 4; optimizer=Ipopt.Optimizer)

## use custom preferences
using ReduceModel
using PowerModels
using Ipopt

rei_opt = REIOptions(ACPPowerModel, build_opf, true, true)

file = joinpath(dirname(@__FILE__), "cases/case300.m")

original_net = parse_file(file)

reduced_net = call_rei(
    file,
    2;
    options = rei_opt,
    optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0),
    export_file = false,
    path = "./examples",
    no_tries = ReduceModel.NO_TRIES,
)

pl = makePlots(original_net, reduced_net)

## Plots
using ReduceModel
using PowerModels
using Ipopt
using Plots
# using Pkg; Pkg.build("GR")
plotlyjs()

rei_opt = REIOptions(ACPPowerModel, build_pf, false, false)

file = joinpath(dirname(@__FILE__), "cases/case_ACTIVSg2000.m")

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
savefig(pl, joinpath(dirname(@__FILE__), "../results/reduced_plot_118.pdf"))


## check active power losses error
sol_red = run_ac_opf(reduced_net, Ipopt.Optimizer)
sol_ori = run_ac_opf(original_net, Ipopt.Optimizer)

red_loss = [loss[2]["pt"]+loss[2]["pf"] for loss in sol_red["solution"]["branch"]]
total_red_loss = sum(red_loss)

ori_loss = [loss[2]["pt"]+loss[2]["pf"] for loss in sol_ori["solution"]["branch"]]
total_ori_loss = sum(ori_loss)

loss_err = ((total_ori_loss - total_red_loss) / total_ori_loss) * 100
