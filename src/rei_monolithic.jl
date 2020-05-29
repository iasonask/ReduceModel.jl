using Revise
using PowerModels
using Ipopt

# REI monolithic script
println("Calculating REI...")

# load network_data
network_data = PowerModels.parse_file("../examples/case5.m")

# create PF model
pm = instantiate_model(network_data, ACPPowerModel, PowerModels.build_pf)


# calculate the admittance matric using PowerModels function !! it returns the conjugate !!
Ybus = calc_admittance_matrix(network_data)

# run OPF
result = optimize_model!(pm, optimizer=Ipopt.Optimizer)


print(pm.model)
print(result["solution"])
