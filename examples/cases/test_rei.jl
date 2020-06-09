## Initializations
using ReduceModel

# load network_data
# network_data = PowerModels.parse_file("/Users/iasonas/Documents/MATLAB/matpower7.0/data/case118_mod.m")
# file = "/Users/iasonas/Documents/MATLAB/matpower7.0/data/case39.m"
file = "data/Matpower/case118.m"
# file = "data/Matpower/case_ACTIVSg500.m"
# file = "data/Matpower/case39.m"
# file = "data/Matpower/case2848rte.m"

case = call_rei(file, 4)
