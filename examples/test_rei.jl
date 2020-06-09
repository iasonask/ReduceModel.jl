## Test the REI functionality
using SpineReduce

# load network_data
# file = "examples/case118.m"
# file = "examples/case_ACTIVSg500.m"
# file = "examples/case_ACTIVSg2000.m"
# file = "examples/case39.m"
file = "examples/case300.m"

# calculate rei and return the result in PowerModel dict, choosing number
# of areas and all other default values
case = call_rei(file, 4)
