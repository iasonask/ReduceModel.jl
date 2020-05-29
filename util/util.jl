using SparseArrays


## admittance matrix calculation
# from: https://kersulis.github.io/2015/11/05/admittance-matrix/

"""
    createYbus(f,t,x [,r,b]) -> Y
Create an admittance matrix for AC power flow.
All inputs are real. The output matrix is real if no line
resistances are provided (DC case), and complex otherwise.
* `f`,`t`: vectors encoding all lines (fi,ti)
* `x`: per-unit reactance xi for all lines
* `r`: per-unit resistance ri for all lines
* `b`: per-unit susceptance bi for all lines
"""
function createYbus(
    f::Vector{Int64},
    t::Vector{Int64},
    x::Vector{Float64},
    r=0.0::Union{Vector{Float64},Float64},
    b=0.0::Union{Vector{Float64},Float64}
    )
    z = r + x * 1im
    y = 1 ./ z
    b = b * 1im
    Y = sparse([f; t; t; f],[t; f; t; f],[-y; -y; y + b./2; y + b./2])

    # for DC power flow, we typically want a matrix with real entries:
    if r == 0
        return imag(Y)
    else
        return Y
    end
end
