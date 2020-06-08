using PyCall
using Clustering
using Plots

# required python packages
# numpy
# matplotlib
# scipy

push!(pyimport("sys")["path"], dirname(@__FILE__))
nl = pyimport("network_layout")


function partition(clusters, ext2int, bus_data, branch_data)
    """ Partition network into areas based on coordinates and Scikit method 'method'. """

    bus_ = zeros(length(bus_data), 2)
    for bus in bus_data
        bus_[ext2int[bus.second["bus_i"]], :] = [bus.second["bus_i"] bus.second["vm"]]
    end

    br_ = zeros(length(branch_data), 3)
    for br in branch_data
        br_[br.second["index"], :] = [br.second["f_bus"] br.second["t_bus"] br.second["br_x"]]
    end

    xy = nl.network_map(bus_[:, 1], bus_[:, 2], br_[:, 1], br_[:, 2], br_[:, 3])

    clustering = kmeans(xy', clusters; maxiter=100)
    display(scatter(xy[:,1], xy[:,2], marker_z=clustering.assignments, color=:lightrainbow, legend=false))

    areas = Dict()
    for area in 1:clusters
        areas[area] = findall(x -> x == area, clustering.assignments)
    end
areas
end
