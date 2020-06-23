"""
     makePlots(original_net, reduced_net) -> p::Plot

Create plots for the original and the reduced network


"""

function makePlots(original_net::Dict{String,Any}, reduced_net::Dict{String,Any};
    width=500,
    scale_factor=8,
    scale_factor_EB=1.4,
    scale_factor_red=8,
    scale_factor_EB_red=1.4,
    font_size=10,
    )

    # Prepare Network layout
    # original network
    areas = reduced_net["areas"].clusters

    border_bus = vcat([
        bs.second["indE"] for bs in reduced_net["areas"].data if bs.first > 0
    ]...)
    plt, xy_full = plot_grid(original_net, areas, border_bus; scale_factor=scale_factor, scale_factor_EB=scale_factor_EB, linealpha=0.5)

    # plot reduced network
    bus_key = reduced_net["areas"].data[0]["mapping"]
    heavy = PyDict(Dict(
        b => [x y]
        for
        (b, x, y) in zip(
            bus_key[border_bus],
            xy_full[border_bus, 1],
            xy_full[border_bus, 2],
        )
    ))
    plt_red, xy_full_red = plot_grid(reduced_net, areas, border_bus, bus_key, heavy; scale_factor=scale_factor_red, scale_factor_EB=scale_factor_EB_red, linealpha=0.1)

    phi = MathConstants.golden

    # combine plots
    p = plot(
        plt,
        plt_red,
        layout = grid(1, 2),
        size = (width * phi, width),
        title = ["Original Network" "Reduced Network"],
        titleloc = :center,
        titlefont = font(font_size),
        grid = false
    )
    xaxis!(p, false)
    yaxis!(p, false)
    p
end

function plot_grid(network::Dict{String, Any}, areas::Dict{Int64, Array{Int64,1}}, border_bus::Array{Int64,1};
    scale_factor=scale_factor,
    scale_factor_EB=scale_factor_EB,
    linealpha=linealpha,
    )

    bus = network["bus"]
    branch = network["branch"]
    ext2int = Dict(
        bus["index"] => i
        for
        (i, (k, bus)) in
        enumerate(sort(network["bus"], by = x -> parse(Int, x)))
    )
    bus_ = zeros(length(bus), 2)
    for bus in bus
        bus_[ext2int[bus.second["bus_i"]], :] =
            [bus.second["bus_i"] bus.second["vm"]]
    end

    br_ = zeros(length(branch), 3)
    for br in branch
        br_[br.second["index"], :] =
            [br.second["f_bus"] br.second["t_bus"] br.second["br_x"]]
    end

    # println("Calculating network map of original network...")
    xy_full = ReduceModel.net_layout.network_map(
        bus_[:, 1],
        bus_[:, 2],
        br_[:, 1],
        br_[:, 2],
        br_[:, 3],
    )

    # set colors
    colors = Plots.palette(:lightrainbow, max(length(areas), 2))
    plt = scatter()

    # add lines
    for br in branch
        xl =
            [xy_full[ext2int[br.second["f_bus"]], 1] xy_full[
                ext2int[br.second["t_bus"]],
                1,
            ]]'
        yl =
            [xy_full[ext2int[br.second["f_bus"]], 2] xy_full[
                ext2int[br.second["t_bus"]],
                2,
            ]]'
        plot!(plt, xl, yl, color = "black", linealpha = linealpha, legend = false)
    end

    graph_area = abs(minimum(xy_full)) * abs(maximum(xy_full))
    marker_size = graph_area / (scale_factor * length(bus))
    # border buses
    for area in areas
        # scatter each area
        # add non_border buses
        non_border = [bs for bs in area.second if ~(bs in border_bus)]
        scatter!(
            xy_full[:, 1][non_border],
            xy_full[:, 2][non_border],
            marker = (:circle, marker_size, 0.7),
            color = colors[area.first],
            legend = false,
        )
        # add border buses
        border = [bs for bs in area.second if (bs in border_bus)]
        scatter!(
            plt,
            xy_full[:, 1][border],
            xy_full[:, 2][border],
            marker = (:diamond, marker_size * scale_factor_EB, 0.7),
            color = colors[area.first],
            legend = false,
        )
    end
    plt, xy_full
end

function plot_grid(network::Dict{String, Any}, areas::Dict{Int64, Array{Int64, 1}}, border_bus::Array{Int64, 1}, bus_key::Array{Int64,1}, heavy::PyDict;
    scale_factor=scale_factor,
    scale_factor_EB=scale_factor_EB,
    linealpha=linealpla,
    )

    # set colors
    colors = Plots.palette(:lightrainbow, max(length(areas), 2))

    bus_ = zeros(length(network["bus"]), 2)
    for bus in network["bus"]
        bus_[bus.second["bus_i"], :] =
            [bus.second["bus_i"] bus.second["vm"]]
    end

    br_ = zeros(length(network["branch"]), 3)
    for br in network["branch"]
        br_[br.second["index"], :] =
            [br.second["f_bus"] br.second["t_bus"] br.second["br_x"]]
    end

    xy_full = ReduceModel.net_layout.network_map(
        bus_[:, 1],
        bus_[:, 2],
        br_[:, 1],
        br_[:, 2],
        br_[:, 3],
        heavy = heavy,
    )

    plt = scatter()
    # add lines
    for br in network["branch"]
        xl =
            [xy_full[br.second["f_bus"], 1] xy_full[
                br.second["t_bus"],
                1,
            ]]'
        yl =
            [xy_full[br.second["f_bus"], 2] xy_full[
                br.second["t_bus"],
                2,
            ]]'
        plot!(plt, xl, yl, color = "black", linealpha = linealpha, legend = false)
    end

    gr_area = abs(minimum(xy_full)) * abs(maximum(xy_full))
    marker_size = gr_area / (scale_factor * length(network["bus"]))

    areas_reduced = Dict(
        [area => [
            b[2]["bus_i"] for b in network["bus"] if b[2]["zone"] == area
        ] for area = 1:length(areas)],
    )
    border_bus = bus_key[border_bus]
    # border buses
    for area in areas_reduced
        # scatter each area
        # add non_border buses
        non_border = [bs for bs in area.second if ~(bs in border_bus)]
        scatter!(
            plt,
            xy_full[non_border, 1],
            xy_full[non_border, 2],
            marker = (:circle, marker_size, 0.7),
            color = colors[area.first],
            legend = false,
        )
        # add border buses
        border = [bs for bs in area.second if (bs in border_bus)]
        scatter!(
            plt,
            xy_full[border, 1],
            xy_full[border, 2],
            marker = (:diamond, marker_size * scale_factor_EB, 0.7),
            color = colors[area.first],
            legend = false,
        )
    end
    plt, xy_full
end
