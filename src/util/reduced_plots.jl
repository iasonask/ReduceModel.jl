"""
     makePlots(original_net, reduced_net) -> p::Plot

Create plots for the original and the reduced network


"""

function makePlots(original_net, reduced_net)

    # Prepare Network layout
    # original network
    bus = original_net["bus"]
    branch = original_net["branch"]
    ext2int = Dict(
        bus["index"] => i
        for
        (i, (k, bus)) in
        enumerate(sort(original_net["bus"], by = x -> parse(Int, x)))
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
    areas = reduced_net["areas"].clusters

    border_bus = vcat([
        bs.second["indE"] for bs in reduced_net["areas"].data if bs.first > 0
    ]...)

    # set colors
    colors = palette(:lightrainbow, length(areas))
    plt = scatter()
    # estimated area of graph
    scale_factor = 8
    scale_factor_EB = 1.4

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
        plot!(plt, xl, yl, color = "black", linealpha = 0.5, legend = false)
    end

    graph_area = abs(minimum(xy_full)) * abs(maximum(xy_full))
    marker_size = graph_area / (scale_factor * length(bus))
    # border buses
    for area in areas
        # scatter each area
        # add non_border buses
        non_border = [bs for bs in area.second if ~(bs in border_bus)]
        # non_border = [bs for bs in area if ~(bs in border_bus)]
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

    # plot reduced network
    bus_key = reduced_net["areas"].data[0]["mapping"]
    # heavy = PyDict(Dict(b - 1 => xy_full[b, :] for b in border_bus))
    heavy = PyDict(Dict(
        b => [x y]
        for
        (b, x, y) in zip(
            bus_key[border_bus],
            xy_full[border_bus, 1],
            xy_full[border_bus, 2],
        )
    ))

    bus_red_ = zeros(length(reduced_net["bus"]), 2)
    for bus in reduced_net["bus"]
        bus_red_[bus.second["bus_i"], :] =
            [bus.second["bus_i"] bus.second["vm"]]
    end

    br_red_ = zeros(length(reduced_net["branch"]), 3)
    for br in reduced_net["branch"]
        br_red_[br.second["index"], :] =
            [br.second["f_bus"] br.second["t_bus"] br.second["br_x"]]
    end

    # println("Calculating network map of reduced network...")
    xy_full_red = ReduceModel.net_layout.network_map(
        bus_red_[:, 1],
        bus_[:, 2],
        br_red_[:, 1],
        br_red_[:, 2],
        br_red_[:, 3],
        heavy = heavy,
    )

    plt1 = scatter()
    # add lines
    for br in reduced_net["branch"]
        xl =
            [xy_full_red[br.second["f_bus"], 1] xy_full_red[
                br.second["t_bus"],
                1,
            ]]'
        yl =
            [xy_full_red[br.second["f_bus"], 2] xy_full_red[
                br.second["t_bus"],
                2,
            ]]'
        plot!(plt1, xl, yl, color = "black", linealpha = 0.1, legend = false)
    end

    gr_area = abs(minimum(xy_full_red)) * abs(maximum(xy_full_red))
    scale_factor_red = 1
    scale_factor_EB_red = 0.5
    marker_size_red = gr_area / (scale_factor_red * length(bus))

    areas_reduced = Dict(
        [area => [
            b[2]["bus_i"] for b in reduced_net["bus"] if b[2]["zone"] == area
        ] for area = 1:length(areas)],
    )
    border_bus_red = bus_key[border_bus]
    # border buses
    for area in areas_reduced
        # scatter each area
        # add non_border buses
        non_border = [bs for bs in area.second if ~(bs in border_bus_red)]
        scatter!(
            plt1,
            xy_full_red[non_border, 1],
            xy_full_red[non_border, 2],
            marker = (:circle, marker_size_red, 0.7),
            color = colors[area.first],
            legend = false,
        )
        # add border buses
        border = [bs for bs in area.second if (bs in border_bus_red)]
        scatter!(
            plt1,
            xy_full_red[border, 1],
            xy_full_red[border, 2],
            marker = (:diamond, marker_size_red * scale_factor_EB_red, 0.7),
            color = colors[area.first],
            legend = false,
        )
    end

    width = 500
    phi = MathConstants.golden

    p = plot(
        plt,
        plt1,
        layout = grid(2, 1),
        size = (width, width * phi),
        title = ["Original Network" "Reduced Network"],
        titleloc = :center,
        titlefont = font(8),
    )
    p
end
