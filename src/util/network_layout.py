# -*- coding: utf-8 -*-
"""
network_layout by Jon Olauson (olauson@gmail.com), version 2018-09-12

Makes a layout of an electric grid
1. Calculate distance matrix
2. Project in two dimensions
3. Possibly make a simple plot
4. Return bus coordinates (and possibly stress and computational time)

For more information and citation:
J.Olauson, M.Marin, L.Söder
"Creating power system network layouts: A fast algorithm suitable for Python and Matlab"
submitted to IEEE transactions on power systems

Note that time.clock() is used for timing
In this implementation, the network needs to be fully connected
Tested with Python 2.7 and 3.6
See network_layout_examples.py for some examples

Main function is network_map() which call other functions:

Required input (lists or numpy arrays)
  bus - bus numbers
  v_base - bus base voltages in kV
  bus1 - branch from buses
  bus2 - branch to buses
  X - branch reactances in p.u.

Optional input
  heavy - {bus:[x, y [,heaviness]], ...} for predefined coordinates (default None)
          heaviness = inf -> not moved (default), >1 -> moved less than other buses
          (see Example 2)
  ex - exponent for X in distance calculation (0.3)
  ev - exponent for V in distance calculation (0.3)
  iterations - iterations in MDS (10)
  plot - make simple plot (True)
  print_results - print time for calculations and system stress (False)
  return_stress_time - if True, also return stress an computational time (False)
  bus_size - size of bus with highest voltage when plotting (5)
  min_d - minimum allowed desired distance as share of system diameter (0.001)
  w_exp - exponent for branch weight;  w = d ^ w_exp (-2 default, -1 for Sammon)
  norm_stress - normalise stress with sum of desired distances (False)
  epsilon - used to calculate smallest step size (1.0)

Output
  N-by-2 matrix with x and y coordinates of buses
  Stress, time_distance, time_layout (if return_stress_time=True)

"""

import numpy as np
import matplotlib.pyplot as plt
from numpy import flatnonzero as find
from numpy import atleast_1d as arr
from scipy.spatial.distance import cdist
from time import clock
from scipy.sparse.csgraph import dijkstra
from scipy.optimize import minimize


def network_map(bus, v_base, bus1, bus2, X, heavy=None, return_stress_time=False,
                ex=0.3, ev=0.3, iterations=10, plot=False, print_results=False,
                bus_size=5, min_d=0.001, w_exp=-2, norm_stress=False, epsilon=1.):
    """ Main function: calculates distance matrix, project in 2D and plot. """

    bus1, bus2, X, V = aggregate_branches(bus, v_base, bus1, bus2, X)  # aggregate parallel branches
    d, t1 = calculate_distance(bus, bus1, bus2, X, V, ex, ev)  # desired distances (N*N matrix)
    xy, t2, S = VSGD(d, heavy, bus, iterations, min_d, w_exp, epsilon, norm_stress)  # calculate layout

    if print_results:
        print('Time for distance matrix: %0.2f sec' % t1)
        print('Time for layout calculation: %0.2f sec' % t2)
        print('Stress: %0.3e' % S)

    if plot:  # Simple plot
        plot_grid(xy, bus, v_base, bus1, bus2, bus_size)

    if return_stress_time:
        return xy, S, t1, t2
    else:
        return xy


def mult_ind(a, b):
    """ Get indices for elements of a in b, returns numpy array."""
    bind = {}
    for i, elt in enumerate(b):
        if elt not in bind:
            bind[elt] = i
    return arr([bind[itm] for itm in a])


def aggregate_branches(bus, v_base, bus1, bus2, x0):
    """ Aggregate parallel branches. """

    v1 = v_base[mult_ind(bus1, bus)]  # Base voltage at from bus
    v2 = v_base[mult_ind(bus2, bus)]  # Base voltage at to bus

    pairs = np.sort(np.column_stack((bus1, bus2)))  # bus pairs (lowest number first)
    unique_pairs, ib = np.unique(pairs, axis=0, return_inverse=True)
    bus1, bus2 = unique_pairs[:, 0], unique_pairs[:, 1]
    x = 0. * bus1  # reactance (parallel coupled)
    v = 0. * bus1  # minimum voltage of bus pair
    for n in range(len(bus1)):
        ind = find(ib == n)
        v[n] = np.minimum(v1[ind[0]], v2[ind[0]])
        temp = x0[ind]
        if sum(temp == 0) > 0:
            x[n] = 0
        else:
            x[n] = sum(1 / temp) ** -1  # parallel coupling

    return bus1, bus2, x, v


def calculate_distance(bus, bus1, bus2, X, V, ex=0.3, ev=0.3):
    """ Calculate distance matrix with Dijkstra's algorithm. """
    ta = clock()

    dist = np.zeros((len(bus), len(bus)))  # branch distances
    ind1 = mult_ind(bus1, bus)
    ind2 = mult_ind(bus2, bus)

    # Handle limitation in djikstra
    # Now minimum distance is 1 except for zero distances which are set to 1e-6
    temp = np.abs(X) ** ex * V ** ev
    temp /= min(temp[temp > 0])
    temp[temp == 0] = 1e-6
    dist[ind1, ind2] = dist[ind2, ind1] = temp

    d = dijkstra(dist, directed=False)  # shortest path between all bus-pairs
    # quick fix -- TO REVIEW!
    max_val = np.max(d)
    if np.isinf(max_val):
        max_val = np.unique(d)[-2]

    d = d.astype(np.float32) / max_val * 100  # normalise so system diameter = 100 units
    t = clock() - ta
    return d, t


def create_sets(N):
    """ Make sets of bus pairs (indices)."""
    sets = []
    for n in range(1, N):
        pairs = np.zeros((N - n, 2), int)  # pairs on diagonal n
        pairs[:, 0] = np.arange(N - n)
        pairs[:, 1] = pairs[:, 0] + n
        mask = np.mod(range(N - n), 2 * n) < n
        s1 = pairs[mask]
        s2 = pairs[~mask]
        if len(s1) > 0:
            sets.append(s1)
        if len(s2) > 0:
            sets.append(s2)
    return sets


def heavy_data(heavy, bus):
    """ Extract indices, coordinates and weights from dict with heavy buses."""
    h_ind = mult_ind(list(heavy.keys()), bus)  # indeces for buses
    h_xy = np.zeros((len(heavy), 2))  # pre-defined coordinates
    h_weight = np.zeros(len(heavy))  # weight (>1 -> moves shorter, np.inf (default) -> no movement)
    for n, (k, v) in enumerate(heavy.items()):
        h_xy[n, :] = v[0:2]
        if len(v) == 3:
            h_weight[n] = v[2]
        else:
            h_weight[n] = np.inf
    return h_ind, h_xy, h_weight


def xy_rescale(X, d):
    """ Find scaling and shift in x and y compatible with predefined coordinates. """

    def OF(s):
        d2 = ((s[0] * (x1 - x2)) ** 2 + (s[1] * (y1 - y2)) ** 2) ** 0.5
        return np.sum(np.abs(d - d2))

    N = len(X)
    if N == 1:  # no rescaling possible
        scale = [1, 1]
    elif N == 2:  # same scaling in x and y
        scale = [d[0, 1] / cdist(X, X)[0, 1]] * 2
    else:  # optimise scaling in x and y
        mask = np.ones((N, N)) == 1 - np.tril(np.ones((N, N)))
        d = d[mask]  # upper triangular as vector
        ind2, ind1 = np.meshgrid(range(N), range(N))
        x1 = X[ind1[mask], 0]
        y1 = X[ind1[mask], 1]
        x2 = X[ind2[mask], 0]
        y2 = X[ind2[mask], 1]
        scale = minimize(OF, [1, 1], method='Nelder-Mead', options={'maxiter': 100}).x
    shift = np.mean(X, axis=0)  # move center of random layout to center of heavy buses

    return scale[0], scale[1], shift[0], shift[1]


def VSGD(d, heavy, bus, iterations, min_d, w_exp, epsilon, norm_stress):
    """ Find layout with multidemensional scaling using vectorised stochastic gradient descent. """
    ta = clock()
    d[d <= np.max(d) * min_d] = np.max(d) * min_d  # not too small distances
    N = len(d)
    mask = np.ones((N, N)) == 1 - np.tril(np.ones((N, N)))  # upper triangular except diagonal
    X = np.random.rand(N, 2).astype(np.float32) * 100 - 50  # random layout with diameter 100

    # Some operations to account for predefined coordinates
    if heavy is not None:
        h_ind, h_xy, h_weight = heavy_data(heavy, bus)
        xs, ys, dx, dy = xy_rescale(h_xy, d[h_ind, :][:, h_ind])  # scale and shift
        X[:, 0] += dx  # center of initial layout moved
        X[:, 1] += dy
        X[h_ind, :] = h_xy
        bw = np.ones(len(bus))  # bus weight
        bw[h_ind] = h_weight

    w = d ** w_exp  # bus-pair weights (lower for distant buses)

    stepmax = 1 / np.min(w[mask])
    stepmin = epsilon / np.max(w[mask])
    lambda1 = -np.log(stepmin / stepmax) / (iterations - 1)  # exponential decay of allowed adjustment
    sets = create_sets(N)  # construct sets of bus pairs
    for iteration in range(iterations):

        step = stepmax * np.exp(-lambda1 * iteration)  # how big adjustments are allowed?
        rand_order = np.random.permutation(N)  # we don't want to use the same pair order each iteration

        for p in sets:
            b1, b2 = rand_order[p[:, 0]], rand_order[p[:, 1]]  # arrays of bus1 and bus2

            # current distances (possibly accounting for rescaling of system)
            if heavy is not None:
                d2 = ((xs * (X[b1, 0] - X[b2, 0])) ** 2 + (ys * (X[b1, 1] - X[b2, 1])) ** 2) ** 0.5
            else:
                d2 = ((X[b1, 0] - X[b2, 0]) ** 2 + (X[b1, 1] - X[b2, 1]) ** 2) ** 0.5
            r = (d[b1, b2] - d2)[:, None] / 2 * (X[b1] - X[b2]) / d2[:, None]  # desired change
            dX1 = r * np.minimum(1, w[b1, b2] * step)[:, None]
            dX2 = -dX1
            if heavy is not None:
                dX1 /= bw[b1, None]  # divide change with heaviness
                dX2 /= bw[b2, None]
            X[b1, :] += dX1  # update position
            X[b2, :] += dX2

    if heavy is None:
        S = sum((w * (d - cdist(X, X)) ** 2)[mask])  # stress
    else:  # if predefined coordinates provided
        S = sum((w * (d - cdist(X * [xs, ys], X * [xs, ys])) ** 2)[mask])
    if norm_stress:
        S /= sum(d[mask])
    t2 = clock() - ta

    return X, t2, S


def plot_grid(xy, bus, v_base, bus1, bus2, bus_size):
    """ Simple plot of resulting layout. """
    fig, ax = plt.subplots()

    # Buses
    ax.scatter(xy[:, 0], xy[:, 1], float(bus_size) * v_base / np.max(v_base))  # size depend on voltage

    # Plot branches on vectorised form (for plotting speed)
    ind1 = mult_ind(bus1, bus)  # indices to get coordinates for "from buses"
    ind2 = mult_ind(bus2, bus)
    xy_branch = np.nan * np.zeros((len(bus1) * 3 - 1, 2))  # line coordinates (separated by nans)
    temp = np.mod(range(len(bus1) * 3 - 1), 3)  # to put things on the right lines
    xy_branch[temp == 0, :] = xy[ind1, :]
    xy_branch[temp == 1, :] = xy[ind2, :]
    ax.plot(xy_branch[:, 0], xy_branch[:, 1], c='r', lw=0.5)

    ax.set_aspect('equal')
    plt.show()
