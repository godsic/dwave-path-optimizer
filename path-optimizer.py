#!/usr/bin/env python3
"""
Copyright (c) 2020 Mykola Dvornik

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

from dwave.system import LeapHybridSampler
from neal import SimulatedAnnealingSampler
from collections import defaultdict
import dimod
import json
import numpy as np
import networkx as nx
import matplotlib.pylab as plt
import argparse

plt.switch_backend("Qt5Agg")

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--L_max", metavar="L_max", type=float,
                    help="deliver range", default=500.)
parser.add_argument("--p_0", metavar="p0", type=int, nargs=2,
                    help="coordinates of origin", default=[70, -80])
parser.add_argument("--l_g", metavar="l_g", type=float,
                    help="weight of the goal function", default=0.005)
parser.add_argument("--l_I", metavar="l_I", type=float,
                    help="weight of number of taken segmets maximizer", default=0.002)
parser.add_argument("--l_II", metavar="l_II", type=float,
                    help="weight of discontinuous routes penalties", default=0.)
parser.add_argument("--l_III", metavar="l_III", type=float,
                    help="weight of two-segmets-per-target constrain", default=0.01)
parser.add_argument("--l_IV", metavar="l_IV", type=float,
                    help="weight of starting point constrain", default=0.1)
parser.add_argument("--l_V", metavar="l_V", type=float,
                    help="weight of discontinuous routes penalties", default=0.1)
parser.add_argument("--simulate", action='store_true',
                    help="toggle simulated annealing")
parser.add_argument("--debug", action='store_true',
                    help="toggle verbose visualization")
parser.add_argument('coordinates_file', nargs='?', type=argparse.FileType(
    'r'), help="2D coordinates filename", default="data/coordinates.dat")

args = parser.parse_args()

# Travel range
L_max = args.L_max

# starting point
p0 = 0
p0_x = args.p_0[0]
p0_y = args.p_0[1]

# Lagrangian multipliers for the constrains
l_g = args.l_g
l_I = args.l_I
l_II = args.l_II
l_III = args.l_III
l_IV = args.l_IV
l_V = args.l_V

# Other
debug = args.debug

data = np.loadtxt(args.coordinates_file)
p_x = data[:, 0]
p_y = data[:, 1]

p_x = np.insert(p_x, 0, p0_x)
p_y = np.insert(p_y, 0, p0_y)
# number of targets
N_p = len(p_x)

max_p_x = 1.1 * abs(max(p_x.min(), p_x.max(), key=abs))
max_p_y = 1.1 * abs(max(p_y.min(), p_y.max(), key=abs))

N_segments = int(N_p * (N_p - 1) / 2)

if debug:
    print("p_x: {}".format(p_x))
    print("p_y: {}".format(p_y))

print("Number of targets: {}".format(N_p))
print("Number of segmets: {}".format(N_segments))

G = nx.Graph()

l = np.zeros(N_segments)
for i in range(N_p):
    idx0 = i * N_p - np.sum(range(i + 1), dtype=int)
    G.add_node(i)
    for j in range(i + 1, N_p):
        idx = idx0 + (j - i - 1)
        d = np.sqrt((p_x[i]-p_x[j])**2 + (p_y[i]-p_y[j])**2) / L_max
        l[idx] = d
        G.add_node(j)
        G.add_edge(i, j, weight=d, n=idx)
        print("{} = ({},{}): {}".format(idx, i, j, d)) if debug else True

print(l) if debug else True

all_edges = np.array([e for e in G.edges])

# Q = Q1 + λ·Q2
def QUBO_sum(Q1, λ, Q2):
    return {k: Q1.get(k, 0.) + λ * Q2.get(k, 0.) for k in set(Q1) | set(Q2)}


if __name__ == "__main__":

    # Goal: find path that does not exceed working range
    Q_goal = defaultdict(float)
    for i in range(N_segments):
        Q_goal[(i, i)] = l[i] * (l[i] - 2.0)
        for j in range(i+1, N_segments):
            Q_goal[(i, j)] = 2.0 * l[i] * l[j]

    # Constrain I: take as many segments as it possible
    Q_cI = defaultdict(float)
    for i in range(N_segments):
        Q_cI[(i, i)] = -1.

    # Constrain II: take only connected segments
    Q_cII_ = defaultdict(float)
    for i in range(N_segments):
        e_i = all_edges[i]
        Q_cII_[(i, i)] = -1.
        for j in range(i+1, N_segments):
            e_j = all_edges[j]
            D_ij = 0. if set(e_i).isdisjoint(e_j) else 1.
            Q_cII_[(i, j)] = (2. - 1. * D_ij)
            for p in range(j+1, N_segments):
                e_p = all_edges[p]
                D_jp = 0. if set(e_j).isdisjoint(e_p) else 1.
                D_ip = 0. if set(e_i).isdisjoint(e_p) else 1.
                Q_cII_[(i, j, p)] = 2. * (D_ij * D_jp + D_ij *
                                          D_ip + D_jp * D_ip - D_ij - D_jp - D_ip)
                for q in range(p+1, N_segments):
                    e_q = all_edges[q]
                    D_pq = 0. if set(e_p).isdisjoint(e_q) else 1.
                    D_jq = 0. if set(e_j).isdisjoint(e_q) else 1.
                    D_iq = 0. if set(e_i).isdisjoint(e_q) else 1.
                    Q_cII_[(i, j, p, q)] = 2. * \
                        (D_ij * D_pq + D_ip * D_jq + D_iq * D_jp)

    bqm = dimod.make_quadratic(Q_cII_, 20.0, dimod.BINARY)
    Q_cII = bqm.to_qubo()[0]

    # Constrain III: minimize the amount of disjoint segments
    Q_cIII_ = [defaultdict(float) for p in range(N_p)]
    targets = [*(range(N_p))]
    targets.remove(p0)

    for t in targets: 
        for ei in G.edges(t):
            i = G.get_edge_data(*ei)['n']
            Q_cIII_[t][(i, i)] = (1. - 2. * 2.)
            for ej in G.edges(t):
                j = G.get_edge_data(*ej)['n']
                if j > i:
                    Q_cIII_[t][(i, j)] = 2.

    # Constrain IV: start from a giving point
    Q_cIV = defaultdict(float)
    for ei in G.edges(p0):
        i = G.get_edge_data(*ei)['n']
        Q_cIV[(i, i)] = -1.
        for ej in G.edges(p0):
            j = G.get_edge_data(*ej)['n']
            if j > i:
                Q_cIV[(i, j)] = 2.

    # Constrain V: penaltize discontinious trips
    Q_cV = defaultdict(float)
    for i in range(N_segments):
        e_i = all_edges[i]
        for j in range(i+1, N_segments):
            e_j = all_edges[j]
            D_ij = 0. if set(e_i).isdisjoint(e_j) else 1.
            Q_cV[(i, j)] = (1.0 - D_ij)

    Q = defaultdict(float)

    Q = QUBO_sum(Q, l_V, Q_cV)
    Q = QUBO_sum(Q, l_IV, Q_cIV)

    for Q_cIII in Q_cIII_:
        Q = QUBO_sum(Q, l_III, Q_cIII)

    Q = QUBO_sum(Q, l_II, Q_cII)
    Q = QUBO_sum(Q, l_I, Q_cI)
    Q = QUBO_sum(Q, l_g, Q_goal)

    if args.simulate:
        solution = SimulatedAnnealingSampler().sample_qubo(Q, num_reads=1000)
    else:
        solution = LeapHybridSampler().sample_qubo(Q)

    print(solution.first) if debug else True

    result = np.array([solution.first.sample[key]
                       for key in range(N_segments)])
    print(result) if debug else True

    traveled_distance = np.sum(result * l) * L_max

    print("Took {} segments with total travel distance of {}".format(
        int(np.sum(result)), traveled_distance))

    visited_edges = all_edges[np.where(result == 1.0)]

    edges, weights = zip(*nx.get_edge_attributes(G, 'weight').items())

    pos = {key: [p_x[key], p_y[key]] for key in range(N_p)}
    print(pos) if debug else True

    nx.draw(G, pos, node_color='b', edgelist=edges,
            edge_color=weights, width=2*debug, edge_cmap=plt.cm.Greys)

    if debug:
        edge_labels = nx.get_edge_attributes(G, "n")
        nx.draw_networkx_labels(G, pos, font_color='w')
        nx.draw_networkx_edge_labels(
            G, pos, font_color='r', edge_labels=edge_labels)

    nx.draw_networkx_edges(G, pos, edgelist=visited_edges,
                           edge_color='g', alpha=0.5, width=4.0)

    plt.axhline(y=0, ls='--', c='orange', lw=0.5)
    plt.axvline(x=0, ls='--', c='orange', lw=0.5)

    plt.xlim(-max_p_x, max_p_x)
    plt.ylim(-max_p_y, max_p_y)

    plt.text(p0_x, p0_y, "▲", fontsize=12, c='w', ha='center', va='center')
    plt.show()
