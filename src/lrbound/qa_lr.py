from math import factorial
import copy
import networkx as nx
import numpy as np
import qutip as qt


numInt = [
    1,
    1,
    9,
    17,
    345,
    897,
    30177,
    97185,
    4719249,
    17854785,
    1156537305,
    4979241009,
    408846429225,
    1960817257665,
    196942827937905,
    1036262358096705,
    124004783953852065,
    707693614922363265,
    98874185992337058345,
    606583437653875935825,
]


def nestedInt(r, T):
    return numInt[r - 1] * T**r / factorial(2 * r)


def lr_alpha(r, d, T, alpha):
    k = r / 2
    t1 = (2 / alpha) ** k * (2 * (d - 1) ** k) * nestedInt(r, T)
    t2 = (2 / alpha) ** (k + 1) * (2 * (d - 1) ** k) * nestedInt(r + 1, T)
    t3 = (
        (2 / alpha) ** (k + 1)
        * (2 * (d - 1) ** k * (d * (k + 1) + 1))
        * nestedInt(r + 2, T)
    )
    t4 = (
        (2 / alpha) ** (k + 2)
        * (2 * (d - 1) ** k * (d * (k + 2)))
        * nestedInt(r + 3, T)
    )
    t5 = (
        (2 / alpha) ** (k + 2)
        * (
            2
            * (d - 1) ** k
            * (d**2 * ((k + 1) ** 2 + 3) / 2 + d * (3 * k + 2) / 2 + k)
        )
        * nestedInt(r + 4, T)
    )
    t6 = (
        (2 / alpha) ** (k + 3)
        * (
            2
            * (d - 1) ** k
            * (d**2 * ((k + 2) ** 2 + 3) / 2 + d * (k + 1) / 2 + k - 1)
        )
        * nestedInt(r + 5, T)
    )
    t7 = (4 / alpha * d) ** (k + 3) * nestedInt(r + 6, T)
    return t1 + t2 + t3 + t4 + t5 + t6 + t7


def commutativity_graph(g):
    com_g = copy.deepcopy(g)
    n = g.number_of_nodes()
    m = g.number_of_edges()

    old_edges = g.edges()

    new_nodes = [i for i in range(n, n + m)]
    com_g.add_nodes_from(new_nodes)
    i = 0
    for a, b in old_edges:
        com_g.remove_edge(a, b)
        com_g.add_edges_from([(a, new_nodes[i]), (new_nodes[i], b)])
        i += 1
    return n, com_g


def num_path(g, com_g, n, k):
    l_paths_1 = [
        p
        for p in nx.all_simple_paths(
            com_g, n, [x for x in g.nodes() if g.degree(x) == 2], cutoff=2 * k
        )
    ]
    l_paths_1 = [p for p in l_paths_1 if len(p) > 2 * k - 2]
    l_paths_2 = [
        p
        for p in nx.all_simple_paths(
            com_g, n, [x for x in g.nodes() if g.degree(x) == 1], cutoff=2 * k
        )
    ]
    l_paths_2 = [p for p in l_paths_2 if len(p) > 2 * k - 2]
    return len(l_paths_1) + 2 * len(l_paths_2)


def loc_lr_p3(g, T, alpha):
    r = 8
    k = r / 2
    n, com_g = commutativity_graph(g)
    npg = num_path(g, com_g, n, k)
    npg1 = num_path(g, com_g, n, k + 1)
    npg2 = num_path(g, com_g, n, k + 2)

    t1 = (2 / alpha) ** k * npg * nestedInt(r, T)
    t2 = (2 / alpha) ** (k + 1) * npg * nestedInt(r + 1, T)
    t3 = (2 / alpha) ** (k + 1) * (npg * (2 + 14) + npg1) * nestedInt(r + 2, T)
    t4 = (2 / alpha) ** (k + 2) * (npg * (2 + 16) + npg1) * nestedInt(r + 3, T)
    t5 = (
        (2 / alpha) ** (k + 2)
        * (npg * (2**2 + 2 * 17 + 113) + npg2)
        * nestedInt(r + 4, T)
    )
    t6 = (
        (2 / alpha) ** (k + 3)
        * (npg * (2**2 + 2 * 19 + 144) + npg2)
        * nestedInt(r + 5, T)
    )
    t7 = (4 * 3 / alpha) ** (k + 3) * nestedInt(r + 6, T)
    return t1 + t2 + t3 + t4 + t5 + t6 + t7


def HB(n, alpha):
    """
    Create H0=-sum(sigma_x)/alpha
    """
    h = np.zeros((1 << n, 1 << n))
    for i in range(1 << n):
        bs = bin(i)[2:].zfill(n)
        for b in range(n):
            j_bs = list(bs)
            j_bs[b] = "0"
            if bs[b] == "0":
                j_bs[b] = "1"
            j_bs = "".join(j_bs)
            j = int(j_bs, 2)

            h[i, j] = -1
    return qt.Qobj(h / alpha)


def Hmc(graph):
    """
    Create H1=-sum(1-sigma_z sigma_z)/2
    """
    n = graph.order()
    h = np.zeros((1 << n, 1 << n))
    for i in range(1 << n):
        bs = bin(i)[2:].zfill(n)
        for a, b in graph.edges():
            if bs[int(a)] != bs[int(b)]:
                h[i, i] += -1.00000

    return qt.Qobj(h)


def H0_coeff(t, args):
    return 1 - t / args["tmax"]


def H1_coeff(t, args):
    return t / args["tmax"]


def gen_schedule(graph, alpha):
    """
    Generate the time dependent Hamiltonian
    """
    H0 = HB(graph.order(), alpha)
    H1 = Hmc(graph)
    return [[H0, H0_coeff], [H1, H1_coeff]]


def run_annealing(graph, tmax, alpha):
    n = graph.order()
    state = sum(qt.basis(1 << n, i) for i in range(1 << n))
    psi0 = state.unit()
    drive = gen_schedule(graph, alpha)

    return qt.sesolve(
        drive, psi0, [0.0, tmax], args={"tmax": tmax}, options=qt.Options(nsteps=50000)
    )


def edge_energy(graph, tmax, alpha, edge=("0", "1")):
    obs = graph.copy()
    obs.remove_edges_from(graph.edges())
    obs.add_edge(*edge)
    O_X = -Hmc(obs)
    res = run_annealing(graph, tmax, alpha)
    return res.states[-1].dag().data * O_X.full() * res.states[-1].data
