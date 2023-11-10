import numpy as np
import networkx as nx
import copy
from itertools import combinations
from networkx.algorithms.isomorphism import GraphMatcher


class MyGraphMatcher(GraphMatcher):
    def is_isomorphic(self):
        """Returns True if G1 and G2 are isomorphic graphs."""

        try:
            x = next(self.isomorphisms_iter())
            return True
        except StopIteration:
            return False


def is_isomorphic_same_degrees(g1, g2, edge_match):
    return MyGraphMatcher(g1, g2, edge_match=edge_match).is_isomorphic()


def checknadd(bs, nb):
    em = nx.algorithms.isomorphism.categorical_edge_match(
        "marked", nx.get_edge_attributes(nb, "marked")
    )
    for b in bs:
        if is_isomorphic_same_degrees(nb, b, edge_match=em):
            return 0
    bs.append(nb)


def hkeys(G):
    return tuple(
        sorted(
            [
                round(x, ndigits=6)
                for x in nx.eigenvector_centrality(G, max_iter=200).values()
            ]
        )
    )


def save_ball(balls, p, j):
    balls = [b for lb in balls for b in lb]
    for i, b in enumerate(balls):
        nx.write_adjlist(b, f"balls/{p}/p{p}ball{j+1}_{i+1}")


def update(b, nbs):
    ld = hkeys(b)
    l = nbs.get(ld, [])
    checknadd(l, b)
    nbs[ld] = l


def initialization(all_balls, j):
    ori_ball = all_balls[j]
    balls = {hkeys(ori_ball): [ori_ball]}
    node_to_comp = np.where(np.array(list(dict(nx.degree(ori_ball)).values())) < 3)[0]
    return balls, node_to_comp


def complete_one_node(node, balls):
    new_balls = {}
    for ball_bucket in balls.values():
        for ball in ball_bucket:
            if nx.degree(ball)[node] == 3:
                update(ball, new_balls)
                continue

            if nx.is_regular(ball) and nx.degree(ball, node) == 3:
                update(ball, new_balls)
                continue

            new_neigh = set(
                x for x in ball.nodes if x != node and nx.degree(ball)[x] < 3
            )
            new_neigh.discard(0)
            new_neigh.discard(1)
            if 3 - ball.degree[node] == 1:
                for n1 in [*new_neigh, -1]:
                    new_ball = copy.deepcopy(ball)
                    if n1 == -1:
                        new_ball.add_edge(node, len(new_ball))
                    else:
                        new_ball.add_edge(node, n1)

                    update(new_ball, new_balls)

            if 3 - ball.degree[node] == 2:
                for ns in set(combinations([*new_neigh, -1, -1], 2)):
                    new_ball = copy.deepcopy(ball)
                    for n1 in ns:
                        if n1 == -1:
                            new_ball.add_edge(node, len(new_ball))
                        else:
                            new_ball.add_edge(node, n1)

                    update(new_ball, new_balls)

    return new_balls
