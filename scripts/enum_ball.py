import os
import networkx as nx
from tqdm import tqdm
from lrbound.enumeration import *

BALLS_DIRECTORY = "./balls"
# Creating directories to store all the enumerated balls
if not os.path.isdir(BALLS_DIRECTORY):
    os.mkdir(BALLS_DIRECTORY)
    for p in range(4):
        os.mkdir(os.path.join(BALLS_DIRECTORY, f"{p}"))


# Generating the only ball of radius 0
edge = nx.Graph()
edge.add_edge(0, 1)
nx.set_edge_attributes(edge, {n: n == (0, 1) for n in edge.edges()}, "marked")

all_balls = [edge]
new_balls = []
for p in range(1, 4):
    print(">>>>   Enumerating all balls of radius p =", p)
    for j in tqdm(range(len(all_balls))):
        balls, node_to_comp = initialization(all_balls, j)
        for node in node_to_comp:
            balls = complete_one_node(node, balls)
        save_ball(balls.values(), p, j)
        new_balls.extend(ball for bucket in balls.values() for ball in bucket)
    all_balls = new_balls
