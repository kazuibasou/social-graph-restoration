import graph
import numpy as np
from collections import defaultdict

class SampledData():
    def __init__(self, v):
        self.index = v
        self.nlist = []

def query(G: graph.Graph, v):
    data = SampledData(v)
    data.nlist = list(G.nlist[v])

    return data

def select_seed(G: graph.Graph):

    return np.random.choice(G.nodes)

def random_walk(G, sample_size, seed):

    if sample_size <= 0 or sample_size > G.N:
        print("Error: The number of nodes to be queried must be not less than 0 and not more than " + str(G.N) + ".")
        exit()

    v = seed
    samplinglist = []
    queried = defaultdict(int)

    while len(queried) < sample_size:
        data = query(G, v)
        queried[v] = 1
        samplinglist.append(data)
        v = np.random.choice(data.nlist)

    return samplinglist