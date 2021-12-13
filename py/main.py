import graph
import sampling
import generation
import time
import os
import sys

def read_graph(graphname):
    datadir = "../data/"

    graphpath = str(datadir) + str(graphname) + ".txt"
    f = open(graphpath, 'r')
    line = f.readline()
    G = graph.Graph()
    space = ' '

    while line:
        data = line[:-1].split(space)
        u = int(data[0])
        v = int(data[1])
        G.nlist[u].append(v)
        G.nlist[v].append(u)
        line = f.readline()

    f.close()

    G.nodes = list(G.nlist.keys())
    G.N = len(G.nodes)

    m = 0
    for v in G.nodes:
        d = int(len(G.nlist[v]))
        m += d
        if d > G.maxd:
            G.maxd = d

    G.M = int(m/2)

    print('Read ' + str(graphname) + " graph.")
    print("Number of nodes: " + str(G.N))
    print("Number of edges: " + str(G.M) + "\n")

    return G

def write_edge_list(graphname,genG):

    if not os.path.exists("../gen_graph/"):
        os.mkdir("../gen_graph/")

    f_path = "../gen_graph/" + str(graphname) + ".txt"
    f = open(f_path, 'w')

    for v in genG.nodes:
        for w in genG.nlist[v]:
            if v <= w:
                f.write(str(v) + " " + str(w) + "\n")

    f.close()

    print("Wrote the generated graph at social-graph-restoration/gen_graph/.\n")

    return 0

if __name__ == '__main__':

    args = sys.argv
    graphname = str(args[1])
    sample_size = int(args[2])

    G = read_graph(graphname)

    seed = sampling.select_seed(G)
    samplinglist = sampling.random_walk(G,sample_size,seed)

    start = time.time()
    genG = generation.graph_restoration_method(samplinglist,False)
    end = time.time()
    print("Generation time [sec]: " + str(end-start))

    write_edge_list(graphname, genG)
