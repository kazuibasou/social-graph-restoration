import networkx as nx
import igraph
from collections import Counter
from collections import defaultdict
from collections import deque
import math
import numpy as np
import warnings
warnings.simplefilter("ignore", np.ComplexWarning)

class Graph():

    def __init__(self):
        self.nodes = set() # set of nodes
        self.qry_nodes = set() # set of queried nodes
        self.vis_nodes = set() # set of visible nodes
        self.nlist = defaultdict(list) # lists of neighbors

        # graph properties
        self.N = 0 # number of nodes
        self.M = 0 # number of edges
        self.maxd = 0 # maximum degree
        self.aved = 0 # average degree
        self.acc = 0 # average clustering coefficient
        self.apl = 0 # average shortest path length
        self.diameter = 0 # diameter
        self.lambda_1 = 0
        self.dd = defaultdict(float) # degree distribution
        self.num_deg = defaultdict(int) # number of nodes with degrees
        self.jdd = defaultdict(lambda: defaultdict(float)) # joint degree distribution
        self.knn = defaultdict(float) # neighbor connectivity
        self.num_tri = defaultdict(int) # number of triangles of nodes
        self.ddcc = defaultdict(float) # degree-dependent clustering coefficient
        self.cnd = defaultdict(float)  # common neighbor distribution
        self.spld = defaultdict(float) # shortest path length distribution
        self.ddbc = defaultdict(float) # degree-dependent betweeness centrality

def convert_to_Graph_of_networkx(G: Graph):
    nxG = nx.Graph()

    for v in G.nodes:
        for w in G.nlist[v]:
            nxG.add_edge(v, w)

    return nxG

def convert_to_MultiGraph_of_networkx(G: Graph):
    nxG = nx.MultiGraph()

    for v in G.nodes:
        for w in G.nlist[v]:
            nxG.add_edge(v, w)

    return nxG

def convert_to_Graph_of_igraph(G: Graph):
    iG = igraph.Graph()
    edges = []

    for v in G.nodes:
        for w in G.nlist[v]:
            if w >= v:
                edges.append([v, w])

    iG.add_vertices(len(G.nodes))
    iG.add_edges(edges)

    return iG

def is_connected(G: Graph):
    Q = deque()
    V = list(G.nodes)
    visit = defaultdict(int)

    v = V[0]
    Q.append(v)
    visit[v] = 1
    n = 0

    while len(Q) > 0:
        v = Q.popleft()
        n += 1

        for w in G.nlist[v]:
            if visit[w] == 0:
                visit[w] = 1
                Q.append(w)

    return n == len(V)

def largest_connected_component(G: Graph):
    search = {v: 0 for v in list(G.nodes)}
    LCC_nodes = []
    LCC_nlist = defaultdict(list)
    n = 0
    N = len(G.nodes)

    for v in G.nodes:
        if sum(list(search.values())) >= N - n:
            break

        if search[v] == 1:
            continue

        Q = deque()
        visit = defaultdict(int)
        visit[v] = 1
        Q.append(v)
        CC_nodes = []
        CC_nlist = defaultdict(list)

        while len(Q) > 0:
            u = Q.popleft()
            CC_nodes.append(u)
            search[u] = 1
            for w in G.nlist[u]:
                CC_nlist[u].append(w)
                if visit[w] == 0:
                    visit[w] = 1
                    Q.append(w)

        if len(CC_nodes) > n:
            n = len(CC_nodes)
            LCC_nodes = list(CC_nodes)
            LCC_nlist = defaultdict(CC_nlist)

    LCC = Graph()
    LCC.nodes = list(LCC_nodes)
    LCC.N = len(LCC.nodes)
    LCC.nlist = defaultdict(LCC_nlist)

    m = 0
    for v in LCC.nodes:
        d = int(len(LCC.nlist[v]))
        m += d
        if d > LCC.maxd:
            LCC.maxd = d
    LCC.M = int(m)/2

    return LCC

def add_edge(G: Graph, u, v):
    G.nlist[u].append(v)
    G.nlist[v].append(u)

    return G

def remove_edge(G: Graph, u, v):
    G.nlist[u].remove(v)
    G.nlist[v].remove(u)

    return G

def calc_dd(G: Graph):
    V = list(G.nodes)
    n = len(V)
    degrees = {v: int(len(G.nlist[v])) for v in V}
    num_deg = Counter(list(degrees.values()))
    G.dd = defaultdict(float)
    for d in num_deg:
        G.dd[d] = float(num_deg[d])/n

    return

def calc_num_deg(G: Graph):
    V = list(G.nodes)
    degrees = {v: int(len(G.nlist[v])) for v in V}
    G.num_deg = Counter(list(degrees.values()))

    return

def calc_jdd(G: Graph):
    G.jdd = defaultdict(lambda: defaultdict(float))
    V = list(G.nodes)

    for v in V:
        k = int(len(G.nlist[v]))
        for w in list(G.nlist[v]):
            l = int(len(G.nlist[w]))
            G.jdd[k][l] += 1

    for k in G.jdd:
        for l in G.jdd[k]:
            G.jdd[k][l] = float(G.jdd[k][l])/(2*G.M)

    return

def calc_knn(G: Graph):
    V_k = defaultdict(list)

    for v in list(G.nodes):
        k = len(G.nlist[v])
        V_k[k].append(v)

    G.knn = defaultdict(float)
    for k in V_k:
        if k*len(V_k[k]) == 0:
            continue
        for v in V_k[k]:
            for w in G.nlist[v]:
                G.knn[k] += len(G.nlist[w])
        G.knn[k] = float(G.knn[k])/(k*len(V_k[k]))

    return

def calc_num_tri(G: Graph):
    G.num_tri = defaultdict(int)
    V = list(G.nodes)

    for v in V:
        d = int(len(G.nlist[v]))
        if d == 0 or d == 1:
            continue

        n_t = 0
        for i in range(0, d-1):
            x = G.nlist[v][i]
            for j in range(i+1, d):
                y = G.nlist[v][j]
                if v != x and x != y and y != v:
                    n_t += G.nlist[x].count(y)

        G.num_tri[d] += n_t

    return

def calc_clustering(G: Graph):
    V_d = defaultdict(int)
    sum_lcc_d = defaultdict(float)
    sum_lcc = 0

    V = list(G.nodes)

    for v in V:
        d = int(len(G.nlist[v]))
        V_d[d] += 1

        if d == 0 or d == 1:
            continue

        lcc = 0
        for i in range(0, d-1):
            x = G.nlist[v][i]
            for j in range(i+1, d):
                y = G.nlist[v][j]
                if v != x and x != y and y != v:
                    lcc += 2*G.nlist[x].count(y)

        lcc = float(lcc)/(d*(d-1))
        sum_lcc_d[d] += lcc
        sum_lcc += lcc

    G.ddcc = defaultdict(float)
    for d in V_d:
        if V_d[d] > 0:
            G.ddcc[d] = float(sum_lcc_d[d])/V_d[d]

    N = len(V)
    G.acc = float(sum_lcc)/N

    return

def calc_common_neighbor_distribution(G:Graph):
    G.cnd = defaultdict(float)

    for i in G.nodes:
        for j in G.nlist[i]:
            if j <= i:
                continue
            m = 0
            for k in G.nlist[i]:
                if k == i and k == j:
                    continue
                m += G.nlist[j].count(k)
            G.cnd[m] += 1

    norm = sum(list(G.cnd.values()))
    for m in G.cnd:
        G.cnd[m] = float(G.cnd[m])/norm

    return

def calc_shortest_path_properties(G:Graph):
    #Note: calculate shortest path properties of the largest connected component of a given graph.

    iG = convert_to_Graph_of_igraph(G)
    igraph_path_length_hist = iG.path_length_hist(directed=False)

    G.spld = defaultdict(float)
    num_all = igraph_path_length_hist.n
    bins = tuple(igraph_path_length_hist.bins())

    for (i, j, k) in bins:
        if j != i+1:
            print("Error.")
            exit(0)
        G.spld[i] = float(k)/num_all

    G.diameter = max(list(dict(G.spld).keys()))
    G.apl = sum([l*G.spld[l] for l in G.spld])

    return

def calc_betweenness(G: Graph):
    # Note: calculate shortest path of the largest connected component of a given graph.

    iG = convert_to_Graph_of_igraph(G)
    degrees = iG.degree(list(range(0, len(G.nodes))))
    bc = iG.betweenness(directed=False)
    n = int(iG.vcount())

    G.ddbc = defaultdict(float)
    V_d = defaultdict(int)
    for i in range(0, len(degrees)):
        d = degrees[i]
        G.ddbc[d] += float(bc[i])/((n-1)*(n-2))
        V_d[d] += 1

    for d in G.ddbc:
        G.ddbc[d] = float(G.ddbc[d])/V_d[d]

    return

def calc_largest_eigenvalue(G:Graph):
    iG = convert_to_Graph_of_igraph(G)
    L = iG.laplacian(normalized=True)
    eigenvalues = np.linalg.eigvals(L)
    G.lambda_1 = float(max(eigenvalues))

    return

def calc_properties(G: Graph):
    G.N = len(G.nodes)

    G.M = 0
    for i in G.nodes:
        G.M += len(G.nlist[i])
    G.M = int(G.M)/2

    G.aved = float(2*G.M)/G.N

    calc_dd(G)

    calc_knn(G)

    calc_clustering(G)

    calc_common_neighbor_distribution(G)

    calc_shortest_path_properties(G)

    calc_betweenness(G)

    calc_largest_eigenvalue(G)

    return

def calc_relative_error_for_scalar_property(G_property, genG_property):

    return float(math.fabs(G_property - genG_property))/G_property

def calc_normalized_L1_distance_for_distribution(G_property, genG_property):

    keys = set(G_property.keys()) | set(genG_property.keys())
    dist = 0
    norm = sum(list(G_property.values()))
    for key in keys:
        dist += math.fabs(G_property[key] - genG_property[key])

    return float(dist)/norm

def calc_error_of_each_proeprty(G:Graph, genG:Graph):
    calc_properties(G)
    calc_properties(genG)

    print("Normalized L1 distance of each property of a generated graph.")

    print("Number of nodes:", calc_relative_error_for_scalar_property(G.N, genG.N))

    print("Average degree:", calc_relative_error_for_scalar_property(G.aved, genG.aved))

    print("Degree distribution:", calc_normalized_L1_distance_for_distribution(G.dd, genG.dd))

    print("Neighbor connectivity:", calc_normalized_L1_distance_for_distribution(G.knn, genG.knn))

    print("Average local clustering coefficient:", calc_relative_error_for_scalar_property(G.acc, genG.acc))

    print("Degree-dependent clustering coefficient:", calc_normalized_L1_distance_for_distribution(G.ddcc, genG.ddcc))

    print("Common neighbor distribution:", calc_normalized_L1_distance_for_distribution(G.cnd,genG.cnd))

    print("Average shortest path length:", calc_relative_error_for_scalar_property(G.apl, genG.apl))

    print("Shortest path length distribution:", calc_normalized_L1_distance_for_distribution(G.spld,genG.spld))

    print("Diameter:", calc_relative_error_for_scalar_property(G.diameter,genG.diameter))

    print("Degree-dependent betweenness centrality:", calc_normalized_L1_distance_for_distribution(G.ddbc,genG.ddbc))

    print("Largest eigenvalue:", calc_relative_error_for_scalar_property(G.lambda_1,genG.lambda_1))
    
    return

