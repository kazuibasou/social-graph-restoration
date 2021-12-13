import graph
from collections import defaultdict
import random
import math

def rewirable(u, v, x, y, k_u, k_v, k_x, k_y):
    if u == v:
        return False
    if u == x:
        return False
    if u == y:
        return False
    if v == x:
        return False
    if v == y:
        return False
    if x == y:
        return False
    if k_u == k_x:
        return True
    if k_u == k_y:
        return True
    if k_v == k_x:
        return True
    if k_v == k_y:
        return True

    return False

def calculate_number_of_triangles_to_add(genG:graph.Graph, node_degree, u, v, y, num_tri_to_add):

    for w in genG.nlist[u]:
        if node_degree[w] <= 1 or u == w:
            continue

        if v != w and node_degree[v] > 1:
            t_minus = genG.nlist[v].count(w)
            num_tri_to_add[node_degree[u]] -= t_minus
            num_tri_to_add[node_degree[v]] -= t_minus
            num_tri_to_add[node_degree[w]] -= t_minus

        if y != w and node_degree[y] > 1:
            t_plus = genG.nlist[y].count(w)
            num_tri_to_add[node_degree[u]] += t_plus
            num_tri_to_add[node_degree[y]] += t_plus
            num_tri_to_add[node_degree[w]] += t_plus

    return num_tri_to_add

def rewiring_random_edge_pair_preserving_joint_degree_matrix(genG:graph.Graph, node_degree, rewirable_edges):
    i_e1 = random.randrange(0, len(rewirable_edges))
    i_e2 = random.randrange(0, len(rewirable_edges))
    e1 = rewirable_edges[i_e1]
    e2 = rewirable_edges[i_e2]

    u = e1[0]
    v = e1[1]
    x = e2[0]
    y = e2[1]

    k_u = node_degree[u]
    k_v = node_degree[v]
    k_x = node_degree[x]
    k_y = node_degree[y]

    while not rewirable(u,v,x,y,k_u,k_v,k_x,k_y):
        i_e1 = random.randrange(0, len(rewirable_edges))
        i_e2 = random.randrange(0, len(rewirable_edges))
        e1 = rewirable_edges[i_e1]
        e2 = rewirable_edges[i_e2]

        u = e1[0]
        v = e1[1]
        x = e2[0]
        y = e2[1]

        k_u = node_degree[u]
        k_v = node_degree[v]
        k_x = node_degree[x]
        k_y = node_degree[y]

    num_tri_to_add = defaultdict(int)
    rewiring_case = -1

    if k_u == k_x or k_v == k_y:
        num_tri_to_add = calculate_number_of_triangles_to_add(genG,node_degree,u,v,y,num_tri_to_add)
        num_tri_to_add = calculate_number_of_triangles_to_add(genG,node_degree,x,y,v,num_tri_to_add)

        if node_degree[v] > 1 and node_degree[y] > 1:
            t_minus = genG.nlist[v].count(y)
            num_tri_to_add[node_degree[u]] -= t_minus
            num_tri_to_add[node_degree[v]] -= 2 * t_minus
            num_tri_to_add[node_degree[x]] -= t_minus
            num_tri_to_add[node_degree[y]] -= 2 * t_minus

        if node_degree[u] > 1 and node_degree[x] > 1:
            t_minus = genG.nlist[x].count(u)
            if node_degree[v] > 1:
                num_tri_to_add[node_degree[u]] -= t_minus
                num_tri_to_add[node_degree[v]] -= t_minus
                num_tri_to_add[node_degree[x]] -= t_minus
            if node_degree[y] > 1:
                num_tri_to_add[node_degree[u]] -= t_minus
                num_tri_to_add[node_degree[x]] -= t_minus
                num_tri_to_add[node_degree[y]] -= t_minus

        rewiring_case = 0

    elif k_u == k_y or k_v == k_x:
        num_tri_to_add = calculate_number_of_triangles_to_add(genG,node_degree,u,v,x,num_tri_to_add)
        num_tri_to_add = calculate_number_of_triangles_to_add(genG,node_degree,y,x,v,num_tri_to_add)

        if node_degree[v] > 1 and node_degree[x] > 1:
            t_minus = genG.nlist[v].count(x)
            num_tri_to_add[node_degree[u]] -= t_minus
            num_tri_to_add[node_degree[v]] -= 2 * t_minus
            num_tri_to_add[node_degree[x]] -= 2 * t_minus
            num_tri_to_add[node_degree[y]] -= t_minus

        if node_degree[u] > 1 and node_degree[y] > 1:
            t_minus = genG.nlist[y].count(u)
            if node_degree[v] > 1:
                num_tri_to_add[node_degree[u]] -= t_minus
                num_tri_to_add[node_degree[v]] -= t_minus
                num_tri_to_add[node_degree[y]] -= t_minus
            if node_degree[x] > 1:
                num_tri_to_add[node_degree[u]] -= t_minus
                num_tri_to_add[node_degree[x]] -= t_minus
                num_tri_to_add[node_degree[y]] -= t_minus

        rewiring_case = 1

    else:
        print("Error: Selected edge pair destroys the present joint degree matrix.")
        exit()

    return [num_tri_to_add, i_e1, i_e2, rewiring_case]

def calc_L1_distance(tgt_ddcc, cur_ddcc):
    degrees = set(tgt_ddcc.keys()) | set(cur_ddcc.keys())
    dist = 0
    norm = sum(list(tgt_ddcc.values()))

    for d in degrees:
        dist += math.fabs(tgt_ddcc[d] - cur_ddcc[d])

    return [dist, norm]

def targeting_rewiring_for_clustering(genG:graph.Graph, tgt_ddcc, rewirable_edges, R_C=500):
    node_degree = defaultdict(int)
    n_d = defaultdict(int)
    for v in genG.nodes:
        d = len(genG.nlist[v])
        if d > 1:
            node_degree[v] = d
            n_d[d] += 1

    const_coeff = defaultdict(float)
    for d in n_d:
        if d > 1:
            const_coeff[d] = float(2)/(d*(d-1))
            const_coeff[d] = float(const_coeff[d])/n_d[d]

    graph.calc_clustering(genG)
    cur_ddcc = genG.ddcc.copy()

    [dist, norm] = calc_L1_distance(tgt_ddcc, cur_ddcc)

    R = R_C*len(rewirable_edges)

    for r in range(0, R):

        rewired_ddcc = cur_ddcc.copy()
        rewired_dist = dist

        [num_tri_to_add, i_e1, i_e2, rewiring_case] = rewiring_random_edge_pair_preserving_joint_degree_matrix(genG,node_degree,rewirable_edges)

        for d in num_tri_to_add:
            if d > 1:
                rewired_ddcc[d] += float(num_tri_to_add[d]*const_coeff[d])
                rewired_dist += math.fabs(tgt_ddcc[d] - rewired_ddcc[d]) - math.fabs(tgt_ddcc[d] - cur_ddcc[d])

        delta_dist = rewired_dist - dist

        if delta_dist >= 0:
            continue

        u = rewirable_edges[i_e1][0]
        v = rewirable_edges[i_e1][1]
        x = rewirable_edges[i_e2][0]
        y = rewirable_edges[i_e2][1]

        if rewiring_case == 0:
            graph.remove_edge(genG,u,v)
            graph.remove_edge(genG,x,y)
            graph.add_edge(genG,u,y)
            graph.add_edge(genG,v,x)

            tmp = rewirable_edges[i_e1][1]
            rewirable_edges[i_e1][1] = rewirable_edges[i_e2][1]
            rewirable_edges[i_e2][1] = tmp

        elif rewiring_case == 1:
            graph.remove_edge(genG, u, v)
            graph.remove_edge(genG, x, y)
            graph.add_edge(genG, u, x)
            graph.add_edge(genG, v, y)

            tmp = rewirable_edges[i_e1][1]
            rewirable_edges[i_e1][1] = rewirable_edges[i_e2][0]
            rewirable_edges[i_e2][0] = tmp

        cur_ddcc = rewired_ddcc.copy()
        dist = rewired_dist

    return genG

def check_update_num_tri(genG:graph.Graph,rewirable_edges,R = 1000):

    node_degree = defaultdict(int)
    n_d = defaultdict(int)
    for v in genG.nodes:
        d = len(genG.nlist[v])
        if d > 1:
            node_degree[v] = d
            n_d[d] += 1

    const_coeff = defaultdict(float)
    for d in n_d:
        if d > 1:
            const_coeff[d] = float(2) / (d * (d - 1))
            const_coeff[d] = float(const_coeff[d]) / n_d[d]

    graph.calc_num_tri(genG)
    cur_num_tri = genG.num_tri.copy()

    for r in range(0, R):

        [num_tri_to_add, i_e1, i_e2, rewiring_case] = rewiring_random_edge_pair_preserving_joint_degree_matrix(genG, node_degree, rewirable_edges)

        for d in num_tri_to_add:
            if d > 1:
                cur_num_tri[d] += num_tri_to_add[d]

        u = rewirable_edges[i_e1][0]
        v = rewirable_edges[i_e1][1]
        x = rewirable_edges[i_e2][0]
        y = rewirable_edges[i_e2][1]

        if rewiring_case == 0:
            graph.remove_edge(genG, u, v)
            graph.remove_edge(genG, x, y)
            graph.add_edge(genG, u, y)
            graph.add_edge(genG, v, x)

            tmp = rewirable_edges[i_e1][1]
            rewirable_edges[i_e1][1] = rewirable_edges[i_e2][1]
            rewirable_edges[i_e2][1] = tmp

        elif rewiring_case == 1:
            graph.remove_edge(genG, u, v)
            graph.remove_edge(genG, x, y)
            graph.add_edge(genG, u, x)
            graph.add_edge(genG, v, y)

            tmp = rewirable_edges[i_e1][1]
            rewirable_edges[i_e1][1] = rewirable_edges[i_e2][0]
            rewirable_edges[i_e2][0] = tmp

        graph.calc_num_tri(genG)
        correct_num_tri = genG.num_tri.copy()
        degrees = set(cur_num_tri.keys()) | set(correct_num_tri.keys())
        for d in degrees:
            if cur_num_tri[d] != correct_num_tri[d]:
                print("Error: Failed to correctly update numbers of triangles to which nodes with each degree belong.")
                exit()

    return 0
