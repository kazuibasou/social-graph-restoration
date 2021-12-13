import numpy as np
from collections import defaultdict
import math
from functools import cmp_to_key
import random
import graph
import estimation
import rewiring

#Note: We use round to even.

def check_deg_vec(deg_vec):

    sum_deg_vec = 0
    for k in deg_vec:
        if deg_vec[k] < 0:
            print("Error: The number of nodes with degree " + str(k) + " is less than zero.")
            return False
        sum_deg_vec += k*deg_vec[k]

    if sum_deg_vec % 2 != 0:
        print("Error: The sum of degrees is not an even number.")
        return False

    return True

def check_jnt_deg_mat(jnt_deg_mat, deg_vec):

    for k in jnt_deg_mat:
        for l in jnt_deg_mat[k]:
            if k == l and jnt_deg_mat[k][l] % 2 != 0:
                print("Error: The number of edges between nodes with degree " + str(k) + " and nodes with degree " + str(l) + " is not an even number.")
                return False
            if jnt_deg_mat[k][l] != jnt_deg_mat[l][k]:
                print("Error: The number of edges between nodes with degree " + str(k) + " and nodes with degree " + str(l) + " is not symmetrical.")
                return False
            if jnt_deg_mat[k][l] < 0:
                print("Error: The number of edges between nodes with degree " + str(k) + " and nodes with degree " + str(l) + " is less than zero.")
                return  False
        sum_k = sum(list(jnt_deg_mat[k].values()))
        if sum_k != k*deg_vec[k]:
            print("Error: The sum of numbers of edges between nodes with degree " + str(k) + " and nodes with degree l for all l is not equal to " + str(k) + " times the number of nodes with degree " + str(k) + ".")
            return False

    return True

def select_random_key_with_smallest_value(dic):
    # Select a key with the smallest value in a given dictionary
    # If If there are two or more keys with the same value, we uniformly and randomly select between the keys.

    if len(dic) == 0:
        print("Error: The size of a given object is zero.")
        exit()

    min_value = float("inf")
    keys = set()
    for key in dic:
        if dic[key] < min_value:
            min_value = dic[key]
            keys = set()
            keys.add(key)
        elif dic[key] == min_value:
            keys.add(key)

    return np.random.choice(list(keys))

def select_min_key_with_smallest_value(dic):
    # Select a key with the smallest value in a given dictionary
    # If If there are two or more keys with the same value, we select the smallest key.

    if len(dic) == 0:
        print("Error: The size of a given object is zero.")
        exit()

    min_value = float("inf")
    keys = set()
    for key in dic:
        if dic[key] < min_value:
            min_value = dic[key]
            keys = set()
            keys.add(key)
        elif dic[key] == min_value:
            keys.add(key)

    return min(keys)

def initialize_target_degree_vector(est_n, est_dd):
    # Initialize the target degree vector

    tgt_deg_vec = defaultdict(int)
    for k in est_dd:
        tgt_deg_vec[k] = max(round(est_dd[k]*est_n), 1)

    return tgt_deg_vec

def adjust_target_degree_vector(est_n, est_dd, tgt_deg_vec):
    # Adjust the target degree vector

    sum_deg = sum([k*tgt_deg_vec[k] for k in tgt_deg_vec])

    if sum_deg % 2 == 0:
        return tgt_deg_vec

    degree_candidates = {}
    for k in tgt_deg_vec:
        if k % 2 == 0:
            continue
        x = est_dd[k]*est_n
        y = tgt_deg_vec[k]

        if x != 0:
            delta_e = float(math.fabs(x-y-1))/x - float(math.fabs(x-y))/x
        else:
            delta_e = float("inf")

        degree_candidates[k] = delta_e

    if len(degree_candidates) > 0:
        d = select_min_key_with_smallest_value(degree_candidates)
        tgt_deg_vec[d] += 1
    else:
        tgt_deg_vec[1] += 1

    return tgt_deg_vec

def cmp(a:list, b:list):
    if a[1] < b[1]:
        return -1
    elif a[1] > b[1]:
        return 1
    else:
        if a[0] < b[0]:
            return 1
        else:
            return -1

def modify_target_degree_vector(subG: graph.Graph, est_n, est_dd, tgt_deg_vec):
    # Assign the target degree of each node in the subgraph.
    # In parallel with this assignment process, modify the target degree vector.

    subG_deg_vec = defaultdict(int)
    tgt_node_deg = {}
    for v in subG.qry_nodes:
        subG_d = len(subG.nlist[v])
        tgt_node_deg[v] = subG_d
        subG_deg_vec[subG_d] += 1

    for d in subG_deg_vec:
        if tgt_deg_vec[d] < subG_deg_vec[d]:
            tgt_deg_vec[d] = subG_deg_vec[d]

    visible_node_pairs = []
    for v in subG.vis_nodes:
        visible_node_pairs.append([v, len(subG.nlist[v])])
    visible_node_pairs.sort(key=cmp_to_key(cmp))

    for visible_node_pair in visible_node_pairs:
        v = visible_node_pair[0]
        subG_d = visible_node_pair[1]

        degree_candidates = []
        for k in tgt_deg_vec:
            if k >= subG_d and tgt_deg_vec[k] > subG_deg_vec[k]:
                for i in range(0, tgt_deg_vec[k] - subG_deg_vec[k]):
                    degree_candidates.append(k)

        if len(degree_candidates) > 0:
            tgt_node_deg[v] = np.random.choice(list(degree_candidates))
        else:
            degree_to_add_candidates = {}
            for k in est_dd:
                if k < subG_d:
                    continue

                x = est_n*est_dd[k]
                y = float(tgt_deg_vec[k])
                if x != 0:
                    delta_e = float(math.fabs(x-y-1))/x - float(math.fabs(x-y))/x
                else:
                    delta_e = float("inf")

                degree_to_add_candidates[k] = delta_e

            if len(degree_to_add_candidates) > 0:
                tgt_node_deg[v] = select_min_key_with_smallest_value(degree_to_add_candidates)
            else:
                tgt_node_deg[v] = subG_d

            tgt_deg_vec[tgt_node_deg[v]] += 1

        subG_deg_vec[tgt_node_deg[v]] += 1

    tgt_deg_vec = adjust_target_degree_vector(est_n, est_dd, tgt_deg_vec)

    return [subG_deg_vec, tgt_deg_vec, tgt_node_deg]

def initialize_target_joint_degree_matrix(est_n, est_aved, est_jdd):
    # Initialize the target joint degree mat

    tgt_jnt_deg_mat = defaultdict(lambda: defaultdict(int))
    for k in est_jdd:
        for l in est_jdd[k]:
            if est_jdd[k][l] <= 0:
                continue
            x = round(est_n*est_aved*est_jdd[k][l])
            if k != l:
                tgt_jnt_deg_mat[k][l] = max(x, 1)
            else:
                if x % 2 == 0:
                    tgt_jnt_deg_mat[k][l] = max(x, 2)
                else:
                    y = est_n*est_aved*est_jdd[k][l]
                    if math.fabs(y-x+1) <= math.fabs(y-x-1):
                        tgt_jnt_deg_mat[k][l] = max(x-1,2)
                    else:
                        tgt_jnt_deg_mat[k][l] = max(x+1, 2)

    return tgt_jnt_deg_mat

def adjust_target_joint_degree_matrix(est_n, est_aved, est_jdd, tgt_deg_vec, min_jnt_deg_mat, tgt_jnt_deg_mat):
    # Adjust the target joint degree matrix

    degree_k1_set = set(tgt_deg_vec.keys())
    if 1 not in degree_k1_set:
        degree_k1_set.add(1)
    degree_k1_set = sorted(list(degree_k1_set), reverse=True)

    for k1 in degree_k1_set:
        target_sum = k1*tgt_deg_vec[k1]
        present_sum = sum(tgt_jnt_deg_mat[k1].values())
        diff = target_sum - present_sum

        if diff == 0:
            continue

        degree_k2_set = set([k2 for k2 in degree_k1_set if k2 <= k1])

        if k1 == 1 and abs(target_sum - present_sum) % 2 != 0:
            tgt_deg_vec[1] += 1
            target_sum += 1

        while target_sum != present_sum:
            if target_sum > present_sum:
                degree_k2_candidate = {}
                for k2 in degree_k2_set:
                    if present_sum == target_sum - 1 and k2 == k1:
                        continue

                    x = est_jdd[k1][k2]*est_n*est_aved
                    y = float(tgt_jnt_deg_mat[k1][k2])

                    if x == 0:
                        delta_e = float("inf")
                    else:
                        if k2 != k1:
                            delta_e = float(math.fabs(x-y-1))/x - float(math.fabs(x-y))/x
                        else:
                            delta_e = float(math.fabs(x-y-2))/x - float(math.fabs(x-y))/x
                    degree_k2_candidate[k2] = delta_e

                k2 = select_random_key_with_smallest_value(degree_k2_candidate)

                tgt_jnt_deg_mat[k1][k2] += 1
                tgt_jnt_deg_mat[k2][k1] += 1

                if k1 != k2:
                    present_sum += 1
                else:
                    present_sum += 2

            else:
                degree_k2_candidate = {}
                for k2 in degree_k2_set:
                    if tgt_jnt_deg_mat[k1][k2] <= min_jnt_deg_mat[k1][k2]:
                        continue
                    if present_sum == target_sum + 1 and k2 == k1:
                        continue

                    x = est_jdd[k1][k2] * est_n * est_aved
                    y = float(tgt_jnt_deg_mat[k1][k2])

                    if x == 0:
                        delta_e = float("inf")
                    else:
                        if k2 != k1:
                            delta_e = float(math.fabs(x - y + 1))/x - float(math.fabs(x - y))/x
                        else:
                            delta_e = float(math.fabs(x - y + 2))/x - float(math.fabs(x - y))/x

                    degree_k2_candidate[k2] = delta_e

                if len(degree_k2_candidate) > 0:
                    k2 = select_random_key_with_smallest_value(degree_k2_candidate)
                    tgt_jnt_deg_mat[k1][k2] -= 1
                    tgt_jnt_deg_mat[k2][k1] -= 1

                    if k1 != k2:
                        present_sum -= 1
                    else:
                        present_sum -= 2
                else:
                    if k1 > 1:
                        target_sum += k1
                        tgt_deg_vec[k1] += 1
                    else:
                        target_sum += 2
                        tgt_deg_vec[1] += 2

    return [tgt_jnt_deg_mat, tgt_deg_vec]

def modify_target_joint_degree_matrix(subG:graph.Graph, est_n, est_aved, est_jdd, tgt_node_deg, tgt_deg_vec, tgt_jnt_deg_mat):

    degree_set = set(tgt_deg_vec.keys())
    if 1 not in degree_set:
        degree_set.add(1)
    degree_set = set(sorted(list(degree_set)))

    subG_jnt_deg_mat = defaultdict(lambda: defaultdict(int))
    for v in subG.nodes:
        k1 = tgt_node_deg[v]
        for w in subG.nlist[v]:
            k2 = tgt_node_deg[w]
            subG_jnt_deg_mat[k1][k2] += 1

    for k1 in subG_jnt_deg_mat:
        for k2 in subG_jnt_deg_mat[k1]:
            while subG_jnt_deg_mat[k1][k2] > tgt_jnt_deg_mat[k1][k2]:
                tgt_jnt_deg_mat[k1][k2] += 1
                tgt_jnt_deg_mat[k2][k1] += 1

                degree_k3_candidates = {}
                for k3 in degree_set:
                    if k3 == k1 or tgt_jnt_deg_mat[k2][k3] <= subG_jnt_deg_mat.get(k2, {}).get(k3, 0):
                        continue

                    x = est_jdd[k2][k3]*est_n*est_aved
                    y = tgt_jnt_deg_mat[k2][k3]

                    if x == 0:
                        delta_e = float("inf")
                    else:
                        if k2 != k3:
                            delta_e = float(math.fabs(x-y+1))/x - float(math.fabs(x-y))/x
                        else:
                            delta_e = float(math.fabs(x-y+2))/x - float(math.fabs(x-y))/x

                    degree_k3_candidates[k3] = delta_e

                k3 = -1
                if len(degree_k3_candidates) > 0:
                    k3 = select_random_key_with_smallest_value(degree_k3_candidates)
                    tgt_jnt_deg_mat[k2][k3] -= 1
                    tgt_jnt_deg_mat[k3][k2] -= 1

                degree_k4_candidates = {}
                for k4 in degree_set:
                    if k4 == k2 or tgt_jnt_deg_mat[k1][k4] <= subG_jnt_deg_mat.get(k1, {}).get(k4, 0):
                        continue

                    x = est_jdd[k1][k4] * est_n * est_aved
                    y = tgt_jnt_deg_mat[k1][k4]

                    if x == 0:
                        delta_e = float("inf")
                    else:
                        if k1 != k4:
                            delta_e = float(math.fabs(x - y + 1)) / x - float(math.fabs(x - y)) / x
                        else:
                            delta_e = float(math.fabs(x - y + 2)) / x - float(math.fabs(x - y)) / x

                    degree_k3_candidates[k4] = delta_e

                if len(degree_k4_candidates) > 0:
                    k4 = select_random_key_with_smallest_value(degree_k4_candidates)
                    tgt_jnt_deg_mat[k4][k1] -= 1
                    tgt_jnt_deg_mat[k1][k4] -= 1

                    if k3 > 0:
                        tgt_jnt_deg_mat[k3][k4] += 1
                        tgt_jnt_deg_mat[k4][k3] += 1

    [tgt_jnt_deg_mat,tgt_deg_vec] = adjust_target_joint_degree_matrix(est_n,est_aved,est_jdd,tgt_deg_vec,subG_jnt_deg_mat,tgt_jnt_deg_mat)
    
    return [subG_jnt_deg_mat, tgt_jnt_deg_mat, tgt_deg_vec]

def graph_restoration_method(samplinglist, test=True):
    genG = graph.Graph()

    # (1) construct the subgraph
    node_index = {}
    edges_to_add = []

    # (1-1) index queried nodes
    i = 0
    for data in samplinglist:
        v = data.index
        if v not in node_index:
            node_index[v] = i
            i += 1
    genG.qry_nodes = set(range(0, i))

    # (1-2) index visible nodes
    marked = set()
    for data in samplinglist:
        v = data.index
        if v in marked:
            continue

        marked.add(v)
        for w in data.nlist:
            if w not in node_index:
                node_index[w] = i
                i += 1
            if w not in marked:
                edges_to_add.append([node_index[v], node_index[w]])

    genG.vis_nodes = set(range(len(genG.qry_nodes), i))

    # (1-3) Construct the subgraph
    genG.nodes = genG.qry_nodes | genG.vis_nodes
    for [v, w] in edges_to_add:
        graph.add_edge(genG, v, w)

    # (1-4) If there are no visible nodes, return the subgraph.
    if len(genG.vis_nodes) == 0:
        genG.N = len(genG.nodes)
        genG.M = 0
        genG.maxd = 0
        for v in genG.nodes:
            d = len(genG.nlist[v])
            genG.M += d
            if d > genG.maxd:
                genG.maxd = d
        genG.M = int(genG.M/2)

        return genG

    # (2) Construct the target degree vector
    est_n = estimation.size_estimator(samplinglist)
    est_dd = estimation.degree_distribution_estimator(samplinglist)

    tgt_deg_vec = initialize_target_degree_vector(est_n, est_dd)
    tgt_deg_vec = adjust_target_degree_vector(est_n, est_dd, tgt_deg_vec)
    [subG_deg_vec, tgt_deg_vec, tgt_node_deg] = modify_target_degree_vector(genG,est_n,est_dd,tgt_deg_vec)

    # (3) Construct the target joint degree matrix
    est_aved = estimation.average_degree_estimator(samplinglist)
    est_jdd = estimation.JDD_estimator_hybrid(samplinglist,est_n,est_aved)

    tgt_jnt_deg_mat = initialize_target_joint_degree_matrix(est_n,est_aved,est_jdd)
    min_jnt_deg_mat = defaultdict(lambda: defaultdict(int))
    [tgt_jnt_deg_mat, tgt_deg_vec] = adjust_target_joint_degree_matrix(est_n,est_aved,est_jdd,tgt_deg_vec,min_jnt_deg_mat,tgt_jnt_deg_mat)
    [subG_jnt_deg_mat, tgt_jnt_deg_mat, tgt_deg_vec] = modify_target_joint_degree_matrix(genG,est_n,est_aved,est_jdd,tgt_node_deg,tgt_deg_vec,tgt_jnt_deg_mat)

    if test:
        if not check_deg_vec(tgt_deg_vec):
            print("Error: Target degree vector does not satisfy realization conditions.")
            exit()
        if not check_jnt_deg_mat(tgt_jnt_deg_mat, tgt_deg_vec):
            print("Error: Target joint degree matrix does not satisfy realization conditions.")
            exit()
        for d in subG_deg_vec:
            if tgt_deg_vec[d] < subG_deg_vec[d]:
                print("Error: The target number of nodes with degree " + str(d) + "is smaller than the number of nodes with target degree " + str(d) + " in the subgraph.")
                exit()
        for k1 in subG_jnt_deg_mat:
            for k2 in subG_jnt_deg_mat[k1]:
                if tgt_jnt_deg_mat[k1][k2] < subG_jnt_deg_mat[k1][k2]:
                    print("Error: The target number of edges between nodes with degree " + str(k1) + " and nodes with degree " + str(k2) + " is smaller than the number of edges between nodes with target degree " + str(k1) + " and nodes with target degree " + str(k2) + " in the subgraph.")
                    exit()

    # (4) Construct a graph that preserves the target degree vector and the target joint degree matrix

    # (4-1) Determine the target number of nodes
    tgt_N = sum(list(tgt_deg_vec.values()))
    subG_N = len(genG.nodes)
    for v in range(subG_N, tgt_N):
        genG.nodes.add(v)

    # (4-2) Assign the target degree to each added node
    deg_seq = []
    for d in tgt_deg_vec:
        for i in range(0, tgt_deg_vec[d]-subG_deg_vec[d]):
            deg_seq.append(d)

    cur_deg_vec = defaultdict(int)
    for d in tgt_deg_vec:
        cur_deg_vec[d] = subG_deg_vec[d]

    random.shuffle(deg_seq)
    for v in range(subG_N, tgt_N):
        d = deg_seq.pop()
        tgt_node_deg[v] = d
        cur_deg_vec[d] += 1

    if test:
        for d in tgt_deg_vec:
            if tgt_deg_vec[d] != cur_deg_vec[d]:
                print("Error: A generated graph does not preserve the target degree vector.")
                exit()

        for d in cur_deg_vec:
            if cur_deg_vec[d] != tgt_deg_vec[d]:
                print("Error: A generated graph does not preserve the target degree vector.")
                exit()

    # (4-3) Make stub list
    stub_list = defaultdict(list)
    for v in genG.nodes:
        d = tgt_node_deg[v]
        subG_d = len(genG.nlist[v])
        for i in range(0, d-subG_d):
            stub_list[d].append(v)

    for d in stub_list:
        random.shuffle(stub_list[d])

    # (4-4) Connect each free stub of nodes with degree k1 and degree k2 uniformly at random
    cur_jnt_deg_mat = defaultdict(lambda: defaultdict(int))
    for k1 in tgt_jnt_deg_mat:
        for k2 in tgt_jnt_deg_mat[k1]:
            cur_jnt_deg_mat[k1][k2] = subG_jnt_deg_mat[k1][k2]

    for k1 in tgt_jnt_deg_mat:
        for k2 in tgt_jnt_deg_mat[k1]:
            while cur_jnt_deg_mat[k1][k2] != tgt_jnt_deg_mat[k1][k2]:
                u = stub_list[k1].pop()
                v = stub_list[k2].pop()
                graph.add_edge(genG, u, v)
                cur_jnt_deg_mat[k1][k2] += 1
                cur_jnt_deg_mat[k2][k1] += 1

    if test:
        for k1 in tgt_jnt_deg_mat:
            for k2 in tgt_jnt_deg_mat[k1]:
                if tgt_jnt_deg_mat[k1][k2] != cur_jnt_deg_mat[k1][k2]:
                    print("Error: A generated graph does not preserve the target joint degree matrix.")
                    exit()
        for k1 in cur_jnt_deg_mat:
            for k2 in cur_jnt_deg_mat[k1]:
                if tgt_jnt_deg_mat[k1][k2] != cur_jnt_deg_mat[k1][k2]:
                    print("Error: A generated graph does not preserve the target joint degree matrix.")
                    exit()

    genG.M = 0
    genG.maxd = 0
    for v in genG.nodes:
        d = len(genG.nlist[v])
        genG.M += d
        if d > genG.maxd:
            genG.maxd = d
    genG.M = int(genG.M/2)

    # (5) Targeting-rewiring process
    est_ddcc = estimation.degree_dependent_clustering_coefficient_estimator(samplinglist)

    rewirable_edges = []
    for v in range(len(genG.qry_nodes), len(genG.nodes)):
        for w in genG.nlist[v]:
            if w >= v and w >= len(genG.qry_nodes):
                rewirable_edges.append([v, w])

    if test:
        rewiring.check_update_num_tri(genG,rewirable_edges)

    genG = rewiring.targeting_rewiring_for_clustering(genG,est_ddcc,rewirable_edges)

    return genG





