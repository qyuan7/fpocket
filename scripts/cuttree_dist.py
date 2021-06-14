#!/usr/bin/env python

import Pycluster
from Pycluster import *
from scipy.spatial import distance_matrix
from collections import Counter
import argparse

def pdb_to_data(fname):
    data = []
    with open(fname) as f:
        for line in f:
            if line.startswith("ATOM"):
                lst = [float(s) for s in line.split()[5:8]]
                data.append(lst)
    return data


def assign_id(fname, cids, default_id, fout, default_fout):
    data = pdb_to_data(fname)
    assert(len(data) == len(cids))
    with open(fout, "w") as f, open(default_fout, "w") as f2:
        for i in range(len(data)):
            f.write(f"ATOM  {i+1:5d}  C   PTH  {cids[i]:3d}     {data[i][0]:8.3f}{data[i][1]:8.3f}{data[i][2]:8.3f}{'0.00':>6}{'0.00':>6}\n")
            if cids[i] == default_id:
                f2.write(f"ATOM  {i+1:5d}  C   PTH  {cids[i]:3d}     {data[i][0]:8.3f}{data[i][1]:8.3f}{data[i][2]:8.3f}{'0.00':>6}{'0.00':>6}\n")
    return


def assign_cluster_to_children(cid, node_id, tree, cluster_assign):
    if node_id < 0:
        assign_cluster_to_children(cid, tree[node_id*(-1)-1].left, tree, cluster_assign)
        assign_cluster_to_children(cid, tree[node_id*(-1)-1].right, tree, cluster_assign)
    else:
        if(cluster_assign[node_id][0]!=1):
            cluster_assign[node_id][0] = 1
            cluster_assign[node_id][1] = cid
    return


def cuttree_distance(nelements, tree, distance):
    cluster_assign = [[-1,0] for _ in range(nelements)]
    cluster_id = 1
    for i in range(nelements-2, -1, -1):
        if tree[i].distance > distance:
            if tree[i].left >= 0:
                assign_cluster_to_children(cluster_id, tree[i].left, tree, cluster_assign)
                cluster_id += 1
            if tree[i].right >= 0:
                assign_cluster_to_children(cluster_id, tree[i].right, tree, cluster_assign)
                cluster_id += 1
        else:
            cluster_id += 1
            assign_cluster_to_children(cluster_id, tree[i].left, tree, cluster_assign)
            assign_cluster_to_children(cluster_id, tree[i].right, tree, cluster_assign)
    return cluster_assign


def reindex_clusters(cids):
    reindex = {}
    res = []
    newid = 0
    for cid in cids:
        if cid in reindex:
            res.append(reindex[cid])
        else:
            newid += 1
            reindex[cid] = newid
            res.append(newid)
    return res


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Clustering of pockets with certain isovalues.")
    parser.add_argument("--input", action="store", dest="infn", help="Input file containing gridpoints with isovalue cutoff")
    parser.add_argument("--output", action="store", dest="outfn", help="Output file with clustered gridpoints" )
    parser.add_argument("--bigpock", action="store", dest="bigpock", help="File for biggest single pocket detected by clustering")
    arg_dict = vars(parser.parse_args())
    infn_, outfn_, bigfn_ = arg_dict.values()
    data = pdb_to_data(infn_)

    #print(len(data), len(data[0]))
    dist_mtx = distance_matrix(data, data)
    tree = treecluster( data=None, distancematrix=dist_mtx, method='s',dist='e')
    clusters = cuttree_distance(len(data), tree, distance=1.74)
    cids = [f[1] for f in clusters]
    newids = reindex_clusters(cids)

    size = Counter(newids)
    default_id = size.most_common()[0][0] # id for the biggest pocket - analyse by default
    #print(size.most_common())
    assign_id(infn_, newids, default_id, outfn_, bigfn_)

    #assign_id("drug_freq_iso_0_5.pdb", newids, default_id, "drug_freq_iso_0_5_cluster.pdb", "drug_freq_iso_0_5_biggest.pdb")
