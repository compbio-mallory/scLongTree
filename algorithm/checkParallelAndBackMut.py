# Author: Rituparna Khan

import argparse
import numpy as np
import math

''' Define the tree structure. '''
class Node:
    def __init__(self, id_):
        self.id = id_  # node or child id
        self.pID = -2 # parent node
        self.mutations = [] # mutations incoming on this node
        self.edges = "" # incoming edges
        self.children = [] # child nodes
        self.timepoint = 0
        self.cells = [] # save the cells belonging to this node
        self.type = "subclone"

def getParallelMut(Tree):
    """
    For each mutation x, collect all edges p->c where x is present in child but absent in parent.
    Return:
      - parallel_mut: list of mutation ids that appear as a 'new gain' on >= 2 edges
      - parallelMut_edges: dict {mut_id -> [ "p_c", ... ]} of ALL candidate new-gain edges
    """
    mutation_node = {}      # {mut -> [child_node_ids where mut newly appears]}
    parallelMut_edges = {}  # {mut -> ["p_c", ...]}

    # Walk edges  (note: assumes Tree indices match node ids)
    for j in range(len(Tree)):
        if j == 0:
            continue  # skip root
        nodeId = Tree[j].id
        p = int(Tree[nodeId].pID)
        if p == -2:
            continue  # removed node
        parent_muts = set(Tree[p].mutations)
        child_muts  = set(Tree[nodeId].mutations)
        new_mut = child_muts - parent_muts
        for mut in new_mut:
            # robust against bad/empty markers
            if mut == '' or mut is None:
                continue
            # record child node where this mut newly appears
            mutation_node.setdefault(int(mut), []).append(nodeId)
            parallelMut_edges.setdefault(int(mut), []).append(f"{p}_{nodeId}")

    # mutations that show up as a 'new gain' on multiple edges
    parallel_mut = [m for m, nodes in mutation_node.items() if len(nodes) > 1]

    # Ensure lists are unique & stable
    for m, edges in list(parallelMut_edges.items()):
        # dedup while preserving content
        parallelMut_edges[m] = list(dict.fromkeys(edges))

    print("Mutation nodes ", mutation_node)
    print(" Parallel mutations ", parallel_mut, " edges ", parallelMut_edges)
    return parallel_mut, parallelMut_edges


''' Get back mutation count. '''
def getBackMutCount(Tree):
    # back mutations are ones where a mutation present in parent is absent in child
    backMut = []
    backMut_edges = {}
    print(" BACK MUTATIONS =================================== ")
    for j in range(len(Tree)):
        if j == 0: # Skip root node
            continue
        p = int(Tree[j].pID)
        if p == -2: # This node doesn't exist after correction
            continue
        new_mut = set(Tree[j].mutations) - set(Tree[p].mutations)
        c_otherThan_newMut = set(Tree[j].mutations) - new_mut
        absent_mut = set(Tree[p].mutations) - c_otherThan_newMut  
        #print(j," New mutation ",new_mut," other than new mut ",c_otherThan_newMut," absent mut ",absent_mut)
        if len(absent_mut) != 0:
            for am in absent_mut:
                if am == '':
                    continue
                if am in backMut_edges:
                    tempList = backMut_edges[am]
                    tempList.append(Tree[j].edges)
                    backMut_edges[am] = tempList
                else:
                    backMut_edges[am] = [Tree[j].edges]
        #print(nodeId," New mutations ",new_mut)
    
    noOfEdges = 0
    for bm, edges in backMut_edges.items():
        noOfEdges = noOfEdges+len(edges)

    print(" Back mutations ",backMut_edges," no. of edges ",noOfEdges)
    print(" ================================================= ")
    return backMut_edges

# Input is tree, node to get its children, a nodes_list to keep all the nodes of a subtree.
# Output is the final nodes_list returning all nodes of a subtree.
''' Get the subtree given the child node as the root. '''
def nodeChildren(Tree, node, nodes_list):
    for j in range(len(Tree)):
        nodeId = Tree[j].id
        if nodeId == node:
            #print(" ID ",Tree[j].id," Children ",Tree[j].children)
            temp_list = []
            temp_list.extend(Tree[nodeId].children)
            nodes_list.extend(Tree[nodeId].children)
            #print(" Nodes list here ",nodes_list)
            if temp_list != []:
                for c in temp_list:
                    nodeChildren(Tree, c, nodes_list)
    return nodes_list

# Input is tree, root is the node from where we get the subtree, parallel mutation x.
# Output is all nodes of a subtree after checking for loss of mutation x.
''' Get the list of subtree. '''
def subtree(Tree, root, x):
    nodes_list = nodeChildren(Tree, root, [])
    for node in nodes_list:
        if x not in Tree[node].mutations: # If there is loss of mutation x then remove that node
            nodes_list.remove(node)
    return nodes_list

# Input is tree, parallel mutation x, edges of the parallel mutation x A_x, noisy D matrix, consensus genotype G, \alpha and \beta for timepoints
# Output is the edge with the maximum ratio R_e to place the parallel mutation.
''' Select the best edge to place the parallel mutation. '''
def selectTopEdgesParallel(Tree, x, A_x, D_matrix, timepoint_FP, timepoint_FN, max_parallel_extra):
    """
    Select up to (1 + max_parallel_extra) edges for mutation x, ranked by summed log-likelihood ratios
    across the child edge's subtree (parallel model).
    """
    edge_ratios = {}
    for e in A_x:
        child_node = int(e.split('_')[1])
        subtree_nodes = subtree(Tree, child_node, x)

        child_node_tp = Tree[child_node].timepoint
        node_cells = Tree[child_node].cells
        R_e = calculateRatio(node_cells, x, timepoint_FP[child_node_tp], timepoint_FN[child_node_tp], D_matrix, "parallel")

        for n in subtree_nodes:
            if Tree[n].type == "unobserved":
                continue
            n_cells = Tree[n].cells
            n_tp = Tree[n].timepoint
            R_n = calculateRatio(n_cells, x, timepoint_FP[n_tp], timepoint_FN[n_tp], D_matrix, "parallel")
            R_e = R_e + R_n

        edge_ratios[e] = R_e

    # pick top (1 + max_parallel_extra) by ratio
    keep_count = min(1 + max_parallel_extra, len(edge_ratios))
    sorted_edges = sorted(edge_ratios.items(), key=lambda kv: kv[1], reverse=True)
    selected = [e for (e, _) in sorted_edges[:keep_count]]

    print("Parallel mutation edge ratios:", edge_ratios)
    print(f"Selected {keep_count} edges for mut {x}: {selected}")
    return selected


# Input is tree, back mutation x, edges of the back mutation x B_x, noisy D matrix, consensus genotype G, \alpha and \beta for timepoints, k = no. of back mutation edges
# Output is the top k edges with the maximum summation of ratio R_e to allow back mutations.
''' Select the top k edges for the back mutation. '''
def selectBackMutEdges(Tree, x, B_x, D_matrix, timepoint_FP, timepoint_FN, k):
    if len(B_x) <= k:
        return B_x
    edge_ratios = {} # output dictionary with edge as key and ratio as values
    for e in B_x:
        child_node = int(e.split('_')[1])
        subtree_nodes = nodeChildren(Tree, child_node, [])
        #subtree_nodes = subtree(Tree, child_node, x)
        #print(" Subtree nodes ",subtree_nodes)

        child_node_tp = Tree[child_node].timepoint
        #child_node_tp = int(getNodeTimepoint(Tree, child_node))
        node_cells = Tree[child_node].cells
        #print(" Timepoint ",child_node_tp)
        R_e = calculateRatio(node_cells, x, timepoint_FP[child_node_tp], timepoint_FN[child_node_tp], D_matrix, "back")
        for n in subtree_nodes:
            if Tree[n].type == "unobserved": # skip calculating for the unobserved subclones as they have 0 cells.
                continue
            n_cells = Tree[n].cells
            n_tp = Tree[n].timepoint
            #n_tp = int(getNodeTimepoint(Tree, n))
            #print(" Timepoint ",n_tp)
            R_n = calculateRatio(n_cells, x, timepoint_FP[n_tp], timepoint_FN[n_tp], D_matrix, "back")
            R_e = R_e + R_n

        #print(" Ratio ",R_e)
        edge_ratios[e] = R_e
    
    # Select top k Edges. First sort the dictionary by their ratio in decreasing order and choose top k
    sorted_edge_ratios = dict(sorted(edge_ratios.items(), key=lambda item:item[1], reverse=True))
    kEdges = list(sorted_edge_ratios.keys())[:k]
    print("Back mutation sorted_edge_ratios ",sorted_edge_ratios)
    print(" k edges ",kEdges)
    return kEdges

# Input is cells in a node, parallel mutation x, \alpha and \beta of their respective timepoints, noisy D_matrix, mutation type
# Output is the ratio R_e for all cells in the node.
def calculateRatio(cells, x, alpha, beta, D_matrix, type_):
    R_cells = 0
    #print("FP rate ",alpha,"FN rate ",beta)
    for c in cells:
        #print(" Cell ",c," Mutation ",x)
        D_ix = D_matrix[int(c)][int(x)]
        if D_ix == 3:
            continue
        if D_ix == 0:
            P_Dix_Gx_1 = beta
            P_Dix_Gx_0 = 1 - alpha
        if D_ix == 1:
            P_Dix_Gx_1 = 1 - beta
            P_Dix_Gx_0 = alpha

        #print(P_Dix_Gx_1 / P_Dix_Gx_0)
        if type_ == "parallel":
            if (P_Dix_Gx_1 / P_Dix_Gx_0) == 0:
                R_cells = R_cells + np.log(float(1e-6))
            else:
                R_cells = R_cells + np.log(P_Dix_Gx_1 / P_Dix_Gx_0)
        else:
            if (P_Dix_Gx_0 / P_Dix_Gx_1) == 0:
                R_cells = R_cells + np.log(float(1e-6))
            else:
                R_cells = R_cells + np.log(P_Dix_Gx_0 / P_Dix_Gx_1)
    #print("R_cells ",R_cells)
    if R_cells == 0:
        return float(1e-6)
    else:
        return R_cells

# Save the selected edges for each parallel mutation
# Input is parallel mutations and their edges
# Output is a dictionary with parallel mutation as key and edge as value
def finalParallelMutEdges(Tree, parallel_mut, parallelMut_edges, D_matrix, tp_alpha, tp_beta, max_parallel_extra):
    """
    For each parallel mutation, keep up to (1 + max_parallel_extra) edges.
    Returns: dict {mutation -> [edge, edge, ...]}
    """
    selected_edges = {}
    for pm in parallel_mut:
        edges = selectTopEdgesParallel(Tree, pm, parallelMut_edges[pm], D_matrix, tp_alpha, tp_beta, max_parallel_extra)
        selected_edges[pm] = edges
    return selected_edges

# Save the selected edges for each back mutation
# Input is back mutations and their edges
# Output is a dictionary with back mutation as key and list of edges as value
def finalBackMutEdges(Tree, backMut_edges, D_matrix, tp_alpha, tp_beta, k):
    back_mut = list(backMut_edges.keys())
    selected_edges = {}
    for bm in back_mut:
        edges = selectBackMutEdges(Tree, bm, backMut_edges[bm], D_matrix, tp_alpha, tp_beta, k)
        selected_edges[bm] = edges
    return selected_edges

''' Place the selected parallel mutation on correct edges and update the Tree. '''
# Input is Tree, parallel edges for each mutation, selected parallel edges, back mutation edges, allowed back mutation edges
def correctParallelMut(Tree, parallel_edges, mut_edges):
    """
    Correct only mutations that were actually identified as parallel.

    parallel_edges:
        dict {mut -> [candidate occurrence edges]}
        This dictionary may contain ordinary single-occurrence mutations too.

    mut_edges:
        dict {parallel_mut -> [selected occurrence edges to keep]}
        Only mutations present in this dictionary should be corrected.
    """

    def find_all_occurrence_edges(mut):
        out = []

        for j in range(1, len(Tree)):
            p = Tree[j].pID

            if p == -2:
                continue

            parent_muts = set(Tree[p].mutations)
            child_muts = set(Tree[j].mutations)

            # Mutation is acquired on p -> j
            if (mut in child_muts) and (mut not in parent_muts):
                out.append(f"{p}_{j}")

        return out

    # IMPORTANT:
    # Iterate only over mutations that were actually selected/identified
    # as parallel mutations.  parallel_edges can contain ordinary
    # single-occurrence mutations and those must not be removed.
    for mut, selected in mut_edges.items():

        full_edges = set(find_all_occurrence_edges(mut))
        candidate_edges = set(parallel_edges.get(mut, []))

        # Include all observed occurrence edges for safety.
        all_edges = full_edges | candidate_edges

        keep_edges = set(selected)

        for e in all_edges:

            if e in keep_edges:
                continue

            child_node = int(e.split('_')[1])

            # Remove the unwanted parallel occurrence from this child.
            if mut in Tree[child_node].mutations:
                Tree[child_node].mutations.remove(mut)

            # Preserve the existing special-case behavior.
            pNode = (
                int(list(keep_edges)[0].split('_')[0])
                if keep_edges
                else None
            )

            if (
                pNode is not None
                and Tree[pNode].pID == child_node
            ):
                continue

            # Remove the unwanted occurrence throughout its subtree.
            subtree_nodes = nodeChildren(Tree, child_node, [])

            for n in subtree_nodes:
                if mut in Tree[n].mutations:
                    Tree[n].mutations.remove(mut)

    return Tree


''' Place the mutation on the edges where it shouldn't be lost. '''
def correctBackMut(Tree, backMut_edges, backEdges):
    # Correcting back mutations will mean checking the whole subtree and including the mutations there.
    for mut, edges in backMut_edges.items():
        edgeList = backEdges[mut]
        #print("Edges ",edges," edge list ",edgeList)
        
        for e in edges:
            if e in edgeList: # allow mutation loss
                continue
            parent_node = int(e.split('_')[0])
            child_node = int(e.split('_')[1])
            #print(" Child node ",child_node," mut ",mut)
            if mut not in Tree[child_node].mutations and mut in Tree[parent_node].mutations:
                Tree[child_node].mutations.append(mut)
            # add this mutation in the subtree 
            subtree_nodes = nodeChildren(Tree, child_node, [])
            for n in subtree_nodes:
                pID = Tree[n].pID
                if mut not in Tree[n].mutations and mut in Tree[pID].mutations:
                    Tree[n].mutations.append(mut)
    return Tree

# Input child node id and parent node id
# Output is an incoming edge of type string
def getEdges (Tree, nodeid, pID):
     Tree[nodeid].edges = str(pID)+"_"+str(nodeid)

# Input is tree
# Update the tree with its children
def getChildren(Tree):
    for j in range(1,len(Tree)):
        p = Tree[j].pID
        #print(j," parent ",p)
        Tree[p].children.append(j)
        #print(Tree[p].children)
        getEdges(Tree, j, p)

# Input is Tree, a dictionary cluster_cells where cluster is the key and cells [] value, and dictonary cluster_node
# Function updates the cells for each node
''' Get cells from the clusters. '''
def getCells (Tree, cluster_cells, cluster_node):
    for cluster, c in cluster_cells.items():
        nodeId = cluster_node[cluster]
        Tree[nodeId].cells = c

# Input is tree and unobserved node 
# Updates list of cells in its immediate descendants as ones for unobserved node
def getCellsForUnobservedNode (Tree, nodeId):
    children = Tree[nodeId].children
    for c in children:
        Tree[nodeId].cells.extend(Tree[c].cells)

''' Print the contents of the Tree after updating. '''
def printTree(Tree):
    print("=============== TREE =====================")
    for j in range(len(Tree)):
        print(" ID ",Tree[j].id," PID ",Tree[j].pID," Timepoint ",Tree[j].timepoint," Mutations ",Tree[j].mutations," Type ",Tree[j].type," Edges ",Tree[j].edges," Children ",Tree[j].children," Cells ",len(Tree[j].cells))
    print(" ====================================== ")

''' Read the input tree from the file and include it in the Tree structure defined here. '''
def readInputTree(treeFile):
    tf = open(treeFile,'r')
    tf_line = tf.readline().rstrip("\n")
    Tree = []
    unobserved_nodes = []
    while(tf_line != ""):
        tf_arr = tf_line.split('\t')
        nodeid = int(tf_arr[0])
        Tree.append(Node(nodeid))
        Tree[nodeid].id_ = int(nodeid)
        Tree[nodeid].timepoint = float(tf_arr[1])
        Tree[nodeid].pID = int(tf_arr[2])
        Tree[nodeid].mutations = tf_arr[4].split(',')
        Tree[nodeid].type = tf_arr[3]
        if "unobserved" in Tree[nodeid].type:
            unobserved_nodes.append(nodeid)
        #Tree[nodeid].edges = [tf_arr[6]]
        tf_line = tf.readline().rstrip("\n")

    return Tree, unobserved_nodes

''' Update the tree. '''
def updateTree(Tree, cluster_cells, cluster_node, unobserved_nodes):
    getChildren(Tree)
    getCells (Tree, cluster_cells, cluster_node)
    for node in reversed(unobserved_nodes):
        getCellsForUnobservedNode(Tree, node)


''' Save the D_matrix. '''
def readDMatrix(Dfile):
    dfile = open(Dfile,'r')
    df_line = dfile.readline().rstrip('\n')
    D_matrix = []
    while(df_line != ""):
        mut = df_line.split('\t')
        mut = [int(i) for i in mut]
        D_matrix.append(mut)
        df_line = dfile.readline().rstrip('\n')
    print(" No. of cells ",len(D_matrix)," mut ",len(D_matrix[0]))
    return D_matrix
