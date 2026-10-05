# Author: Rituparna Khan

import argparse
import numpy as np
import pydot
import random
import os
import copy
from compareClusters import getNewClusters, updateTpClusterGen, checkClusterGenotype, reassignCells
from calculateClusterProb import calc_cluster_prob, timepoint_missingRate
from treeAlgorithm import Node, buildTree, printTree, connectSameNodes, find_unobservedSubclone, connect_remainingNodes, getEdges
from checkParallelAndBackMut import getParallelMut, getBackMutCount, finalParallelMutEdges, finalBackMutEdges, correctParallelMut, correctBackMut
from visualizeTree import drawTree
#from visualizeTree_large import drawTree

''' In a given timepoint calculate the likelihood of the cells for the clusters. '''
# Input is noisy D matrix, alpha beta values for each timepoint, cells to be reassigned, clusters and their genotype in each timepoint
def cell_cluster_likelihood(D_matrix, alpha, beta, cells, cluster_genotype, ignore_cluster, cluster_cells):
    # If there are no more subclones left to be merged in a timepoint then don't check any further
    if len(cluster_genotype) == 1 and ignore_cluster in cluster_genotype:
        print("Inside Done condition ",ignore_cluster)
        return "Done", "Done"

    del(cluster_cells[ignore_cluster]) # Remove the spurious cluster
    del(cluster_genotype[ignore_cluster])
   
    for i in cells:
        p_Di_Ck = 0 # i = cell, j = mutation, k = subclone
        total_mutations = len(D_matrix[0])
        cluster_prob = {}
        for cluster, gen in cluster_genotype.items():
            if cluster == ignore_cluster:
                continue
            ignore_mut = 0
            for j in range(len(gen)):
                if gen[j] == 3 or D_matrix[i][j] == 3:
                    ignore_mut = ignore_mut+1
                    continue
                if D_matrix[i][j] == 0 and gen[j] == 1:
                    L_j = beta
                if D_matrix[i][j] == 1 and gen[j] == 1:
                    L_j = 1 - beta
                if D_matrix[i][j] == 1 and gen[j] == 0:
                    L_j = alpha
                if D_matrix[i][j] == 0 and gen[j] == 0:
                    L_j = 1 - alpha
                p_Di_Ck = p_Di_Ck + np.log(L_j)
            if (total_mutations - ignore_mut) == 0:
                avg_prob = 0
            else:
                avg_prob = p_Di_Ck / (total_mutations - ignore_mut)
            cluster_prob[cluster] = avg_prob
        # select the cluster with max likelihood
        #print(" TP cluster ",cluster_prob," cell ",i)
        max_cluster = max(cluster_prob,key=cluster_prob.get)
        #print("Reassigned cell ",i," max cluster ",max_cluster)
        # Reassigning the cells to the selected cluster
        if i not in cluster_cells[max_cluster]:
            cluster_cells[max_cluster].append(i)

    return cluster_cells, cluster_genotype

''' Get the cluster genotype from BnpC's cell genotype without calculating max likelihood. '''
def getBnpcClonalGenotype(bnpc_cells_genotype, tp_cluster_cells):
    print("Inside BnpC Clonal genotype ============== ")
    bnpc_cluster_gen = {}
    #print(bnpc_cells_genotype)
    #print(tp_cluster_cells)
    tpClusterCells = copy.deepcopy(tp_cluster_cells)
    for tp, clusters in tpClusterCells.items():
        bnpc_cluster_gen[tp] = clusters
        for cluster, cells in clusters.items():
            cell_genotype = bnpc_cells_genotype[cells[0]]
            bnpc_cluster_gen[tp][cluster] = cell_genotype

    for tp, clusters in bnpc_cluster_gen.items():
        print("Timepoint ",tp)
        for cl, gen in clusters.items():
            print(cl," ",gen)
    #print(bnpc_cluster_gen)
    return bnpc_cluster_gen

''' After removing spurious subclones, reassigning cells we need to update cluster genotype. '''
def correctClustering(bnpc_cells_genotype, D_matrix, tp_cluster_cells, tp_alpha, tp_beta, tp_cluster_gen, initialFlag):
    if initialFlag:
        #updateTpClusterGen(D_matrix, tp_cluster_cells, tp_alpha, tp_beta)
        tp_cluster_gen = getBnpcClonalGenotype(bnpc_cells_genotype, tp_cluster_cells)
    print(" TP cluster gen ",tp_cluster_gen)
    tp_mergedClusters, tp_updatedCG = checkClusterGenotype(tp_cluster_gen)
    print(" TP merged clusters ",tp_mergedClusters)
    print(" TP updated CG ",tp_updatedCG)
    print(" TP cluster cells ",tp_cluster_cells)
    tp_reassignedCells = reassignCells(tp_mergedClusters, tp_cluster_cells)
    print(" TP reassigned cells ",tp_reassignedCells)
    return tp_reassignedCells, tp_updatedCG

''' Infer the tree using timepoint clusters, cells and genotype. '''
# Input is cells in each cluster and cluster genotype
def getTree(tp_cluster_cells, tp_cluster_genotype):
    print("INSIDE GET TREE ===========================")
    print(" TP cluster cells ",tp_cluster_cells)
    Tree, clone_node, tp_nodes = buildTree(tp_cluster_cells, tp_cluster_genotype)
    #print("Clone nodes ",clone_node)
    Tree, connectedNodes = connectSameNodes(Tree)
    Tree, tp_usc_nodes = find_unobservedSubclone(Tree, tp_nodes, connectedNodes)
    print("Unobserved subclones in each timepoint ",tp_usc_nodes)
    Tree = connect_remainingNodes(Tree, tp_nodes, tp_usc_nodes)
    Tree = getEdges(Tree)
    printTree(Tree)
    return Tree, clone_node

''' Save the Tree in a .csv file. '''
def saveTree(Tree, out_f):
    f = open(out_f, "w")
    f.write("\t".join(['ID','PID','Timepoint','Type','Mutations','Edges','Children','Cells'])+ "\n")
    for i in range(len(Tree)):
        #print(",".join(map(str,G[i].parentID)))
        #print(" ID ",Tree[j].id," PID ",Tree[j].pID," Timepoint ",Tree[j].timepoint," Mutations ",Tree[j].mutations," Type ",Tree[j].type," Edges ",Tree[j].edges," Children ",Tree[j].children," Cells ",len(Tree[j].cells))
        f.write("\t".join([str(Tree[i].id), str(Tree[i].pID), str(Tree[i].timepoint), str(Tree[i].type), ",".join(map(str,Tree[i].mutations)), str(Tree[i].edges), ",".join(map(str,Tree[i].children)), ",".join(map(str,Tree[i].cells))]) + "\n")
    f.close()

''' Get details from the Tree and generate plots. '''
def plotTree(Tree, plMut_edges, backMut_edges, opFile, name, sample):
    tp_nodes = {}
    id_cells = {}
    id_newMut = {}
    id_plMut = {}
    id_backMut = {}
    id_edges = {}
    usc_nodes = []
    skip_nodes = []
    for i in range(len(Tree)):
        nodeId = Tree[i].id
        if Tree[nodeId].pID == -2: # Don't consider this node for plotting
            skip_nodes.append(nodeId)
            continue
        if Tree[nodeId].type == "unobserved":
            cList = Tree[nodeId].children
            for c in cList:
                #if Tree[c].pID == nodeId:
                usc_nodes.append(Tree[nodeId].id)
                #else:
                #    skip_nodes.append(nodeId)
                    
        #if nodeId in skip_nodes:
        #    continue
        #print("Skip nodes ",skip_nodes)
        # Add the nodes in each timepoint
        if Tree[nodeId].timepoint in tp_nodes:
            tp_nodes[Tree[nodeId].timepoint].append(i)
        else:
            tp_nodes[Tree[nodeId].timepoint] = [nodeId]
        id_cells[nodeId] = Tree[nodeId].cells
        id_edges[nodeId] = Tree[nodeId].edges
        # Get the new and back mutation for this node
        p_mut = Tree[Tree[nodeId].pID].mutations
        new_mut = set(Tree[nodeId].mutations) - set(p_mut)
        id_newMut[nodeId] = list(new_mut)

    for m, edges in backMut_edges.items():
        for e in edges:
            nodeId = int(e.split('_')[1])
            if nodeId in skip_nodes: # skip adding nodes already removed
                continue
            if nodeId in id_backMut:
                id_backMut[nodeId].append(m)
            else:
                id_backMut[nodeId] = [m]

    #if name == "Final Tree":
    #    print("Back mutation ",id_backMut)

    for p, edges in plMut_edges.items():
        for e in edges:
            if e == '':
                continue
            nodeId = int(e.split('_')[1])
            if nodeId in id_plMut:
                id_plMut[nodeId].append(p)
            else:
                id_plMut[nodeId] = [p]

    plottedTree = drawTree(tp_nodes, id_cells, id_newMut, id_plMut, id_backMut, id_edges, usc_nodes, name, sample)
    #print(plottedTree)
    plottedTree.write_pdf(opFile)
    #print(" Timepoint nodes ",tp_nodes)
    #print(" ID cells ",id_cells)
    #print(" ID new mut ",id_newMut)
    #print(" ID edges ",id_edges)
    #print(" ID back mut ",id_backMut)

''' Given the probability of all the clones multiply their probabilities. '''
def calculate_treeProb(tp_cluster_prob):
    tree_prob = 0
    # p exists as log(cluster_prob) so we sum log of probabilities to prevent underflow
    for c,p in tp_cluster_prob.items():
        tree_prob = tree_prob + p
    return tree_prob

''' Select only those subclones < upper bound for each timepoint and sort them increasingly with probabilities. '''
# Input is the cluster prob threshold for each timepoint and prob of clusters in each timepoint
def getSpuriousSubclones(tp_thresholds, cluster_prob):
    selected_clusters = {}
    for cluster, prob in cluster_prob.items():
        tp = int(cluster.split("_")[0])
        if prob < tp_thresholds[tp]:
            selected_clusters[cluster] = prob

    # Then sort the selected clusters with increasing probability values
    sort_selected_clusters = dict(sorted(selected_clusters.items(), key=lambda x:x[1]))
    #sort_selected_clusters = dict(sort_selected_clusters)
    #sort_selected_clusters = list(sort_selected_clusters.keys())
    print(" Subclones < upper bound and sorted by increasing probabilities.")
    #print(" Cluster prob ",cluster_prob)
    print(" Sort selected clusters ",sort_selected_clusters)
    return sort_selected_clusters

''' Get weighted thresholds for all timepoints. '''
# Input is prob of clusters in each timepoint and cells in each timepoint
def clusterUpperBound(tp_cluster_prob, tpCells, tp_cluster_cells):
    tp_thresholds = {}
    tp_prob = {}
    
    for tp_cluster, prob in tp_cluster_prob.items():
        tp = int(tp_cluster.split('_')[0])
        cluster = int(tp_cluster.split('_')[1])
        cells = tp_cluster_cells[tp][cluster]
        cluster_prob = prob * len(cells)

        if tp in tp_prob:
            tp_prob[tp] = tp_prob[tp] + cluster_prob
        else:
            tp_prob[tp] = cluster_prob

    for i in range(len(tpCells)):
        #print(" Timepoint ",i," cells ",len(tpCells[i]))
        tp_thresholds[i] = tp_prob[i] / len(tpCells[i])

    print(" Timepoint thresholds ",tp_thresholds)
    # Find the max weighted avg of all time points which will be the upper bound. Return the upper bound.
    #upper_bound = max(tp_thresholds.values())
    #print(" Upper bound of the threshold ",upper_bound)
    return tp_thresholds

''' After correction some subclones might need merging for having same set of mutations. '''
def mergeClones(Tree):
    # Get the time point nodes. For each node in same timepoint check and merge with cells reassigned.
    tp_nodes = {}
    # First save the observed nodes in each timepoint
    for i in range(len(Tree)):
        nodeId = Tree[i].id
        if Tree[nodeId].type == "unobserved" or Tree[nodeId].pID == -2:
            continue
        tp = Tree[nodeId].timepoint
        if tp in tp_nodes:
            nodes = tp_nodes[tp]
            nodes.append(nodeId)
            tp_nodes[tp] = nodes
        else:
            tp_nodes[tp] = [nodeId]

    print("Timepoint nodes ",tp_nodes)
    # Commpare mutations of each node in a timepoint
    for i in range(len(Tree)):
        nodeId = Tree[i].id
        if Tree[nodeId].type == "unobserved" or Tree[nodeId].pID == -2:
            continue
        tp = Tree[nodeId].timepoint
        timepoint_nodes = tp_nodes[tp]

        for n in timepoint_nodes:
            if n == nodeId or Tree[nodeId].pID == -2:
                break
            if set(Tree[nodeId].mutations) == set(Tree[n].mutations):
                print("Checking nodes ",nodeId," ",n)
                Tree[n].children.extend(Tree[nodeId].children)
                Tree[n].cells.extend(Tree[nodeId].cells)
                Tree[Tree[nodeId].pID].children.remove(nodeId)
                print("Removed node ",nodeId)
                for c in Tree[nodeId].children:
                    Tree[c].pID = n
                    Tree[c].edges = str(n)+"_"+str(c)
                Tree[nodeId].pID = -2
                Tree[nodeId].children = []
                Tree[nodeId].edges = ""

    #print("Tree after merging ===================")
    #printTree(Tree)
    return Tree

''' Re-check and correct the Tree to remove or merge unobserved subclones after fixing parallel and back mutations. '''
def recheckTree(Tree):
    #printTree(Tree)
    # Check if at any time point two subclones have same set of mutations and merge them.
    #Tree = mergeClones(Tree)
    # Then check after correction if any unobserved subclones are unnecessary
    #nodesToRemove = set()
    for i in range(len(Tree)):
        #print(" Tree i ",i)
        nodeId = Tree[i].id
        if Tree[nodeId].pID == -2: # Don't check for already removed clones
            continue
        nodeId = Tree[i].id
        if Tree[nodeId].type == "unobserved":
            #if Tree[i].pID == 0:
            # Check for dangling nodes
            for c in Tree[nodeId].children:
                if Tree[c].pID == -2 or Tree[c].pID != nodeId:
                    Tree[nodeId].children.remove(c)
            if len(Tree[nodeId].children) == 1: # If there is just 1 children then remove this unobserved subclone
                Tree[Tree[nodeId].children[0]].pID = Tree[nodeId].pID
                Tree[Tree[nodeId].children[0]].edges = str(Tree[nodeId].pID)+"_"+str(Tree[nodeId].children[0])
                print("NodeId ",nodeId," pID ",Tree[nodeId].pID)
                Tree[Tree[nodeId].pID].children.remove(nodeId)
                Tree[Tree[nodeId].pID].children.extend(Tree[nodeId].children)
                Tree[nodeId].pID = -2
                Tree[nodeId].children = []
                Tree[nodeId].edges = ""
            elif set(Tree[Tree[nodeId].pID].mutations) == set(Tree[nodeId].mutations):
                #print(Tree[Tree[i].pID].mutations," ",Tree[i].mutations)
                Tree[Tree[nodeId].pID].children.remove(nodeId)
                Tree[Tree[nodeId].pID].children.extend(Tree[nodeId].children)
                for c in Tree[nodeId].children:
                    Tree[c].pID = Tree[nodeId].pID
                    Tree[c].edges = str(Tree[nodeId].pID)+"_"+str(c)
                print(" Nodes to remove ",nodeId)
                Tree[nodeId].pID = -2
                Tree[nodeId].children = []
                Tree[nodeId].edges = ""
                #nodesToRemove.add(i)
    return Tree

''' Given a tree check the unobserved subclones in same timepoint. If there can be a parental relationship between them then do that. '''
def checkUnobservedSubclones(Tree):
    # First get the unobserved subclones in each timepoint
    tp_usc = {}
    for i in range(len(Tree)):
        nodeId = Tree[i].id
        # Don't check for already removed clones or observed subclones
        if Tree[nodeId].pID == -2 or Tree[nodeId].type != "unobserved":
            continue
        tp = Tree[nodeId].timepoint
        if tp in tp_usc:
            tp_usc[tp].append(nodeId)
        else:
            tp_usc[tp] = [nodeId]

    # Note the parent child relationships.
    parent_child = {}
    for t, nodes in tp_usc.items():
        for n1 in nodes:
            n1_mut = Tree[n1].mutations
            for n2 in nodes:
                if n1 == n2:
                    continue
                n2_mut = Tree[n2].mutations
                # If parent is subset then connect them
                if set(n2_mut).issubset(set(n1_mut)):
                    if n2 in parent_child:
                        parent_child[n2].append(n1)
                    else:
                        parent_child[n2] = [n1]
    print("Connect unobserved subclones ",parent_child)

    # Before connecting the unobserved subclones first check if they are already connected
    for p, cIDs in parent_child.items():
        for c in cIDs:
            # If already connected to same parent or to another unobserved subclone (both parent or child) then skip
            c_child = Tree[c].children
            skip = False
            for cc in c_child:
                if Tree[cc].type == "unobserved":
                    skip = True
            if skip or Tree[Tree[c].pID].type == "unobserved" or Tree[c].pID == p:
                print("Skipped node ",c)
                continue

            Tree[c].pID = p
            Tree[c].edges = str(p)+"_"+str(c)
            Tree[p].children.append(c)
    return Tree

''' Recheck final tree before plotting. '''
def recheckFinalTree(Tree):
    # Merge clones having similar set of mutations.
    # Then check for unobserved subclones.
    Tree = checkUnobservedSubclones(Tree)
    Tree = mergeClones(Tree)
    Tree = recheckTree(Tree)
    print("Tree after merging and rechecking ===================")
    printTree(Tree)
    return Tree

''' Update the genotype of the clusters from the corrected Tree. '''
# Input is Tree, cluster genotype in each timepoint, clusters to nodes in Tree.
def genotypeFromTree(Tree, tp_cluster_gen, clone_node, noOfMutations):
    print("Clone node ",clone_node)
    print(" tp_cluster_gen ",tp_cluster_gen)
    for c,n in clone_node.items():
        tp = int(c.split('_')[0])
        cluster = int(c.split('_')[1])
        mut = Tree[n].mutations
        if mut == []:
            genotype = [0] * noOfMutations
        else:
            genotype = [0] * noOfMutations
            for n in mut:
                genotype[n] = 1
        tp_cluster_gen[tp][cluster] = genotype
    print("Genotype ",tp_cluster_gen)
    return tp_cluster_gen

#''' If any parallel edge appear as back mutation then keep that as back mutation. '''
# Input is parallel and mutation egdes
#def filterParallelEdges(parallelEdges, backEdges):
#    for bm in backEdges:
#        if bm in parallelEdges:
#            del(parallelEdges[bm])
#    return parallelEdges

''' Prints the mutation from genotype '''
def printMutFromGenotype(genotype):
    for c,g in genotype.items():
        mut = []
        for i in range(len(g)):
            if g[i] == 1:
                mut.append(i)
        print(str(c)," ",mut)

''' Note the diff between the old and new FP FN rates. '''
def diffFPFN(old_tp_FP, old_tp_FN, new_tp_FP, new_tp_FN):
    tp_diffFP = {}
    tp_diffFN = {}
    for tp in old_tp_FP:
        print(" Timepoint ",tp," Old FP ",old_tp_FP[tp]," Old FN ",old_tp_FN[tp])
        print(" Timepoint ",tp," New FP ",new_tp_FP[tp]," New FN ",new_tp_FN[tp])
        diff_FP = new_tp_FP[tp] - old_tp_FP[tp]
        diff_FN = new_tp_FN[tp] - old_tp_FN[tp]
        print(" Timepoint ",tp," Diff FP ",diff_FP," Diff FN ",diff_FN)
        tp_diffFP[tp] = diff_FP
        tp_diffFN[tp] = diff_FN
    return tp_diffFP, tp_diffFN

''' Check diff between FP FN rates. If it is too high then we don't drop this spurious subclone. '''
# Input is FP FN in each timepoint, timepoint, FP FN rate thresholds for each timepoint
# If error rates are higher than given threshold then return True to indicate that diff is high and change is not allowed.
def check_FPFN(tp_FP, tp_FN, tp, FPrates, FNrates):
    result = False
    #if tp == "all": # check for diff for all timepoints 
    #    for tp in tp_diffFP:
    #        FP = tp_diffFP[tp]
    #        FN = tp_diffFN[tp]
    #        if abs(FP) >= userRate or abs(FN) >= userRate:
    #            result = True
    #else:
    FP = tp_FP[tp]
    FN = tp_FN[tp]
    print("Timepoint FP ",tp_FP," FN ",tp_FN)
    print("FPrates ",FPrates,"FN rates ",FNrates)
    if FP > FPrates[tp] or FN > FNrates[tp]:
        result = True
    return result

''' All the clustering results will are here: Re-assigning cells after dropping clusters, correcting genotype, calculating prob after merging. 
Returns the updated cluster genotypes and cells assigned to each cluster. '''
def clusteringResults(D_matrix, tp, tp_cluster_cells, tp_cluster_genotype, tp_alpha, tp_beta, tp_MR, cells, cluster, noOfiter, bnpc_tp_FP, bnpc_tp_FN):
    print("Error rates used for updating all changes",tp_alpha," ",tp_beta)
    #print("Cells to be reassigned ",cells)
    # Reassign cells and update the cluster genotype after removing spurious subclone
    cell_cluster, cluster_genotype = cell_cluster_likelihood(D_matrix, bnpc_tp_FP[tp], bnpc_tp_FN[tp], cells, tp_cluster_genotype[tp], cluster, tp_cluster_cells[tp]) 
    # If nothing to reassign then we break with Tree with highest prob
    if cell_cluster == "Done" and cluster_genotype == "Done":
        return "NA", "NA"
    tp_cluster_cells[tp] = cell_cluster
    # Update new cluster genotypes
    tp_cluster_genotype[tp] = cluster_genotype
    print(" Cluster Genotype before reassigning cells =============== ")
    printMutFromGenotype(tp_cluster_genotype[tp])
    # After cell reassigning get the updated corrected cluster genotype 
    tp_cluster_cells, tp_cluster_genotype = correctClustering(D_matrix, D_matrix, tp_cluster_cells, tp_alpha, tp_beta, tp_cluster_genotype, False)

    print("Cluster Genotype after reassigning cells =============== ")
    printMutFromGenotype(tp_cluster_genotype[tp])

    # Prob after merging clones and reassigning cells
    #afterMergeProb = calc_cluster_prob(tp_cluster_cells, tp_cluster_genotype, D_matrix, new_tp_FP, new_tp_FN, tp_MR)
    #print("Tp_cluster_prob ",tp_cluster_prob)
    #print(" After merging prob ",afterMergeProb)
    return tp_cluster_genotype, tp_cluster_cells

''' Correct the parallel and back mutations and return new Tree and corrected cluster genotype. '''
def correctParallelAndBackMut(D_matrix, tp_cluster_cells, tp_cluster_genotype,
                              tp_alpha, tp_beta, Tree, clone_node, k,
                              plotOp, noOfiter, sample, max_parallel_extra):

    # ============================================================
    # 1. BACK MUTATIONS FIRST
    # ============================================================
    #
    # A mutation classified as a back mutation gets precedence.
    # After it has been classified as back-mutated, it must not
    # subsequently be treated as a parallel mutation.
    #
    backMut_edges = getBackMutCount(Tree)

    backEdges = finalBackMutEdges(
        Tree,
        backMut_edges,
        D_matrix,
        tp_alpha,
        tp_beta,
        k
    )

    print("Selected back edges ", backEdges)

    # Keep selected loss edges and repair only unselected losses.
    Tree = correctBackMut(
        Tree,
        backMut_edges,
        backEdges
    )

    print("Tree after correcting back mutations ")
    printTree(Tree)

    # Only mutations with an actually retained/allowed back-mutation
    # edge get precedence over parallel-mutation correction.
    # Mutations whose losses were repaired (e.g. when k=0) must remain
    # eligible for parallel correction.
    back_mut_set = {
        mut
        for mut, edges in backEdges.items()
        if len(edges) > 0
    }

    # ============================================================
    # 2. PARALLEL MUTATIONS SECOND
    # ============================================================
    #
    # Detect parallel mutations only after back-mutation correction.
    #
    parallel_mut, parallelMut_edges = getParallelMut(Tree)

    # Do not allow a mutation already classified as a back mutation
    # to be treated as parallel.
    parallel_mut = [
        mut
        for mut in parallel_mut
        if mut not in back_mut_set
    ]

    parallelMut_edges = {
        mut: edges
        for mut, edges in parallelMut_edges.items()
        if mut not in back_mut_set
    }

    print("Parallel mutations after excluding back mutations ",
          parallel_mut)

    print("Parallel edges after excluding back mutations ",
          parallelMut_edges)

    parallelEdges = finalParallelMutEdges(
        Tree,
        parallel_mut,
        parallelMut_edges,
        D_matrix,
        tp_alpha,
        tp_beta,
        max_parallel_extra
    )

    print("Selected parallel edges ", parallelEdges)

    Tree = correctParallelMut(
        Tree,
        parallelMut_edges,
        parallelEdges
    )

    print("Tree after correcting parallel mutations ")
    printTree(Tree)

    # ============================================================
    # 3. FINAL TREE CLEANUP
    # ============================================================

    Tree = recheckTree(Tree)

    print(
        "========== TREE after correcting BACK and "
        "PARALLEL MUTATIONS =========== "
    )
    printTree(Tree)

    # Update clone genotypes using the final corrected tree.
    new_tp_cluster_genotype = genotypeFromTree(
        Tree,
        tp_cluster_genotype,
        clone_node,
        len(D_matrix[0])
    )

    return (
        Tree,
        new_tp_cluster_genotype,
        parallelEdges,
        backEdges
    )



''' Map the cluster no. to nodes to get their probability. '''
def getNodeProb(spurious_subclone, cluster_nodes):
    node_prob = {}
    for ss, prob in spurious_subclone.items():
        if ss in cluster_nodes:
            node_prob[cluster_nodes[ss]] = prob
    print("Node prob ",node_prob)

''' Get the initial tree along with its estimated FP, FN rates.'''
# Input is D_matrix, cluster probabilities, cells in each cluster, cluster genotype.
# Output is initial Tree, clone node, Tree prob, initial FP FN rates
def getInitialTree(D_matrix, tp_cluster_prob, tp_cluster_cells, tp_cluster_genotype):
    Tree, clone_node = getTree(tp_cluster_cells, tp_cluster_genotype)
    tree_prob = calculate_treeProb(tp_cluster_prob)
    initial_tp_FP, initial_tp_FN = timepointFPFN(D_matrix, tp_cluster_cells, tp_cluster_genotype)
    return Tree, clone_node, tree_prob, initial_tp_FP, initial_tp_FN


# Return all mutations currently represented by at least one clone.
def getPresentMutations(tp_cluster_genotype):
    present = set()

    for tp, clusters in tp_cluster_genotype.items():
        for cluster, genotype in clusters.items():
            for mut, value in enumerate(genotype):
                if value == 1:
                    present.add(mut)

    return present


# Protect clusters that are the sole remaining carrier of a mutation.
def getMutationPreservingClusters(tp_cluster_genotype):
    mutation_holders = {}

    for tp, clusters in tp_cluster_genotype.items():
        for cluster, genotype in clusters.items():
            for mut, value in enumerate(genotype):
                if value == 1:
                    mutation_holders.setdefault(mut, []).append((tp, cluster))

    protected = {}

    for mut, holders in mutation_holders.items():
        if len(holders) == 1:
            holder = holders[0]
            protected.setdefault(holder, []).append(mut)

    return protected

''' All the steps done for each cluster to decide if it should be part of the Tree. '''
def iterationSteps(D_matrix, tp_cluster, tp_cluster_genotype, tp_cluster_cells, tp_MR, noOfiter, tp_alpha, tp_beta, tp_FPrate_threshold, tp_FNrate_threshold, preserve_mutations=False):
    tp = int(tp_cluster.split('_')[0])
    cluster = int(tp_cluster.split('_')[1])
    print(" Timepoint ",tp," checking if we can ignore cluster ",cluster)
    print(" Cluster genotype before merging ")
    if cluster not in tp_cluster_cells[tp]:
        return None, tp_cluster_genotype, tp_cluster_cells, None, None

    # Optional mutation-preserving mode:
    # skip this pruning attempt if this cluster is currently
    # the sole remaining carrier of one or more mutations.
    if preserve_mutations:
        protected_clusters = getMutationPreservingClusters(tp_cluster_genotype)
        key = (tp, cluster)

        if key in protected_clusters:
            print(
                "MUTATION-PRESERVING MODE: cluster",
                cluster,
                "in timepoint",
                tp,
                "will not be pruned; sole carrier of mutations",
                protected_clusters[key]
            )
            return None, tp_cluster_genotype, tp_cluster_cells, None, None

    printMutFromGenotype(tp_cluster_genotype[tp])
    cells = tp_cluster_cells[tp][cluster]
    # old_tp_FP, old_tp_FN are FP FN rates from previous iterations which means it is the updated FP FN rates
    old_tp_FP, old_tp_FN = timepointFPFN(D_matrix, tp_cluster_cells, tp_cluster_genotype)

    # Since these variables are changed in other classes I keep a copy to revert them
    prev_tp_cluster_genotype = copy.deepcopy(tp_cluster_genotype)
    prev_tp_cluster_cells = copy.deepcopy(tp_cluster_cells)
    #print("Error rates used for updating all changes ",tp_alpha," ",tp_beta)

    # check here if FP FN error rate diff is large. If yes then don't drop the subclone and continue to the next.
    #tp_cluster_genotype, tp_cluster_cells = clusteringResults(D_matrix, tp, tp_cluster_cells, tp_cluster_genotype, tp_alpha, tp_beta, tp_MR, cells, cluster, noOfiter)
    tp_cluster_genotype, tp_cluster_cells = clusteringResults(D_matrix, tp, tp_cluster_cells, tp_cluster_genotype, old_tp_FP, old_tp_FN, tp_MR, cells, cluster, noOfiter, tp_alpha, tp_beta)

    # Move to the next cluster since there are no more clusters or dropping this subclone will result in higher error rate
    if tp_cluster_genotype == "NA" or tp_cluster_cells == "NA":
        #print("Cluster not dropped ",cluster," in timepoint ",tp)
        tp_cluster_genotype = prev_tp_cluster_genotype
        tp_cluster_cells = prev_tp_cluster_cells
        return None, tp_cluster_genotype, tp_cluster_cells, old_tp_FP, old_tp_FN

    # Definitive mutation-preservation check.
    # Even if this cluster was not the sole carrier before pruning,
    # merging/reassignment must not cause any currently represented
    # mutation to disappear completely.
    if preserve_mutations:
        before_mutations = getPresentMutations(prev_tp_cluster_genotype)
        after_mutations = getPresentMutations(tp_cluster_genotype)
        lost_mutations = sorted(before_mutations - after_mutations)

        if lost_mutations:
            print(
                "MUTATION-PRESERVING MODE: pruning cluster",
                cluster,
                "in timepoint",
                tp,
                "rejected; would completely remove mutations",
                lost_mutations
            )

            tp_cluster_genotype = prev_tp_cluster_genotype
            tp_cluster_cells = prev_tp_cluster_cells

            return None, tp_cluster_genotype, tp_cluster_cells, old_tp_FP, old_tp_FN

    print(" Cluster genotype after merging ")
    printMutFromGenotype(tp_cluster_genotype[tp])

    new_tp_FP, new_tp_FN = timepointFPFN(D_matrix, tp_cluster_cells, tp_cluster_genotype)
    print(" New FP FN rates ",new_tp_FP," ",new_tp_FN)
    # If dropping cluster results in FP or FN > threshold then DON'T drop the cluster
    if check_FPFN(new_tp_FP, new_tp_FN, tp, tp_FPrate_threshold, tp_FNrate_threshold):
        print("Cluster not dropped ",cluster," in timepoint ",tp)
        tp_cluster_genotype = prev_tp_cluster_genotype
        tp_cluster_cells = prev_tp_cluster_cells
        #print("Timepoint assigned cells ",tp_cluster_cells)
        return None, tp_cluster_genotype, tp_cluster_cells, old_tp_FP, old_tp_FN

    return tp_cluster, tp_cluster_genotype, tp_cluster_cells, new_tp_FP, new_tp_FN

''' Select the tree with the highest probability. '''
# Input is D_matrix, alpha and beta in each timepoint, cluster probabilities in each timepoint, cells in each cluster, cluster genotypes
def selectOptimalTree(D_matrix, bnpc_cells_genotype, tp_alpha, tp_beta, tp_MR, tp_cluster_prob, tp_cluster_cells, tp_cluster_genotype, tpCells, k, tp_FPrate_threshold, tp_FNrate_threshold, max_parallel_extra,plotOp, bnpcRun, sample, preserve_mutations=False):
    # Get the initial Tree from the clustering results
    initialTree, initialCloneNode, tree_prob, initial_tp_FP, initial_tp_FN = getInitialTree(D_matrix, tp_cluster_prob, tp_cluster_cells, tp_cluster_genotype)

    finalTree = copy.deepcopy(initialTree) # initialize final tree
    noOfMutations = len(D_matrix[0])
    print("No of mutations ",noOfMutations)
    finalClusterCells = copy.deepcopy(tp_cluster_cells)
    finalClusterGen = copy.deepcopy(tp_cluster_genotype)
    final_backMut_edges = getBackMutCount(finalTree)

    #initial_tp_FP, initial_tp_FN = tp_alpha, tp_beta
    final_tp_FP, final_tp_FN = initial_tp_FP, initial_tp_FN
    print(" ==== Estimated FP FN before eliminating any spurious subclones ==== ")
    print(" Estimated FP ",final_tp_FP," Estimated FN ",final_tp_FN)
    finalCloneNode = initialCloneNode

    eliminated_clusters = [] # List to maintain the clusters eliminated in each iteration
    refresh_tp_cluster_prob = {}
    refresh_index = 0
    print("Initial Tree ===============================")
    printTree(finalTree)
    # plotTree(finalTree, {}, {}, "plots/initialTree.pdf", "Initial Tree", sample)
    print("Spurious subclones ",tp_cluster_prob)
    getNodeProb(tp_cluster_prob, finalCloneNode)
    noOfiter = 1
    for tp_cluster in tp_cluster_prob:
        print("Spurious subclones ",tp_cluster_prob)
        print("Eliminated clusters ",eliminated_clusters)
        print("Refreshed list ",refresh_tp_cluster_prob)
        # Since we refresh the cluster prob list we might come across clusters already eliminated so we prevent checking them.
        if tp_cluster in eliminated_clusters:
            continue
        print(" Iterations ============================== ",noOfiter)
        # If the cluster prob list is already updated then we iterate over the refresh list
        if refresh_tp_cluster_prob != {}:
            print("Refreshed list ",refresh_tp_cluster_prob)
            getNodeProb(refresh_tp_cluster_prob, clone_node)
            for tp_cluster in refresh_tp_cluster_prob:
                if tp_cluster in eliminated_clusters:
                    continue
                print("Inside refreshed list ")
                print(" Refreshed iterations ============================== ",noOfiter)
                #print("Timepoint cluster cells ",tp_cluster_cells)
                
                # Check for each cluster if it should be part of the Tree.
                cluster, tp_cluster_genotype, tp_cluster_cells, new_tp_FP, new_tp_FN = iterationSteps(D_matrix, tp_cluster, tp_cluster_genotype, tp_cluster_cells, tp_MR, noOfiter, tp_alpha, tp_beta, tp_FPrate_threshold, tp_FNrate_threshold, preserve_mutations)
                if cluster == None:
                    continue

                # Add the cluster eliminated
                eliminated_clusters.append(cluster)
                # Get the Tree
                Tree, clone_node = getTree(tp_cluster_cells, tp_cluster_genotype)
                print(" Before correction iteration ",noOfiter," Clone node ",clone_node)
                
                new_backMut_edges = getBackMutCount(Tree)
                
                # Get the new cluster prob
                new_tp_cluster_prob = calc_cluster_prob(tp_cluster_cells, tp_cluster_genotype, D_matrix, new_tp_FP, new_tp_FN, tp_MR)
                print("After correction Iteration ",noOfiter," Clone node ",clone_node)
                print("Tp_cluster_prob ",new_tp_cluster_prob)

                new_tree_prob = calculate_treeProb(new_tp_cluster_prob)
                print("Tree prob ",tree_prob," New tree prob ",new_tree_prob)
                noOfiter = noOfiter+1

                if new_tree_prob > tree_prob:
                    tree_prob = new_tree_prob
                    finalTree = copy.deepcopy(Tree)
                    finalCloneNode = clone_node
                    finalClusterCells = tp_cluster_cells
                    finalClusterGen = tp_cluster_genotype
                    final_backMut_edges = getBackMutCount(Tree)
                    final_tp_FP = new_tp_FP
                    final_tp_FN = new_tp_FN
                    print("FINAL TREE ITERATION INSIDE REFRESH ",noOfiter-1)
                    printTree(finalTree)
                    # Refresh the cluster prob list and continue checking clusters with this list
                    refresh_tp_cluster_prob = dict(sorted(new_tp_cluster_prob.items(), key=lambda item: item[1]))
                    break
                else:
                    print("=========== TREE NOT SELECTED =============")
                    printTree(Tree)
                    continue
        else:
            print("TP cluster gen ",tp_cluster_genotype)
            print("TP cluster cells ",tp_cluster_cells)
            print("Tp_cluster_prob ",tp_cluster_prob)

            # Check for each cluster if it should be part of the Tree.
            cluster, tp_cluster_genotype, tp_cluster_cells, new_tp_FP, new_tp_FN = iterationSteps(D_matrix, tp_cluster, tp_cluster_genotype, tp_cluster_cells, tp_MR, noOfiter, tp_alpha, tp_beta, tp_FPrate_threshold, tp_FNrate_threshold, preserve_mutations)
            if cluster == None:
                continue
            
            # Add the cluster eliminated
            eliminated_clusters.append(cluster)

            # Get the Tree 
            Tree, clone_node = getTree(tp_cluster_cells, tp_cluster_genotype)
            print(" Before correction iteration ",noOfiter," Clone node ",clone_node)
        
            new_backMut_edges = getBackMutCount(Tree)

            # Get the new cluster prob
            new_tp_cluster_prob = calc_cluster_prob(tp_cluster_cells, tp_cluster_genotype, D_matrix, new_tp_FP, new_tp_FN, tp_MR)
            print("After correction Iteration ",noOfiter," Clone node ",clone_node)
            print("Tp_cluster_prob ",new_tp_cluster_prob)

            new_tree_prob = calculate_treeProb(new_tp_cluster_prob)
            print("Tree prob ",tree_prob," New tree prob ",new_tree_prob)
            noOfiter = noOfiter+1
        
            if new_tree_prob > tree_prob:
                tree_prob = new_tree_prob
                finalTree = copy.deepcopy(Tree)
                finalCloneNode = clone_node
                finalClusterCells = tp_cluster_cells
                finalClusterGen = tp_cluster_genotype
                final_backMut_edges = getBackMutCount(Tree)
                final_tp_FP = new_tp_FP
                final_tp_FN = new_tp_FN
                print("FINAL TREE ITERATION FROM OUTER ",noOfiter)
                # Refresh the list to get only updated cluster probabilites
                refresh_tp_cluster_prob = dict(sorted(new_tp_cluster_prob.items(), key=lambda item: item[1]))
            else:
                print("=========== TREE NOT SELECTED =============")
                printTree(Tree)
                continue
    # the saved Tree will be the one with highest prob
    print(" ================= FINAL TREE =================== ")
    printTree(finalTree)
    #plotTree(finalTree, {}, {},"plots/dummy.pdf", "Final Tree ", sample)
    # Fix the parallel and back mutation. After correction the tp_cluster_genotype will change.
    finalTree = recheckFinalTree(finalTree) # In this method try to fix the parallel and back mutations
    # plotTree(finalTree, {}, {}, "plots/dummy.pdf", "Tree before parallel and back mutation correction", sample)
    treeBeforeCorrection = copy.deepcopy(finalTree)


    #finalTree, cg = correctParallelAndBackMut(D_matrix, finalClusterCells, finalClusterGen, tp_alpha, tp_beta, finalTree, finalCloneNode, k, plotOp, noOfiter, sample)
    #finalTree, cg = correctParallelAndBackMut(D_matrix, finalClusterCells, finalClusterGen, final_tp_FP, final_tp_FN, finalTree, finalCloneNode, k, plotOp, noOfiter, sample)
    finalTree, cg, selParallel, selBack = correctParallelAndBackMut(D_matrix, finalClusterCells, finalClusterGen,final_tp_FP, final_tp_FN,finalTree, finalCloneNode,k, plotOp, noOfiter, sample,max_parallel_extra)

    final_tp_FP, final_tp_FN = timepointFPFN(bnpc_cells_genotype, finalClusterCells, cg)
    print("Final FP after correction ",final_tp_FP," Final FN after correction ",final_tp_FN)
    #plotTree(finalTree, {}, {},"plots/"+plotOp+"/"+bnpcRun+"_Final_pbm.pdf", "Final Tree ", sample)
    finalTree = recheckFinalTree(finalTree)
    printTree(finalTree)
    plot_backMut_edges = getBackMutCount(finalTree)
    # plotTree(finalTree, selParallel, plot_backMut_edges, "plots/"+plotOp+"/"+bnpcRun+"_Final_ac.pdf", "Final Tree", sample)
    #plotTree(finalTree, {}, plot_backMut_edges,"plots/"+plotOp+"/"+bnpcRun+"_Final_ac.pdf", " ", sample)
    print("Final tree prob ",tree_prob)
    return finalTree, plot_backMut_edges, tree_prob, final_tp_FP, final_tp_FN, treeBeforeCorrection

# Input is cells in each timepoint file.
# Output a dictionary with timepoint as key and cells as value.
def cellTimepoints(tpCellsFile):
    tcf = open(tpCellsFile, 'r')
    tcf_lines = tcf.readline().rstrip('\n')
    tpCells = []
    while(tcf_lines != ""):
        #tp = (tcf_lines.split('\t')[0]).split('_')[1]
        cells = (tcf_lines.split('\t')[1]).split(';')
        print(cells)
        tpCells.append(cells)
        tcf_lines = tcf.readline().rstrip('\n')
    #print(tpCells)
    return tpCells

# Input is the D matrix file.
# Output is a list of the genotypes with indices corresponding to cells.
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

# Input is the BnpC's log file.
# Output is the FP FN values at each timepoint.
def readFPFNvalues(logF):
    tp_alpha = {}
    tp_beta = {}
    epsilon = float(1e-6)
    print(" Epsilon ",epsilon)
    for i in range(len(logF)):
        lfile = open(logF[i],'r')
        lf_line = lfile.readline().rstrip('\n')
        lf_line = lfile.readline().rstrip('\n')
        lf_arr = lf_line.split('\t')
        fn_rate = lf_arr[3]
        fp_rate = lf_arr[5]
        tp_alpha[i] = float(fp_rate) + epsilon
        tp_beta[i] = float(fn_rate) + epsilon
    print(" TP alpha ",tp_alpha)
    print(" TP Beta ",tp_beta)
    return tp_alpha, tp_beta

''' Get the missing rate for each time point. '''
# Input is cells in each time point and noisy D matrix
def calcFPFN(D,G):
    FN = 0
    FP = 0
    TP = 0
    TN = 0
    #print(" Mut ",len(D[0]))
    #print(" Mut ",len(G[0]))
    for i in range(len(G)):
        for j in range(len(G[i])):
            if D[i][j] == 0 and G[i][j] == 1:
                FN = FN + 1
            if D[i][j] == 1 and G[i][j] == 0:
                FP = FP + 1
            if D[i][j] == 1 and G[i][j] == 1:
                TP = TP + 1
            if D[i][j] == 0 and G[i][j] == 0:
                TN = TN + 1
    print("FP ",FP,"TN ",TN,"FN ",FN,"TP ",TP)
    if FP+TN == 0 or FP == 0:
        FPR = 1e-06
    else:
        FPR = FP/(FP+TN)

    if FN+TP == 0 or FN == 0:
        FNR = 1e-06
    else:
        FNR = FN/(FN+TP)
    return FPR, FNR

''' Calculate estimated FP FN rates for each timepoint after getting the final tree. '''
def timepointFPFN(D_matrix, finalClusterCells, finalClusterGen):
    print("Inside timepointFPFN ")
    print(finalClusterCells)
    print(finalClusterGen)
    tp_FP = {}
    tp_FN = {}
    for tp, cluster in finalClusterCells.items():
        tp_cluster_G = []
        tp_D = []
        tp_cells = []

        for c in cluster:
            cells = finalClusterCells[tp][c]
            tp_cells.extend(cells)
            for cell in cells:
                tp_cluster_G.append(finalClusterGen[tp][c])
       
        print("Tp cluster cells ",len(tp_cells))
        #print(len(tp_cluster_G))
        for c in tp_cells:
            tp_D.append(D_matrix[c])
        #print(len(tp_D))

        FPR, FNR = calcFPFN(tp_D, tp_cluster_G)
        tp_FP[tp] = FPR
        tp_FN[tp] = FNR
        print(" Timepoint ",tp," FP ",FPR," FN ",FNR)
    return tp_FP, tp_FN

''' Get the initial clustering results. '''
def intialClusterResults(tpClusters, allClusters, tpCells, bnpc_cells_genotype, D_matrix, tp_alpha, tp_beta, tp_MR):
    tp_cluster_cells = getNewClusters(tpClusters, allClusters, tpCells)
    # For initial clustering get the updated cluster genotype based on max likelihood
    print(" Initial clustering error rates used ",tp_alpha," ",tp_beta)
    tp_reassignedCells, tp_updatedCG = correctClustering(bnpc_cells_genotype, D_matrix, tp_cluster_cells, tp_alpha, tp_beta, {}, True)
    tp_cluster_prob = calc_cluster_prob(tp_reassignedCells, tp_updatedCG, D_matrix, tp_alpha, tp_beta, tp_MR)
    print(" Initial Timepoint cluster prob ",tp_cluster_prob)
    sorted_cluster_prob = dict(sorted(tp_cluster_prob.items(), key=lambda item:item[1]))
    print(sorted_cluster_prob)
    return sorted_cluster_prob, tp_reassignedCells, tp_updatedCG

