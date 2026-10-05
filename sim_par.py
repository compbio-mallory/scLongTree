#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
@authors: Rituparna Khan
Contacting email: fan@cs.fsu.edu
"""

import sys
import argparse
import numpy as np
import random
import copy
from collections import OrderedDict

class Edge():
    def __init__(self, p, c):
        self.p = p
        self.c = c

class Node():
    def __init__(self, e, p, c):
        # only the edge above this node is listed here
        self.e = e
        self.p = p
        # this node may have multiple children, c is an array
        self.c = c

# return a n*m matrix with zero entries
def init(n, m):
    ret = []
    for i in range(n):
        ret.append([0]*m)
        for j in range(m):
            ret[i][j] = 0
    return ret
   
def print_matrix(M, n, matrix_file):
    matrix_f = open(matrix_file, "w")
    for i in range(n):
        str_ = [str(j) for j in M[i]]
        print("\t".join(str_), file = matrix_f) 
    matrix_f.close()

# given a leaf ID, and an edge dict that has edge ID as the key, .c and .p as the child and parent node ID, return all the edges above this leaf ID in array. This function is used for doublets. May need it later.
#def retrieve_edges(leafID, n_dict):
#    p = n_dict[leafID].p
#    e = n_dict[leafID].e
#    e_array = [e]
#    while p != "-1":
#        ID = p
#        p = n_dict[ID].p
#        e = n_dict[ID].e
#        e_array.append(e)
#    return e_array
 
def get_MutationNo(lam_value):
    mut_num = 0
    while mut_num == 0: # Prevent mutation number from being 0.
        mut_num = np.random.poisson(lam=lam_value)
    return mut_num

def print_f(out_dict, out_file):
    out_f = open(out_file, 'w')
    for key,value in out_dict.items():
        print("\t".join([str(key),str(value)]), file = out_f)
    out_f.close()

# Distribute mutations based on branch length and Poisson distribution. Replaced the mutation constant (mc) by number of mutations.
def distribute_mutations(tree_f, mc, out_f):
    tree_dict = {}
    file = open(tree_f,"r")
    first_line = file.readline().rstrip('\n') # ignore it as this is header
    line = file.readline().rstrip('\n')
    mut_array = []
    mut_dict = {} # Has mutations saved according to node ID. The mutations are distributed above the edge of the node ID.
    while(line != ""):
        line_a = line.split('\t')
        nodeID = int(line_a[1])
        pID = line_a[2]
        edgeL = float(line_a[3])
        mut = line_a[6]
        tree_dict[nodeID] = str(pID)+";"+str(edgeL)+";"+str(mut)
        line = file.readline().rstrip('\n')

    #print(tree_dict)
    prev_mut = -1
    for nodeId in tree_dict.keys():
        if nodeId == 0:
            continue
        val_list = tree_dict[nodeId].split(';')
        #print(" Val list ",val_list)
        pID = int(val_list[0])
        edgeL = val_list[1]
        mut = val_list[2]
        print(" Node ",nodeId," Parent ID ",pID," Mutation ",mut)

        if mut == 'True':
            lam_value = float(edgeL) * float(mc)
            #print(" Lambda value ",lam_value)
            mut_num = get_MutationNo(lam_value)
            print(" No of mutations ",mut_num)

            #mut_num = round(float(edgeL) * total_mut)
            mut_IDs = str(prev_mut + 1)
            for j in range(prev_mut + 2, prev_mut + mut_num + 1):
                mut_IDs = mut_IDs + ";" + str(j)
            prev_mut = prev_mut + mut_num
            #print(" No of mutations ",mut_num," mutation ID ",mut_IDs)
            mut_dict[nodeId] = mut_IDs
        else:
            # Save mutations of prev parent node and re-assign it for the child node
            child_mut_IDs = mut_dict[pID] 
            mut_dict[nodeId] = child_mut_IDs
            #print(" Mutations of node ",nodeId," ",mut_IDs)
    print(mut_dict)
    #print(" Total no of mutations ",prev_mut)
    print_f(mut_dict, out_f) # Saves the mutations in a file
    return mut_dict, prev_mut+1
 
# Function to distribute SNVcells at each timepoint.
def distribute_SNVcells(cell_n, tree_f, out_f):
    file = open(tree_f,"r")
    first_line = file.readline().rstrip('\n') # ignore it as this is header
    line = file.readline().rstrip('\n')
    tree_dict = OrderedDict() # Has the timepoints along with the nodes and their intervals. OrderedDict will help maintain the order of time.
    while(line != ""):
        line_a = line.split('\t')
        timepoint = float(line_a[0])
        leaf_val = int(line_a[5])
        # Commenting out leaf_val == -1 because we want to sequence cells on dead node.
        if timepoint == 1.0 or not timepoint.is_integer():
            print(timepoint)
            line = file.readline().rstrip('\n')
            continue
        
        nodeID = int(line_a[1])
        perc = float(line_a[4]) # Interval value
        #print(" Timepoint ",timepoint," NodeID ",nodeID," Perc ",perc)
        if timepoint in tree_dict:
            temp_list = tree_dict[timepoint]
            temp_list.append(str(nodeID)+";"+str(perc))
            tree_dict[timepoint] = temp_list
        else:
            temp_list = []
            temp_list.append(str(nodeID)+";"+str(perc))
            tree_dict[timepoint] = temp_list
        line = file.readline().rstrip('\n')
    print(tree_dict)
    # Loop over the tree_dict with a value of i to access the cell_n array. And assign the cells according to that.
    i = 0 # Use this as an index to access cell numbers from cell_n
    SNVcell_dict = {}
    prev_cell = -1
    for key, value in tree_dict.items():
        cell_num = int(cell_n[i])
        total_cells = 0
        print(" No. of cells ",cell_num)
        for nodes in value:
            nodeID = (nodes.split(';'))[0]
            SNVcellP = (nodes.split(';'))[1]
            SNVcell_num = round(cell_num * float(SNVcellP))
            print("Node ID ",nodeID," SNV cell num ",SNVcell_num)

            if SNVcell_num == 0: # If SNVcell_num becomes negative then take 1 cell from previous node and assign it here. 
                prev_node = int(nodeID) - 1
                SNVcell_dict_key = str(int(key))+'_'+str(prev_node)
                #if SNVcell_dict_key in SNVcell_dict: # To prevent key error and prevent looking in previous timepoint.
                prev_node_cells = SNVcell_dict[SNVcell_dict_key].split(";")
                SNVcell_dict[SNVcell_dict_key] = ";".join(prev_node_cells[:-1]) # Update prev node's cellIDs by removing one from them
                SNVcell_dict[str(int(key))+'_'+nodeID] = prev_node_cells[len(prev_node_cells)-1] # And assign that cellID from the prev node to the current node. Getting the last cellID to maintain the order
                SNVcell_num = 1
                #else:
                #SNVcell_num = 1 # Sometimes due to very small percentage this no. was 0 but we still want to assign atleast 1 cell to the node.

            total_cells = total_cells+SNVcell_num
            #print(" After assigning total cells ",total_cells)

            if total_cells > cell_num: # check that cell assignment does not surpass the total available cells.
                diff = total_cells - cell_num
                SNVcell_num = SNVcell_num - diff
            
            print(" NodeID ",nodeID," Cell perc ",SNVcellP," cell no ",SNVcell_num)

            if SNVcell_num > 0:
                cell_IDs = str(prev_cell + 1)
                for j in range(prev_cell + 2, prev_cell + SNVcell_num + 1):
                    cell_IDs = cell_IDs + ";" + str(j)
                prev_cell = prev_cell + SNVcell_num
                print(" Prev cell ",prev_cell)
                SNVcell_dict[str(int(key))+'_'+nodeID] = cell_IDs
                last_node = nodeID
                #print(SNVcell_dict)

        #print(" Total cells ",total_cells," in tp ",key)
        if total_cells < cell_num:
            diff = cell_num - total_cells
            cellIDs = SNVcell_dict[str(int(key))+'_'+last_node]
            #print(" Before updating ",cellIDs)
            cellIDs_list = SNVcell_dict[str(int(key))+'_'+last_node].split(";")
            cellIDs_len = len(cellIDs_list)
            lastCellID = cellIDs_list[cellIDs_len-1]
            #print(" Last cell ID ",lastCellID," Diff ",diff)

            for j in range(int(lastCellID) + 1, int(lastCellID) + diff + 1):
                cellIDs = cellIDs + ";" + str(j)
            prev_cell = int(lastCellID) + diff
            #print(" After updating ",cellIDs)
            SNVcell_dict[str(int(key))+'_'+last_node] = cellIDs
        i = i+1

    print(SNVcell_dict)
    print_f(SNVcell_dict, out_f) # Saves the cells in a file
    print(" Final total no of cells ",prev_cell+1)
    return SNVcell_dict,prev_cell+1

def add_missing(G, n, m, missingP):
    # in G (n*m) matrix, select missingP * n * m entries and flip them to 3. 
    missing_arr = random.sample(range(n * m), int(n * m * missingP))
    for i in missing_arr:
        row = int(i / m)
        column = i % m
        G[row][column] = 3
    return G


def count_total_value(M, n, m, value):
    total = 0
    for i in range(n):
        for j in range(m):
            if M[i][j] == value:
                total = total + 1
    return total


def add_FPFNs(D, n, m, alpha, beta, total_zeros, total_ones):
    FP_arr = random.sample(range(total_zeros), int(total_zeros * alpha))
    FN_arr = random.sample(range(total_ones), int(total_ones * beta))
    D_ = copy.deepcopy(D)
    index_neg = 0
    index_pos = 0
    for i in range(n):
        for j in range(m):
            if D[i][j] == 0:
                if index_neg in FP_arr:
                    # flip it
                    D_[i][j] = 1
                index_neg = index_neg + 1
            elif D[i][j] == 1:
                if index_pos in FN_arr:
                    # flip it
                    D_[i][j] = 0
                index_pos = index_pos + 1
    return D_
    

# from tree_f we know which cells (on leaf nodes) have which mutations (on edges).
def mutation_matrix(mut_dict, SNVcell_dict, tree_f, n, m, missingP, alpha, beta, G_matrix_file, D_matrix_file, D_miss_file):

    file = open(tree_f, "r")
    first_line = file.readline().rstrip('\n')
    line = file.readline().rstrip('\n')

    # node_dict: a dict that has the child node ID as the key, .e and .p as the values, the edge above it and the parent node ID, respectively
    # edge_dict: a dict that has the edge ID as the key, .c and .p as the values, the child and parent node IDs on the two ends, respectivley
    edge_dict = {} 
    node_dict = {}
    # c is a temporary data structure for retrieving all children nodes for all parent nodes for node_dict
    c_dict = {}
    while(line != ""):
        line_a = line.split('\t')
        # line_a[2] is parent, line_a[1] is child
        edge_dict[line_a[1]] = Edge(line_a[1], line_a[2])
        node_dict[line_a[1]] = Node(line_a[1], line_a[2], "NA")
        # record the children for each parent to be added to node_dict
        if line_a[2] not in c_dict.keys():
            c_dict[line_a[2]] = line_a[1]
        else:
            c_dict[line_a[2]] = c_dict[line_a[2]] + ";" + line_a[1]
        line = file.readline().rstrip('\n')

    file.close()

    print(" c_dict ",c_dict)

    for i in node_dict.keys():
        if i in c_dict.keys():
            node_dict[i].c = c_dict[i]
        else:
            # this is a leaf or dead
            node_dict[i].c = "NA"

    for k, v in node_dict.items():
        print(k," --> Edge ",v.e," PID ",v.p," Child ",v.c)

    # a dict from SNV cell ID to the mutation ID array
    SNVcell_mut = {}
    for node, cells in SNVcell_dict.items():
        #timepoint = (node.split("_"))[0]
        nodeID = int((node.split("_"))[1])
        pID = node_dict[str(nodeID)].p
        mutID_array = []
        for mutID in mut_dict[nodeID].split(";"):
            mutID_array.append(mutID)

        for mutID in mut_dict[int(pID)].split(";"): # Add mutations from parent node to children.
            if mutID in mutID_array:
                continue
            mutID_array.append(mutID)
            c_mutID = mut_dict[nodeID]
            c_mutID = c_mutID + ";" + mutID
            mut_dict[nodeID] = c_mutID

        for cell in cells.split(";"):
            SNVcell_mut[cell] = copy.deepcopy(mutID_array)

    #print(SNVcell_mut)
    # construct G matrix. n*m: n(cell), m(mut)
    G = init(n, m)
    for cellID in range(n):
        if str(cellID) in SNVcell_mut.keys():
            mut_array_ = SNVcell_mut[str(cellID)]
            for mutID in range(m):
                #print(cellID," ",mut_array_," ",mutID)
                if str(mutID) in mut_array_:
                    G[cellID][mutID] = 1

    print_matrix(G, n, G_matrix_file)
    total_zeros = count_total_value(G, n, m, 0)
    total_ones = count_total_value(G, n, m, 1)
    print("total zeros for G: " + str(total_zeros) + "; total ones for G: " + str(total_ones))

    # Step a. now make missing data to produce D 
    D_miss = add_missing(G, n, m, missingP) 
    print_matrix(D_miss, n, D_miss_file)

    # Step b, c. flip 0's by \alpha and 1's by \beta 
    total_zeros = count_total_value(D_miss, n, m, 0)
    total_ones = count_total_value(D_miss, n, m, 1)
    print("total zeros after missing: " + str(total_zeros) + "; total ones after missing: " + str(total_ones))
    D_miss_FP_FN = add_FPFNs(D_miss, n, m, alpha, beta, total_zeros, total_ones)
    total_zeros = count_total_value(D_miss_FP_FN, n, m, 0)
    total_ones = count_total_value(D_miss_FP_FN, n, m, 1)
    print("total zeros after FPFN: " + str(total_zeros) + "; total ones after FPFN: " + str(total_ones))

    # store the D matrix
    print_matrix(D_miss_FP_FN, n, D_matrix_file)

# Use this function to just flip the value of G for the given FP, FN and missing rate
def change_FPFN_MR(G, alpha, beta, missingP, D_matrix_file):
    # First read the G file and save it for modification
    print(" FP ",alpha," FN ",beta," MR ",missingP)
    GFile = open(G,'r')
    GFile_line = GFile.readline().rstrip('\n')
    G_matrix = []
    while GFile_line != "":
        gline = GFile_line.split('\t')
        gline = [int(i) for i in gline]
        G_matrix.append(gline)
        GFile_line = GFile.readline().rstrip('\n')

    n = len(G_matrix) # no. of cells
    m = len(G_matrix[0]) # no. of mutations
    print(" No. of cells ",n)
    print(" No. of mutations ",m)
    total_zeros = count_total_value(G_matrix, n, m, 0)
    total_ones = count_total_value(G_matrix, n, m, 1)
    print("total zeros for G: " + str(total_zeros) + "; total ones for G: " + str(total_ones))

    # Step a. now make missing data to produce D 
    D_miss = add_missing(G_matrix, n, m, missingP)
    #print_matrix(D_miss, n, D_miss_file)

    # Step b, c. flip 0's by \alpha and 1's by \beta 
    total_zeros = count_total_value(D_miss, n, m, 0)
    total_ones = count_total_value(D_miss, n, m, 1)
    print("total zeros after missing: " + str(total_zeros) + "; total ones after missing: " + str(total_ones))
    D_miss_FP_FN = add_FPFNs(D_miss, n, m, alpha, beta, total_zeros, total_ones)
    total_zeros = count_total_value(D_miss_FP_FN, n, m, 0)
    total_ones = count_total_value(D_miss_FP_FN, n, m, 1)
    print("total zeros after FPFN: " + str(total_zeros) + "; total ones after FPFN: " + str(total_ones))
    
    # store the D matrix
    print_matrix(D_miss_FP_FN, n, D_matrix_file)

''' Vary the FP, FN and MR for each timepoints. '''
def vary_FPFN_MR_timepoint(G, cellsTp_file, assigned_FP, assigned_FN, assigned_MR, D_matrix_file):
    print(" Inside varying FP, FN and MR ")
    print(" FP rates ",assigned_FP)
    print(" FN rates ",assigned_FN)
    print(" MR rates ",assigned_MR)
    timepoint_cells = {} # Get the cells or indices for each timepoint
    cells_file = open(cellsTp_file,'r')
    cellsFileLine = cells_file.readline().rstrip('\n')

    while(cellsFileLine != ""):
        timepoint = ((cellsFileLine.split('\t')[0]).split('_'))[0]
        cells = (cellsFileLine.split('\t')[1]).split(';')
        #print(" Timepoint ",timepoint," Cells ",cells)
        if timepoint in timepoint_cells:
            cells_list = timepoint_cells[timepoint]
            cells_list.extend(cells)
            timepoint_cells[timepoint] = cells_list
        else:
            timepoint_cells[timepoint] = cells

        cellsFileLine = cells_file.readline().rstrip('\n')
    #print(timepoint_cells)
    
    # Read the G matrix from the file
    GFile = open(G,'r')
    GFile_line = GFile.readline().rstrip('\n')
    G_matrix = []
    while GFile_line != "":
        gline = GFile_line.split('\t')
        gline = [int(i) for i in gline]
        G_matrix.append(gline)
        GFile_line = GFile.readline().rstrip('\n')

    n = len(G_matrix) # no. of cells
    m = len(G_matrix[0]) # no. of mutations
    print(" No. of cells ",n)
    print(" No. of mutations ",m)
    total_zeros = count_total_value(G_matrix, n, m, 0)
    total_ones = count_total_value(G_matrix, n, m, 1)
    print("total zeros for G: " + str(total_zeros) + "; total ones for G: " + str(total_ones))

    index = 0 # Use this index to get the FP, FN and MR values from their list
    final_D_matrix = [] # Combine the matrices from different timepoints
    for timepoint, cells in timepoint_cells.items():
        G_timepoint = []
        # each cell will be the index for G matrix here.
        for cell in cells:
            G_timepoint.append(G_matrix[int(cell)])
        print(len(G_timepoint))
        FP_rate = assigned_FP[index]
        FN_rate = assigned_FN[index]
        MR_rate = assigned_MR[index]
        print(" FP rate ",FP_rate," FN rate ",FN_rate," MR rate ",MR_rate)
       
        n = len(G_timepoint) # no. of cells 
        m = len(G_timepoint[0]) # no. of mutations
        print(" Total cells ",n," Total mutations ",m)
        # Step a. now make missing data to produce D
        D_miss = add_missing(G_timepoint, n, m, MR_rate)
        #print_matrix(D_miss, n, D_miss_file)

        # Step b, c. flip 0's by \alpha and 1's by \beta
        total_zeros = count_total_value(D_miss, n, m, 0)
        total_ones = count_total_value(D_miss, n, m, 1)
        print("total zeros after missing: " + str(total_zeros) + "; total ones after missing: " + str(total_ones))
        D_miss_FP_FN = add_FPFNs(D_miss, n, m, FP_rate, FN_rate, total_zeros, total_ones)
        total_zeros = count_total_value(D_miss_FP_FN, n, m, 0)
        total_ones = count_total_value(D_miss_FP_FN, n, m, 1)
        print("total zeros after FPFN: " + str(total_zeros) + "; total ones after FPFN: " + str(total_ones))
        print(" D matrix len ",len(D_miss_FP_FN)," ",timepoint)
        final_D_matrix.extend(D_miss_FP_FN)
        index = index+1
    print(len(final_D_matrix))
    #print(final_D_matrix)
    print_matrix(final_D_matrix, len(final_D_matrix), D_matrix_file)

# Uniformly sample cells from the following array for each timepoint
def sampleCells(noOfTimepoints):
   cells = [100,300,600,1000]
   assigned_cells = []
   for i in range(1,int(noOfTimepoints)):
       assigned_cells.append(random.choice(cells))
   return assigned_cells

''' Uniformly sample FP, FN and MR for each timepoints. '''
def sampleFPFNMR(noOfTimepoints):
    fp_rate = [0.001, 0.01, 0.03, 0.05]
    fn_rate = [0.1, 0.2, 0.3, 0.4]
    missing_Rate = [0.2, 0.3]
    assigned_FP = []
    assigned_FN = []
    assigned_MR = []
    for i in range(1,int(noOfTimepoints)):
        assigned_FP.append(random.choice(fp_rate))
        assigned_FN.append(random.choice(fn_rate))
        assigned_MR.append(random.choice(missing_Rate))
    return assigned_FP, assigned_FN, assigned_MR

''' Vary with less no. of FP, FN and MR. '''
def lessVariedFPFNMR():
    fp_rate = [0.01, 0.03]
    fn_rate = [0.2, 0.3]
    missing_Rate = [0.2, 0.2]
    return fp_rate, fn_rate, missing_Rate

# beginning of the program
if len(sys.argv) <= 1:
    print("""
    This generates the mutation matrix with the ground truth data. 
    Usage: python sim_par.py -a [alpha] -b [beta] -m [missing-rate] -c [num-cells] -mc [mut_const] -f [input-tree-file] -P [prefix-output-files]
        -a (--alpha)        False positive rate. [0.01]
        -b (--beta)         False negative rate. [0.2]
        -m (--missing_rate) Missing rate in G. [0.2]
        -c (--num_cells)    Total number of SNV cells at each timepoint. [100,500]
        -mc (--mut_const)    Mutation constant for Poisson distribution. [0.2]
        -f (--tree_file)    The input tree structure file. ["NA"]
        -P (--prefix)       Prefix of output files. 
    """)
    sys.exit(0)

parser = argparse.ArgumentParser(description='This script iterates all parameters to generate simulated data.')
parser.add_argument('-a', '--alpha', default=0.01)
parser.add_argument('-b', '--beta', default=0.2)
parser.add_argument('-m', '--missing_rate', default=0.2)
#parser.add_argument('-c', '--num_cells', nargs="*", default=[100,500])
#parser.add_argument('-c', '--num_cells', nargs="*", default=[500,1000,1500])
parser.add_argument('-t','--timepoints', default=3)
parser.add_argument('-FPFN','--fpfn', default=False)
parser.add_argument('-mc', '--mut_const', default=1.7)
#parser.add_argument('-mut', '--mutations', default=50)
#parser.add_argument('-e', '--doublet-rate', default=0)
parser.add_argument('-f', '--tree_file', default="NA")
parser.add_argument('-P', '--prefix', default="NA")
parser.add_argument('-G', '--Gmatrix', default="NA")
parser.add_argument('-cellstp', '--cellstp', default="NA")

args = parser.parse_args()
alpha = float(args.alpha)
beta = float(args.beta)
missing = float(args.missing_rate)
#cell_n = args.num_cells
#cell_n = sampleCells(args.timepoints)
#print(" Assigned cells ",cell_n)
mut_const = float(args.mut_const)
#total_mutations = float(args.mutations)
#eta = float(args.doublet_rate)
tree_f = args.tree_file
prefix = args.prefix

print("For file --> ",tree_f)

if alpha == 0.001 or alpha == 0.05 or beta == 0.1 or beta == 0.3 or beta == 0.4 or missing == 0.3:
    Gfile = args.Gmatrix
    change_FPFN_MR(Gfile, alpha, beta, missing, prefix + ".D.csv")
elif args.fpfn == "True":
    Gfile = args.Gmatrix
    cellsTp_file = args.cellstp
    #assigned_FP, assigned_FN, assigned_MR = sampleFPFNMR(args.timepoints)
    assigned_FP, assigned_FN, assigned_MR = lessVariedFPFNMR()
    # Use the following method to vary FP, FN and MR for each timepoint
    vary_FPFN_MR_timepoint(Gfile, cellsTp_file, assigned_FP, assigned_FN, assigned_MR, prefix + ".D.csv") 
else:
    cell_n = sampleCells(args.timepoints)
    print(" Assigned cells ",cell_n)
    # First, distribute the mutations on the edges. Read tree_f output to prefix + ".mut.csv" 
    mut_dict, mut_n = distribute_mutations(tree_f, mut_const, prefix + ".mut.csv")
    #mut_dict, mut_n = distribute_mutations(tree_f, total_mutations, prefix + ".mut.csv")

    # Second, distribute the cells on the leaf nodes. Read tree_f output to prefix + ".SNVcell.csv" 
    SNVcell_dict, cell_num = distribute_SNVcells(cell_n, tree_f, prefix + ".SNVcell.csv")

    print(" Total cells ",cell_num," total mutations ",mut_n)

    # Last, make mutation matrices G and D. 
    mutation_matrix(mut_dict, SNVcell_dict, tree_f, cell_num, mut_n, missing, alpha, beta, prefix + ".G.csv", prefix + ".D.csv", prefix + ".miss.csv")

    print("============================================================================================================")
