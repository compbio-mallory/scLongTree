#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
@authors: Rituparna Khan
Contacting email: fan@cs.fsu.edu
"""

import sys
import argparse
import numpy as np

# tree is an array of treenode
class treenode():
    def __init__(self, id_):
        # root's ID is 0, root has two children
        self.id = id_
        self.tuple=[]
        # the edge is the one above the node; edge ID is the same as node id. 
        self.edge_length = 0
        # root's parent ID is -1
        self.parentID = -1
        # root's depth is 1
        # Depth is the time point here. For multisplit nodes it will be in between. For e.g, in between t1 and t2 it will be t1.5.
        self.depth_ = -1
        self.perc = -1
        self.if_leaf = 1 # 1 indicates it is a leaf, 0 indicates not a leaf and -1 indicates dead.
        self.mut = True
    def setTuple(self, a, b):
        self.tuple = [a, b]
    def getTuple(self):
        return self.tuple
    def getID(self):
        return self.id
    def getDepth(self):
        return self.depth_
    def getPerc(self):
        return self.perc
    def getEdgeLength(self):
        return self.edge_length
        
# store the tree by edges. Add timepoints and timepoints in between of multisplits.
def save_tree(T, out_f):
    f = open(out_f, "w")
    f.write("\t".join(['TP','ID','PID','EL','Perc','Leaf','Mut'])+ "\n")
    for i in range(len(T)):
        #leaf = 0
        #if T[i].if_leaf:
        #    leaf = 1
        # ID is the node ID
        #print(" Tuple values ",T[i].tuple[0]," ",T[i].tuple[1])
        f.write("\t".join([str(T[i].depth_), str(T[i].id), str(T[i].parentID), str(T[i].edge_length), str(T[i].perc), str(T[i].if_leaf), str(T[i].mut)]) + "\n")
        #f.write("\t".join([str(T[i].depth_), str(T[i].id), str(T[i].parentID), str(T[i].edge_length), str(T[i].tuple[1] - T[i].tuple[0]), str(T[i].if_leaf), str(T[i].mut)]) + "\n")
    f.close()

# Save the tree as a dict
def get_TreeDict(T):
    tree_dict = {} # timepoint -> (index, nodeID, perc)
    for i in range(len(T)):
        timepoint = T[i].depth_
        # we don't consider root, multisplit nodes and dead nodes. Removing dead node condn below because dead nodes can have cells.
        if timepoint == 1 or not timepoint.is_integer():
            continue
        if timepoint in tree_dict:
            temp_list = tree_dict[timepoint]
            temp_list.append(str(T[i].id)+';'+str(T[i].perc))
            tree_dict[timepoint] = temp_list
        else:
            temp_list = []
            temp_list.append(str(T[i].id)+';'+str(T[i].perc))
            tree_dict[timepoint] = temp_list
    return tree_dict

# Normalize the perc and reassign the intervals of the nodes in each timepoint.
def update_tpNodes_interval(tree_dict,T):
    for tp, nodes in tree_dict.items():
        nodes_list = tree_dict[tp]
        node_perc_dict = {}
        for node in nodes_list:
            nodeID = int((node.split(';'))[0])
            perc = float((node.split(';'))[1])
            node_perc_dict[nodeID] = perc
    
        #print(" Before normalizing ",node_perc_dict.values())
        summation = 0
        for t in node_perc_dict.values():
            summation += t
        #print(" Summation ",summation)
        for key, value in node_perc_dict.items():
            updated_perc = float(value)/float(summation)
            node_perc_dict[key] = updated_perc
            T[key].perc = updated_perc   
        #print(" Updated perc ",node_perc_dict)
    return T

# Generalizes to get the depth and edge length of multisplit timepoints and normal timepoints.
def get_depth_edgeLen(node_depth, parent_v, time_values):
    if node_depth.is_integer():
        depth = node_depth + 1
        t1_v = int(time_values[int(node_depth)-1])
        t2_v = int(time_values[int(node_depth)])
        edge_length = abs(t1_v - t2_v) # Assign edge length
        return depth, edge_length
    else:
        depth = node_depth + 0.5
        t1_v = int(time_values[int(node_depth)-1]) # Get the time from prev timepoint
        #print(" Time point ",node.depth_," t1_v ",t1_v)
        t2_v = int(time_values[int(node_depth)]) # Get the time from next timepoint
        #print(" Time point ",depth," t2_v ",t2_v)
        diff_v = abs(t1_v - t2_v)
        #print(" Node ",i," edgeL ",parent_v," Node number ",node_number)
        edge_length = abs(diff_v - parent_v) # Assign edge length
        return depth, edge_length

# Calculate the edge length for multisplit nodes
def get_multiSplit_edgeLength(t1,t2):
    ti = np.random.beta(2,2,2) # Here Alpha, Beta is 2.
    diff_t = abs(int(t2)-int(t1))
    t_1 = diff_t * ti[0]
    t_2 = diff_t * ti[1]
    print(" Multisplit Edges ",t_1, t_2)
    return t_1, t_2

# Checks if all nodes at a timepoint are dead
def checkIfAllNodes_Dead(tp_nodes, Tree):
    noOfNodes = len(tp_nodes)
    deadNodes = 0
    for node in tp_nodes:
        if Tree[node].if_leaf == -1:
            deadNodes = deadNodes+1
    if deadNodes == noOfNodes:
        return True

# First get the root node and split it into two unobserved nodes
def initializeTree(Beta, Alpha, timepoints, time_values, v1, v2, u, theta):
    Tree = []
    Tree.append(treenode(0))
    Tree[0].perc = 1
    Tree[0].parentID = -1
    Tree[0].edge_length = 0 # Root here still have edge length 0
    Tree[0].if_leaf = 0
    Tree[0].depth_ = 1 # Root is at timepoint t = 1
    Tree[0].tuple=[0,1]

    # The multisplit happens here based on the values of v1, v2, u. A root node will always have one multisplit.
    Tree.append(treenode(1))
    Tree.append(treenode(2))

    # set parent ID
    Tree[1].parentID = 0
    Tree[2].parentID = 0
    # set Edge length of children
    Tree[1].edge_length, Tree[2].edge_length = get_multiSplit_edgeLength(time_values[0],time_values[1])
    print(Tree[1].edge_length, Tree[2].edge_length)

    # set perc
    Bi = np.random.beta(float(Alpha+1),float(Beta+1),1) # Beta distribution for the nodes.
    Tree[1].perc = Bi[0]
    Tree[2].perc = 1 - Bi[0]

    Tree[1].if_leaf = 1 # This condition required to split the nodes further
    Tree[2].if_leaf = 1
    # set depth
    Tree[1].depth_ = 1.5
    Tree[2].depth_ = 1.5

    Tree[1].tuple=[0,Bi[0]]
    Tree[2].tuple=[Bi[0],1]

    tp_nodes = {} # This dictionary will store the nodes at any timepoint.
    tp_nodes[1.5] = [1,2]
    node_number = 2

    # ===== Multisplitting of root node ends here in t1. =====

    tp_depth = 1.5 # tp_depth represents the current timepoint. At this stage, we finished for root in t = 1.
    #node_oneSplit = set()
    tp_multiSplit = set()
    return Tree, tp_depth, tp_nodes, node_number, tp_multiSplit

def gen_tree(Beta, Alpha, timepoints, time_values, v1, v2, u, theta):
   
    # add a root (node 0) to the tree
    # edge length (there are at most 2*n - 1))
    #            | CN0
    #          node 0
    #        / CN1   \ CN2
    #    node 1    node 2

    #Contructing the phylogeny
    Tree, tp_depth, tp_nodes, node_number, tp_multiSplit = initializeTree(Beta, Alpha, timepoints, time_values, v1, v2, u, theta)

    while tp_depth < float(timepoints): # Replace the tree_W with depth condition
        print(" Depth of tree ",tp_depth)
        tp_nodes_list = tp_nodes.get(tp_depth)
        print(tp_nodes)

        # Before proceeding any further check if all nodes in a timepoint are already dead before final timepoint. If yes, then restart the process.
        if checkIfAllNodes_Dead(tp_nodes_list, Tree):
            print(" ENCOUNTERED ALL DEAD NODES AT A TIMEPOINT!! RESTARTING THE PROCESS.")
            Tree, tp_depth, tp_nodes, node_number, tp_multiSplit = initializeTree(Beta, Alpha, timepoints, time_values, v1, v2, u, theta)
            continue

        tp_x = np.random.uniform(0.0,1.0) # This value is used to determine if a node at any timepoint will have a multisplit if tp_x < u.

        for i in tp_nodes_list: # Use the tp -> nodes dict to generate the tree
            node = Tree[i]
            node_x = np.random.uniform(0.0,1.0) # node_x is checked against v1 and v2.
            print(" x value of node ",node.getID()," ",node_x," leaf value ",node.if_leaf)
            #print(node.is_dead," ",not node.is_dead)
    
            if node_x < v1 and node.if_leaf != -1: # node dies in this condition. This is not done and needs more thought.
                print(" Inside dead node condition ... ")
                Tree[i].if_leaf = -1
                Tree[i].mut = True
                perc = node.getPerc() # Dead node will have parent node percentage.
                Tree[node_number].perc = perc
                Tree[node_number].tuple = node.getTuple()
                #break
            elif node.if_leaf == 1 and (node_x >= v1 and node_x < v2): # this will have one daughter node with no new mutations
                this_id = node.getID()
                #if (this_id in node_oneSplit): # or (this_id in dead_node):
                #    continue
                #node_oneSplit.add(this_id)

                Tree[i].if_leaf = 0
                print(" Inside one daughter node condition ")
                # add one daughter node
                node_number+=1
                Tree.append(treenode(node_number))
                # set parent id
                Tree[node_number].parentID = this_id
                # set depth and edge length
                depth, edge_length = get_depth_edgeLen(node.depth_, node.getEdgeLength(), time_values)
                Tree[node_number].depth_ = depth
                Tree[node_number].edge_length = edge_length
                tp_depth = depth
                print(" Updated time ",tp_depth)

                if tp_depth in tp_nodes: # Update the tp_nodes dict with new timepoints and nodes.
                    temp_node_list = tp_nodes.get(tp_depth)
                    temp_node_list.append(node_number)
                    tp_nodes[tp_depth] = temp_node_list
                else:
                    temp_node_list = []
                    temp_node_list.append(node_number)
                    tp_nodes[tp_depth] = temp_node_list

                # determine the percentage from the Beta splitting and parents' percentage.
                perc = node.getPerc()
                Tree[node_number].perc = perc # The daughter nodes will have parents' node percentage. 
                Tree[node_number].if_leaf = 1

                # The new intervals are assigned here
                Tree[node_number].tuple = node.getTuple() # Since there is just one daughter node the parents interval gets carried forward. 
                Tree[node_number].mut = False # One daughter node will not have any mutation.
                #break
            elif node.if_leaf == 1 and node_x >= v2: # This condition will split a node
                #print("===============================")
                #print(" Time ",tp_depth," u_value ",tp_x)
                
                this_id = node.getID()
                Tree[i].if_leaf = 0
                # add two more leaves
                node_number+=2
                Tree.append(treenode(node_number-1))
                Tree.append(treenode(node_number))

                # set parent id
                Tree[node_number].parentID = this_id
                Tree[node_number-1].parentID = this_id

                # determine the percentage from the Beta splitting and parents' percentage
                perc = node.getPerc()
                node_Bi = np.random.beta(float(Alpha+1),float(Beta+1),1)
                Tree[node_number - 1].perc = perc * float(node_Bi[0])
                Tree[node_number].perc = perc * (1-float(node_Bi[0]))
                Tree[node_number - 1].if_leaf = 1
                Tree[node_number].if_leaf = 1

                # The new intervals are assigned here
                a,b = node.getTuple()
                middle = float(node_Bi[0])*float((float(b)-float(a)))+float(a)
                Tree[node_number-1].tuple=[a,middle]
                Tree[node_number].tuple=[middle,b]
                #print(" Tuple values of parents and children ",a,b,middle)

                # Check the theta value and decide on mutations
                if node_x > theta:
                    Tree[node_number-1].mut = False

                if tp_x < u and tp_depth not in tp_multiSplit and node.depth_.is_integer(): # Multisplitted nodes are not allowed to multisplit further.
                    print(" One node selected for multisplit in time ",tp_depth)
                    tp_multiSplit.add(tp_depth)

                    # set Edge length of children
                    Tree[node_number].edge_length, Tree[node_number-1].edge_length = get_multiSplit_edgeLength(time_values[int(node.depth_)-1],time_values[int(node.depth_)])
                    print(" Time values ",time_values[int(node.depth_)-1]," ",time_values[int(node.depth_)])

                    # set depth
                    depth = node.depth_ + 0.5
                    Tree[node_number].depth_ = depth
                    Tree[node_number-1].depth_ = depth
                    tp_depth = depth

                    #print("============================")
                else:
                    print(" Normally split a node between timepoints ")
                    # set depth and edge length
                    depth, edge_length = get_depth_edgeLen(node.depth_, node.getEdgeLength(), time_values)
                    Tree[node_number].depth_ = depth
                    Tree[node_number-1].depth_ = depth
                    tp_depth = depth
                    Tree[node_number].edge_length = edge_length
                    Tree[node_number-1].edge_length = edge_length

                if tp_depth in tp_nodes: # Update the tp_nodes dict with new timepoints and nodes.
                    temp_node_list = tp_nodes.get(tp_depth)
                    temp_node_list.append(node_number)
                    temp_node_list.append(node_number-1)
                    tp_nodes[tp_depth] = temp_node_list
                else:
                    temp_node_list = []
                    temp_node_list.append(node_number)
                    temp_node_list.append(node_number-1)
                    tp_nodes[tp_depth] = temp_node_list
                    #break
            else:
                continue
            #elif node.if_leaf and 
                
    return Tree
    
if len(sys.argv) <= 1:
    print("""
    This generates a tree according to Beta splitting model. 
    Usage: python main.py -T [num_of_timepoints] -t_v [Time in years for each timepoint] -v1 [Helps deciding if a node dies] -v2 [Helps deciding how node splits] -u [Threshold for multiple splits] -B [beta] -o [out-file]
        -T (--timepoints)    Number of timepoints. [3]
        -t_v (--time_value)  Time in years for each timepoint. [15, 7, 0]
        -v1 (--v1)  Helps deciding if the node dies. [0.05]
        -v2 (--v2)  Helps deciding how the node splits. [0.2]
        -u (--u)    Threshold for multiple splits. [0.1]
        -B (--Beta)         Beta value in Beta splitting model. 0.5 is evenly distributed. [0.2]
        -theta (--Theta)    Theta value decides if a child node will have new mutations. [0.3]
        -o (--out-file)     Output file. [tree.csv]
    """)
    sys.exit(0)

parser = argparse.ArgumentParser(description='This outputs the longitudinal tree.')
parser.add_argument('-T', '--timepoints', default=3)
parser.add_argument('-t_v', '--time_value', nargs="*", default=[15,7,0])
parser.add_argument('-v1','--v1', default=0.05)
parser.add_argument('-v2','--v2', default=0.2)
parser.add_argument('-u','--u', default=0.1)
parser.add_argument('-B','--Beta', default=0.2)
parser.add_argument('-theta','--Theta', default=0.3)
parser.add_argument('-o', '--out-file', default="tree.csv")

args = parser.parse_args()
timepoints = args.timepoints
time_values = args.time_value
v1 = float(args.v1)
v2 = float(args.v2)
u = float(args.u)
Beta = float(args.Beta)
theta = float(args.Theta)
out_f = args.out_file

T = gen_tree(Beta, 0.5, timepoints, time_values, v1, v2, u, theta)
tree_dict = get_TreeDict(T)
print(tree_dict)
T1 = update_tpNodes_interval(tree_dict,T)
save_tree(T1, out_f)
