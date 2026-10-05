# Author: Rituparna Khan

import argparse

''' Get the parent and child nodes from the tree. '''
def getParentChildNode(tree_f):
    fileName = open(tree_f,'r')
    line_f = fileName.readline().rstrip('\n')
    pID_cID = {}
    firstLine = True
    deadNodes = []
    while(line_f != ""):
        if firstLine:
            firstLine = False
            line_f = fileName.readline().rstrip('\n')
            continue
        line_arr = line_f.split('\t') # column values from the file
        timepoint = line_arr[0]
        pID = line_arr[2]
        cID = line_arr[1]
        key = timepoint+'_'+pID
        deadNode = line_arr[5]
        usc_tp = float(timepoint).is_integer()
        # This is only for unobserved subclones which do not extend any further.
        # Since these nodes don't have any cells sequenced they shouldn't be used for evaluation.
        if deadNode == '-1' and usc_tp == False:
            print(" Inside dead node condn ",usc_tp)
            deadNodes.append(cID)
            line_f = fileName.readline().rstrip('\n')
            continue
        if key in pID_cID:
            temp = pID_cID[key]
            temp.append(cID)
            pID_cID[key] = temp
        else:
            temp = []
            temp.append(cID)
            pID_cID[key] = temp
        line_f = fileName.readline().rstrip('\n')
    print(" Parent child nodes ",pID_cID)
    return pID_cID, deadNodes

''' Get unobserved subclones that do not have children. '''
def getDeadUnobservedSubclones(parent_child, dead_usc_nodes):
    # {'1_-1': ['0'], '1.5_0': ['1', '2'], '2.0_1': ['3'], '2.0_2': ['4'], '2.5_3': ['5', '6'], '3.0_4': ['7']}
    parent_list = set()
    usc_nodes = []
    # If any unobserved child nodes don't appear as parents then they should be removed
    for p, c in parent_child.items():
        timepoint = p.split('_')[0]
        parent = p.split('_')[1]
        parent_list.add(parent)
        if float(timepoint).is_integer() == False:
            usc_nodes.extend(c)

    for un in usc_nodes:
        if un not in parent_list:
            dead_usc_nodes.append(un)
    print(" Dead unobserved subclones ",dead_usc_nodes)
    return dead_usc_nodes

''' Get the mutations on the edges of the tree. '''
def getMutationsOnEdges(mutFile, pID_cID, deadNodes):
    fileName = open(mutFile,'r')
    line_f = fileName.readline().rstrip('\n')
    node_mut = {} # contains only new mutations 
    node_mut_forAncestral = {} # contains node mutations as it appears
    while(line_f != ""):
        line_arr = line_f.split('\t') # column values from the file
        nodeId = line_arr[0]
        if nodeId in deadNodes:
            line_f = fileName.readline().rstrip('\n')
            continue
        mut = line_arr[1].split(';')
        node_mut[nodeId] = mut
        line_f = fileName.readline().rstrip('\n')
    print("Node mutations read ",node_mut)

    mut_edge = {}
    for pID,cIDs in pID_cID.items():
        pid = (pID.split('_'))[1]
        tp = (pID.split('_'))[0]
        if pid == '-1' or pid == '0':
            continue

        p_mut = node_mut[pid]
        print(pid,"Parent mut ",p_mut)
        if pid not in node_mut_forAncestral:
            node_mut_forAncestral[pid] = p_mut

        for cid in cIDs:
            if cid == '0' or cid in deadNodes: # Root node condition
                continue
            c_mut = node_mut[cid]
            if p_mut == [] and float(tp) >= 2:
                print("Inside this condn")
                node_mut[cid] = []
                node_mut_forAncestral[cid] = []
            else:
                node_mut[cid] = list(set(c_mut) - set(p_mut)) # update with only new mutations
                if node_mut[cid] == []:
                    merged_mut = []
                else:
                    if pid in node_mut_forAncestral: # get updated parent mutations
                        p_mut = node_mut_forAncestral[pid]

                    parent_mut = list(set(p_mut))
                    child_mut = list(set(c_mut))
                    merged_mut = list(set(parent_mut + child_mut))
                #print("PID ",pID," Parent mut ",parent_mut)
                #print(cid," Child mut ",child_mut)
                node_mut_forAncestral[cid] = merged_mut # update this dictionary with all parent mutations to check for ancestral relations
                print(node_mut_forAncestral)

            edge = pID+'_'+cid
            for m in c_mut:
                if m in mut_edge:
                    temp = mut_edge[m]
                    temp.append(edge)
                    mut_edge[m] = temp # Saves the edges for each mutation/SNVs
                else:
                    temp = []
                    temp.append(edge)
                    mut_edge[m] = temp
    print(" Mutation edges ")
    print(" Node edges ",node_mut)
    print(" Ancestral Node mutations ",node_mut_forAncestral)
    return mut_edge, node_mut, node_mut_forAncestral

# if pair of SNVs already checked then prevent checking
def checkPairedSNVs(SNV_pairs,mut1,mut2):
    for mut in SNV_pairs:
        mut_list = mut.split('_')
        #print(mut1," ",mut2," ",mut_list)
        if str(mut1) in mut_list and str(mut2) in mut_list:
            return True

''' List SNVpairs that appear on the same branch. '''
def findSNVpairs_OnSameBranches(node_mut):
    SNV_pairs = []
    for node, mut in node_mut.items():
        for i in range(len(mut)):
            for j in range(i+1,len(mut)):
                if mut[i] == mut[j] or checkPairedSNVs(SNV_pairs,mut[i],mut[j]):
                    continue
                SNV_pairs.append(mut[i]+'_'+mut[j])
    #print("SNV pair on same branch ",set(SNV_pairs))
    return set(SNV_pairs)

''' SNVs appearing on ancestral branch. '''
def findSNVpairs_OnAncestralBranches(pID_cID, node_mut, sameSNV_pairs, deadNodes):
    SNV_pairs = []
    for pid,cid in pID_cID.items():
        if '1_-1' in pid or '1.5_0' in pid:
            continue
        parent_nodeId = pid.split('_')[1]
        p_mut = node_mut[parent_nodeId]
        #print(" Parent mut ",p_mut)
        for child in cid:
            if child in deadNodes:
                continue
            child_mut = node_mut[child]
            #print(" Child mut ",child_mut)
            for pmut in p_mut:
                for cmut in child_mut:
                    if pmut == cmut or checkPairedSNVs(SNV_pairs,pmut,cmut) or pmut+'_'+cmut in sameSNV_pairs or cmut+'_'+pmut in sameSNV_pairs:
                        continue
                    SNV_pairs.append(pmut+'_'+cmut)
    #print(" Ancestral SNVs ",set(SNV_pairs))
    return set(SNV_pairs)

''' parallel branches are the ones that originate from a single node. Compare the mutations of those branches here. '''
def findSNVpairs_OnParallelBranches(pID_cID,node_mut):
    # Parallel nodes are the children of a parent.
    parallel_nodes = []
    SNV_pairs = []
    for pid,cid in pID_cID.items():
        if len(cid) > 1:
            parallel_nodes.append(cid)
        #timepoint = pid.split('_')[0]
        #for child in cid:
        #    if timepoint in tp_nodes:
        #        t_node = tp_nodes[timepoint]
        #        t_node.append(child)
        #        tp_nodes[timepoint] = t_node
        #    else:
        #        t_node = []
        #        t_node.append(child)
        #        tp_nodes[timepoint] = t_node
    print(parallel_nodes)

    for pnodes in parallel_nodes:
        n1 = pnodes[0]
        n2 = pnodes[1]

        n1_mut = node_mut[n1]
        n2_mut = node_mut[n2]

        for n1m in n1_mut:
            for n2m in n2_mut:
                if n1m == n2m:
                    continue
                SNV_pairs.append(n1m+'_'+n2m)
    #print(" Parallel SNVs ",SNV_pairs)
    return set(SNV_pairs)

''' Find incomparable SNV pairs. '''
def findSNVpairs_incomparable(node_mut, SNV_pairs_onSameBranch, SNV_pairs_OnAncestralBranch):
    SNV_pairs = []
    same_SNVpairs = list(SNV_pairs_onSameBranch)
    ancestral_SNVpairs = list(SNV_pairs_OnAncestralBranch)
    for node1, mut1 in node_mut.items():
        for node2, mut2 in node_mut.items():
            for m1 in mut1:
                for m2 in mut2:
                    if m1 == m2 or checkPairedSNVs(SNV_pairs,m1,m2) or (m1+'_'+m2 in same_SNVpairs) or (m1+'_'+m2 in ancestral_SNVpairs) or (m2+'_'+m1 in same_SNVpairs) or (m2+'_'+m1 in ancestral_SNVpairs):
                        continue
                    SNV_pairs.append(m1+'_'+m2)
    return set(SNV_pairs)

''' save the SNV pairs to a file. '''
def save_files(SNVset,fileName):
    with open(fileName,'w') as f:
        f.write(str(SNVset))

parser = argparse.ArgumentParser()
parser.add_argument("-tree", "--tree",dest = "tree", help="The ground truth tree.")
parser.add_argument("-mut", "--mut",dest = "mut", help="Mutations file.")
parser.add_argument("-output", "--output",dest ="output", help="Directory to save the SNV pairs")
args = parser.parse_args()

pID_cID, deadNodes = getParentChildNode(args.tree)
print("=====PARENT CHILD ==== ",pID_cID)
deadNodes = getDeadUnobservedSubclones(pID_cID, deadNodes)
mut_edge, node_mut, node_mut_forAncestral = getMutationsOnEdges(args.mut, pID_cID, deadNodes)
SNV_pairs_onSameBranch = findSNVpairs_OnSameBranches(node_mut)
SNV_pairs_OnAncestralBranch = findSNVpairs_OnAncestralBranches(pID_cID, node_mut_forAncestral, SNV_pairs_onSameBranch, deadNodes)
#SNV_pairs_OnParallelBranches = findSNVpairs_OnParallelBranches(pID_cID,node_mut)
SNV_pairs_inComparable = findSNVpairs_incomparable(node_mut, SNV_pairs_onSameBranch, SNV_pairs_OnAncestralBranch)

# save the SNV pairs in a text file
save_files(SNV_pairs_onSameBranch,args.output+'/sameSNVs.txt')
#save_files(SNV_pairs_OnParallelBranches,args.output+'/inComSNVs.txt')
save_files(SNV_pairs_inComparable,args.output+'/inComSNVs.txt')
save_files(SNV_pairs_OnAncestralBranch,args.output+'/ancestralSNVs.txt')
