# Author: Rituparna Khan
# Save the pairwise SNVs for scLongTree inferred tree using this script.

import argparse

''' Get ony the new acquired mutations for each node. '''
def getNewMutations(pID_cID, node_mut):
    new_mutations = {}
    for parent, children in pID_cID.items():

        if parent in new_mutations:
            p_mut = new_mutations[parent]
            old_p_mut = node_mut[parent]
            if len(set(old_p_mut) - set(p_mut)) > 0: # includes all mutations of the parent node
                p_mut = old_p_mut
            #if len(old_p_mut) != len(p_mut):
            #    p_mut.extend(old_p_mut)
        else:
            if parent == '0':
                p_mut = []
            else:
                if ',' in parent:
                    parent_list = parent.split(',')
                    p_mut = []
                    for p in parent_list:
                        p_mut.extend(node_mut[p])
                else:
                    p_mut = node_mut[parent]

        for child in children:
            c_mut = node_mut[child]
            m1 = set(c_mut) - set(p_mut)
            #print("--------------------")
            #print(" parent mut ",p_mut)
            #print("parent ",parent," child ",child," m1 ",m1)
            if m1 == set() and parent in new_mutations:
                c_mut = new_mutations[parent]
            elif len(m1) > 1:
                c_mut = list(m1)
            #print(child," ",c_mut)
            #print("--------------------")
            new_mutations[child] = c_mut
    print(" Updated mutations ",new_mutations)
    return new_mutations

''' Get the parent and child nodes along with mutations of the nodes. '''
def readTreeFile(treeFile):
    pID_cID = {}
    node_mut = {} # Include only new mutations here
    tf = open(treeFile,'r')
    tf_line = tf.readline().rstrip()
    firstLine = True

    while(tf_line != ""):
        if firstLine:
            firstLine = False
            tf_line = tf.readline().rstrip()
            continue

        tline = tf_line.split('\t')
        pID = tline[1]
        if pID == '-1' or pID == '-2': # Skip the line having root node or removed nodes
            tf_line = tf.readline().rstrip()
            continue

        mutation = tline[4].split(',')
        #new_mut_no = tline[5]
        cID = tline[0]
        
        node_mut[cID] = mutation
        if pID in pID_cID: # Get the parent child relationship
            tchild = pID_cID[pID]
            tchild.append(cID)
            pID_cID[pID] = tchild
        else:
            tchild = []
            tchild.append(cID)
            pID_cID[pID] = tchild

        #if pID in node_mut: # Get the mutations. Maybe need to consider the last variable to get the correct # of mutations.
        #    p_mut = node_mut[pID]
        #    if new_mut_no == '0':
        #        node_mut[cID] = p_mut
        #    else:
        #        c_mut = set(mutation) - p_mut
        #        if c_mut == set():
        #            c_mut = set(mutation)
        #        node_mut[cID] = c_mut
        #else:
        #    node_mut[cID] = set(mutation)
        
        #print(pID," ",cID," ",node_mut)
        tf_line = tf.readline().rstrip()
    print(" Parent child relation ",pID_cID)
    print(" Mutations ",node_mut)
    return pID_cID, node_mut
  
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
                if checkPairedSNVs(SNV_pairs,mut[i],mut[j]) or mut[i] == mut[j]:
                    continue
                SNV_pairs.append(mut[i]+'_'+mut[j])
    #print("SNV pair on same branch ",SNV_pairs)
    return set(SNV_pairs)

''' SNVs appearing on ancestral branch. '''
def findSNVpairs_OnAncestralBranches(pID_cID,node_mut,sameSNVs):
    SNV_pairs = []
    for pid,cid in pID_cID.items():
        if pid == '0':
            continue
        parent_nodeId = pid
        if ',' in parent_nodeId:
            parent_list = parent_nodeId.split(',')
            p_mut = []
            for p in parent_list:
                p_mut.extend(node_mut[p])
        else:
            p_mut = node_mut[parent_nodeId]
        #print(" Parent mut ",p_mut)
        for child in cid:
            child_mut = node_mut[child]
            #print(" Child mut ",child_mut)
            for pmut in p_mut:
                for cmut in child_mut:
                    if checkPairedSNVs(SNV_pairs,pmut,cmut) or pmut == cmut or pmut+'_'+cmut in sameSNVs or cmut+'_'+pmut in sameSNVs:
                        continue
                    SNV_pairs.append(pmut+'_'+cmut)
    #print(" Ancestral SNVs ",set(SNV_pairs))
    return set(SNV_pairs)

''' parallel branches are the ones that originate from a single node. Compare the mutations of those branches here. '''
def findSNVpairs_OnParallelBranches(pID_cID,node_mut):
    # Children from a parent are the parallel edges.
    parallel_nodes = []
    SNV_pairs = []
    for pid,cid in pID_cID.items():
        if len(cid) > 1:
            parallel_nodes.append(cid)
    print(parallel_nodes)

    for pnodes in parallel_nodes:
        for i in range(len(pnodes)):
            for j in range(i+1,len(pnodes)):
                for n1m in node_mut[pnodes[i]]:
                    for n2m in node_mut[pnodes[j]]:
                        if checkPairedSNVs(SNV_pairs,n1m,n2m) or n1m == n2m:
                            continue
                        SNV_pairs.append(n1m+'_'+n2m)

    #print(" Parallel SNVs ",(SNV_pairs))
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
                    if checkPairedSNVs(SNV_pairs,m1,m2) or m1 == m2 or (m1+'_'+m2 in same_SNVpairs) or (m1+'_'+m2 in ancestral_SNVpairs) or (m2+'_'+m1 in same_SNVpairs) or (m2+'_'+m1 in ancestral_SNVpairs):
                        continue
                    SNV_pairs.append(m1+'_'+m2)
    return set(SNV_pairs)

''' save the SNV pairs to a file. '''
def save_files(SNVset,fileName):
    with open(fileName,'w') as f:
        f.write(str(SNVset))

parser = argparse.ArgumentParser()
parser.add_argument("-tree", "--tree",dest = "tree", help="Inferred longitudinal tree.")
parser.add_argument("-output", "--output",dest ="output", help="Directory to save the SNV pairs")
args = parser.parse_args()
pID_cID, node_mut = readTreeFile(args.tree)
print(" Node mut ",node_mut)
new_mutations = getNewMutations(pID_cID, node_mut)
print(" New mutations ",new_mutations)
SNV_pairs_onSameBranch = findSNVpairs_OnSameBranches(new_mutations)
SNV_pairs_OnAncestralBranch = findSNVpairs_OnAncestralBranches(pID_cID, node_mut, [])
#SNV_pairs_OnParallelBranches = findSNVpairs_OnParallelBranches(pID_cID,new_mutations)
SNV_pairs_inComp = findSNVpairs_incomparable(new_mutations, SNV_pairs_onSameBranch, SNV_pairs_OnAncestralBranch)

print(" SNVs on the same branch ",len(SNV_pairs_onSameBranch))
print(" SNVs on ancestral branch ",len(SNV_pairs_OnAncestralBranch))
print(" Incomparable SNV pairs ",len(SNV_pairs_inComp))

# save the SNV pairs in a text file
save_files(SNV_pairs_onSameBranch,args.output+'sameSNVs.txt')
#save_files(SNV_pairs_OnParallelBranches,args.output+'inComSNVs.txt')
save_files(SNV_pairs_inComp,args.output+'inComSNVs.txt')
save_files(SNV_pairs_OnAncestralBranch,args.output+'ancestralSNVs.txt')

