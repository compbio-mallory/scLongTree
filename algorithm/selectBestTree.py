# Author: Rituparna Khan

import argparse
import os
import json
import copy
import pandas as pd
import numpy as np
from selectOptimalTree import cellTimepoints, readDMatrix, timepoint_missingRate, readFPFNvalues, intialClusterResults, selectOptimalTree, plotTree, saveTree

# Calculate the threshold of each timepoint given the FP and FN for each timepoint from the BnpC runs
# Input is an array of FP OR FN rates for each timepoint from all BnpC runs
def calc_threshold(all_bnpc_errors):
    tp_error = {}
    for i in range(len(all_bnpc_errors)):
        maxVal = max(all_bnpc_errors[i])
        threshold = maxVal / 0.8
        tp_error[i] = threshold
    return tp_error

# Get the FP and FN rates for each timepoint from all BnpC runs and calculate the threshold based on the formula provided.
# Input is no. of bnpc runs, no. of timepoints, clustering results location
def getFPFN_threshold(m, t, cLoc):
    all_bnpc_FP = []
    all_bnpc_FN = []
    epsilon = float(1e-6)
    # create empty lists for each timepoint to save the FP FN values.
    for j in range(len(t)):
        all_bnpc_FP.append([])
        all_bnpc_FN.append([])

    for i in range(1,m+1):
        for j in range(len(t)):
            # Read the alpha and beta values here for each timepoint.
            # Then assign them to all_bnpc_FP[j], all_bnpc_FN[j]
            lfile = open(cLoc+"/m"+str(i)+"/"+t[j]+"/errors.txt",'r')
            lf_line = lfile.readline().rstrip('\n')
            lf_line = lfile.readline().rstrip('\n')
            lf_arr = lf_line.split('\t')
            fn_rate = float(lf_arr[3]) + epsilon
            fp_rate = float(lf_arr[5]) + epsilon
            all_bnpc_FP[j].append(fp_rate)
            all_bnpc_FN[j].append(fn_rate)
    print("All bnpc FP ",all_bnpc_FP)
    print("All bnpc FN ",all_bnpc_FN)
    tp_FP = calc_threshold(all_bnpc_FP)
    tp_FN = calc_threshold(all_bnpc_FN)
    print(tp_FP," ",tp_FN)
    return tp_FP, tp_FN

''' Save dictionary to JSON file. '''
def saveDictFile(fileObj, fileOp):
    with open(fileOp,'w') as f:
        f.write(json.dumps(fileObj))

''' Get BnpC's genotype by combining all cells genotypes from the timepoint runs. '''
def getBnpCGenotype(tp_genotypes):
    cell_genotype = []
    for genFile in tp_genotypes:
        t_df = pd.read_csv(genFile, sep='\t')
        transposed_df = t_df.T
        #print(transposed_df)
        temp_genotype = transposed_df.values.tolist()[1:]
        cell_genotype.extend(temp_genotype)
    print(len(cell_genotype)," ",len(cell_genotype[0]))
    return cell_genotype

''' Find the difference in error rates between intial BnpC errors and final timepoint errors. '''
def getErrorRatesDiff(bnpc_error, final_error):
    bnpc_error_list = list(bnpc_error.values())
    final_error_list = list(final_error.values())
    errors_diff = {}
    for i in range(len(bnpc_error)):
        diff = abs(bnpc_error_list[i] - final_error_list[i])
        errors_diff[i] = diff
    print("Diff in error ",errors_diff)
    return errors_diff 

''' Get the Standard deviation of the error rates. And eliminate runs having error rates > SD. '''
def selectRuns(diffFPs, diffFNs, noOfRuns, noOfTimepoints):
    # Now diffFPs and diffFNs are dictionaries with error rates in each timepoint.
    # We need to know the no. of BnpC runs along with no. of timepoints.
    # From each BnpC run get the FPs for each timepoints and then get the timepoint SD.
    # Compare with each timepoint SD and then select the best run.
    all_FPs_threshold = []
    all_FNs_threshold = []
    print("Diff FPs ",diffFPs)
    for t in range(len(noOfTimepoints)):
        tp_FP = []
        tp_FN = []
        for m in range(noOfRuns):
            tp_FP.append(diffFPs[m][t])
            tp_FN.append(diffFNs[m][t])
        FP_sd = np.median(tp_FP) + 1 * (np.std(tp_FP))
        FN_sd = np.median(tp_FN) + 1 * (np.std(tp_FN))
        all_FPs_threshold.append(FP_sd)
        all_FNs_threshold.append(FN_sd)
    print(" All timepoint FPs ",all_FPs_threshold)
    print(" All timepoint FNs ",all_FNs_threshold)

    not_selected_runs = set()
    all_runs = set()
    for t in range(len(noOfTimepoints)):
        for m in range(noOfRuns):
            all_runs.add(m+1)
            if diffFPs[m][t] > all_FPs_threshold[t] or diffFNs[m][t] > all_FNs_threshold[t]:
                not_selected_runs.add(m+1)
    print("Not selected runs ",not_selected_runs)
    return list(all_runs - not_selected_runs)

''' From all the selected runs select the BnpC run with highest prob '''
def selectBnpCBestRun(bnpc_runs, bnpc_prob):
    prob_list = []
    if bnpc_runs == []: # If no run selected with median+SD then select run with highest prob
        return max(bnpc_prob, key= lambda x: bnpc_prob[x])
    for i in bnpc_runs:
        prob_list.append(bnpc_prob.get(i))
    maxProb_index = np.argmax(prob_list)
    bnpcRun_maxProb = bnpc_runs[maxProb_index]
    return bnpcRun_maxProb

''' Return the Tree with the highest prob from all the BnpC runs. '''
# Input: m = no. of Bnpc runs, t = no. of timepoints, cLoc = clustering results location
def getTreeWithHighestProb(m, t, cLoc, tp_MR, tpCells, D_matrix, k, tp_fp, tp_fn, plotOp, sample, max_parallel_extra, preserve_mutations=False):
    Tree_prob = {} # Save the trees for each BnpC run
    Tree_backMut = {} # Save back mutations for each BnpC run
    Tree_beforeCorrection = {} # Save the Trees for each run before parallel and back mutation corrections
    # Save the estimated error rates by the tree with the highest prob
    bnpc_FP = {}
    bnpc_FN = {}
    # Save the BnpC runs to get the cell genotype later 
    bnpc_run_prob = {}
    corrected_bnpc_cells = {}
    corrected_bnpc_genotype = {}
    FP_diff = [] # difference between final FP rates with initial BnpC estimated rates
    FN_diff = [] # difference between final FN rates with initial BnpC estimated rates
    for i in range(1,m+1):
    #for i in range(3,4):
        #allClusters = cLoc+"/m"+str(i)+"/all/assignment.txt" # assignment.txt from BnpC run across all timepoints
        print("BnpC run ",i)
        tp_files = [] # save timepoint assignment.txt files
        tp_errors = [] # save timepoint FP FN rates
        tp_genotypes = [] # BnpC genotypes from each timepoint run
        for j in t:
            tp_files.append(cLoc+"/m"+str(i)+"/"+j+"/assignment.txt")
            tp_errors.append(cLoc+"/m"+str(i)+"/"+j+"/errors.txt")
            tp_genotypes.append(cLoc+"/m"+str(i)+"/"+j+"/genotypes_posterior_mean.tsv")
        print("Timepoint file ",tp_files)
        tp_alpha, tp_beta = readFPFNvalues(tp_errors)
        bnpc_cells_genotype = getBnpCGenotype(tp_genotypes)
        print("BnpC cells genotype ",bnpc_cells_genotype)
        print("Max TP alpha ",tp_fp)
        print("Max TP beta ",tp_fn)
        #sorted_cluster_prob, tp_reassignedCells, tp_updatedCG = intialClusterResults(tp_files, allClusters, tpCells, D_matrix, tp_alpha, tp_beta, tp_MR)
        sorted_cluster_prob, tp_reassignedCells, tp_updatedCG = intialClusterResults(tp_files, "", tpCells, bnpc_cells_genotype, D_matrix, tp_alpha, tp_beta, tp_MR)
        corrected_bnpc_cells[i] = copy.deepcopy(tp_reassignedCells)
        corrected_bnpc_genotype[i] = copy.deepcopy(tp_updatedCG)
        print("TP reassigned cells ",tp_reassignedCells)
        print("TP updated CG ",tp_updatedCG)
        #Tree, back_mut, tree_prob, final_tp_FP, final_tp_FN, treeBeforeCorrection = selectOptimalTree(D_matrix, bnpc_cells_genotype, tp_alpha, tp_beta, tp_MR, sorted_cluster_prob, tp_reassignedCells, tp_updatedCG, tpCells, k, tp_fp, tp_fn, plotOp, str(i), sample)
        Tree, back_mut, tree_prob, final_tp_FP, final_tp_FN, treeBeforeCorrection = selectOptimalTree(
                    D_matrix, bnpc_cells_genotype, tp_alpha, tp_beta, tp_MR,
                    sorted_cluster_prob, tp_reassignedCells, tp_updatedCG, tpCells,
                    k, tp_fp, tp_fn,
                    int(max_parallel_extra),
                    plotOp, str(i), sample,
                    preserve_mutations
                )

        bnpc_FP[i] = tp_alpha
        bnpc_FN[i] = tp_beta
        bnpc_run_prob[i] = tree_prob.copy()
        FP_diff.append(getErrorRatesDiff(tp_alpha, final_tp_FP))
        FN_diff.append(getErrorRatesDiff(tp_beta, final_tp_FN))
        #corrected_bnpc_cells[tree_prob.copy()] = corrected_cells
        #corrected_bnpc_genotype[tree_prob.copy()] = corrected_genotype
        Tree_prob[i] = Tree.copy()
        Tree_backMut[i] = back_mut.copy()
        Tree_beforeCorrection[i] = treeBeforeCorrection.copy()
    bnpc_runs = selectRuns(FP_diff, FN_diff, m, t)
    #max_prob = max(Tree_prob)
    #print(list(Tree_prob.keys()))
    #print("Max prob ",max_prob)
    print("Selected BnpC runs ",bnpc_runs," BnpC runs prob ",bnpc_run_prob)
    bnpcRun_maxProb = selectBnpCBestRun(bnpc_runs, bnpc_run_prob)
    print("BnpC run with max prob:",bnpcRun_maxProb)
    finalTree = Tree_prob[bnpcRun_maxProb]
    finalBackMut = Tree_backMut[bnpcRun_maxProb]
    finalTreeBeforeCorr = Tree_beforeCorrection[bnpcRun_maxProb]
    clonalCells = corrected_bnpc_cells[bnpcRun_maxProb]
    clonalGenotype = corrected_bnpc_genotype[bnpcRun_maxProb]
    saveDictFile(clonalCells, "TreeCSVs/"+plotOp+".cells.json")
    saveDictFile(clonalGenotype, "TreeCSVs/"+plotOp+".clone_genotype.json")
    print("Max FP rates ",bnpc_FP[bnpcRun_maxProb]," Max FN rates ",bnpc_FN[bnpcRun_maxProb])
    #plotTree(finalTree, {}, finalBackMut,"plots/"+plotOp+"/smallSA501.pdf", "", sample)
    # plotTree will save the .pdf file ans saveTree will save the Tree in .csv 
    plotTree(finalTree, {}, finalBackMut,"plots/"+plotOp+".pdf", "", sample)
    saveTree(finalTree, "TreeCSVs/"+plotOp+".csv")
    saveTree(finalTreeBeforeCorr, "TreeCSVs/"+plotOp+"_beforeCorr.csv")

parser = argparse.ArgumentParser()
#parser.add_argument("-tp", "--tp",dest = "tp", help="Timepoint cells genotype.")
#parser.add_argument("-cg", "--cg",dest = "cg", help="Consensus genotype from BnpC.")
# Get the timepoint clusters input as an array because we won't know how many timepoints we will have.
# Get the cells assigned in each timepoint instead of using index.
parser.add_argument("-m", "--m", help="No. of BnpC runs.")
parser.add_argument("-t", "--t", nargs="*", help="No. of time points.")
parser.add_argument("-loc", "--loc", help="Location to save the clustering results.")
parser.add_argument("-cells","--cells", help="Cells assigned at each timepoint.")
#parser.add_argument("-cluster", "--cluster",dest="cluster", help="Assignment.txt or file having cluster information")
#parser.add_argument("-f", "--f",dest="f", help="File having FP FN rates estimated from clusters")
parser.add_argument("-D","--D", help="D matrix.")
parser.add_argument("-k", "--k", help="No. of losses allowed")
#parser.add_argument("-FP", "--FP", help="FP error rate expected", default=0.03)
#parser.add_argument("-FN", "--FN", help="FN error rate expected", default=0.35)
parser.add_argument("-op", "--op", help="Path to save the resulting Tree.")
parser.add_argument("-sample", "--sample", help="Sample name to plot the graphs.", default="default")
parser.add_argument("-p", "--parallel", help="Max extra parallel mutations per locus (P). Total allowed occurrences = 1 + P.", default=0)
parser.add_argument(
    "--preserve-mutations",
    action="store_true",
    help="Do not prune a clone if pruning it would completely remove one or more mutations from the current clone genotypes."
)

args = parser.parse_args()

tpCells = cellTimepoints(args.cells)
D_matrix = readDMatrix(args.D)
tp_MR = timepoint_missingRate(tpCells, D_matrix)
tp_FP, tp_FN = getFPFN_threshold(int(args.m), args.t, args.loc)
# Clear the plots saving directory before saving the plots
#if os.path.exists("plots/"+args.op):
#    items = os.listdir("plots/"+args.op)
#    for item in items:
#        item_path = os.path.join("plots/"+args.op, item)
#        os.remove(item_path)
#        if item.startswith('.'):
#            os.remove(item_path)
#else:
#    os.makedirs("plots/"+args.op) # Make the dir to save the plots

getTreeWithHighestProb(int(args.m), args.t, args.loc, tp_MR, tpCells, D_matrix,
                       int(args.k), tp_FP, tp_FN, args.op, args.sample, int(args.parallel),
                       args.preserve_mutations)
