# Author: Rituparna Khan
# This file can be used to annotate simulated data for LACE's input

import argparse

def get_header(mutations):
    header_arr = []
    header_arr.append('""')
    for i in range(len(mutations)):
        header_arr.append('"'+"SNV"+str(i)+'"')
    #print(header_arr)
    return header_arr

def process_input(input_file, output_file):
    file = open(input_file,"r")
    line = file.readline().rstrip('\n')
    noOfMutations = line.split('\t')
    #print(" No. of mutations ",len(noOfMutations))
    modified_lines = []
    header_arr = get_header(noOfMutations)
    header_arr_str = "\t".join(str(x) for x in header_arr)
    #print(" Header ",header_arr_str)
    count=1 # Uncomment this line to process real data.
    while(line != ""):
        #if count == 1:
        #    line = file.readline().rstrip('\n')
        #    count=count+1
        #    continue
        #line = line.replace("\t",",")
        #print(" Line ",line)
        modified_lines.append(line)
        line = file.readline().rstrip('\n')
    with open(output_file,"w+") as f:
        f.write(header_arr_str)
        f.write("\n")
        #print(" No of cells ",len(modified_lines))
        for l in range(1,len(modified_lines)+1):
            f.write('"'+"cell"+str(l)+'"'+'\t'+modified_lines[l-1])
            #f.write(modified_lines[l-1])
            f.write("\n")

def get_cellsInTimepoint(cells_file):
    cfile = open(cells_file,'r')
    cline = cfile.readline().rstrip('\n')
    tp_cells = {}
    while(cline != ""):
        line_a = cline.split('\t')
        timepoint = (line_a[0].split('_'))[0]
        cells = line_a[1].split(';')
        if timepoint in tp_cells:
            temp = tp_cells[timepoint]
            temp.extend(cells)
            tp_cells[timepoint] = temp
        else:
            tp_cells[timepoint] = cells
        cline = cfile.readline().rstrip('\n')
    
    cells_str = ""
    for tp, cells in tp_cells.items():
        noOfCells = len(cells)
        cells_str+= str(noOfCells) + ' '

    print(cells_str)

parser = argparse.ArgumentParser()
parser.add_argument("-input", "--input",dest ="input", help="Input data matrix")
parser.add_argument("-output", "--output",dest ="output", help="Output data matrix")
parser.add_argument("-cells","--cells",dest="cells", help="The file containing all SNV cells")
args = parser.parse_args()

process_input(args.input, args.output)
get_cellsInTimepoint(args.cells)
