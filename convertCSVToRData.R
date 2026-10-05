# As an input we need the following: the input D matrix, cell nos for each time point (make it dynamic so that the no of dataframes to create is decided based on this), output file name
# This is the main structure to follow. But we need to change the simulated csv files to a file with cell IDs and SNVs.
args = commandArgs(trailingOnly=TRUE)
#str(args)
op_fileName = args[length(args)] # Last argument is always a output filename.
#str(op_fileName)

data <- read.csv(args[1], sep="\t", stringsAsFactors = F,head = T, row.names = 1) # We read the D matrix here which is the first argument
data <- as.matrix(data)
str(data)

startIndex = 1
endIndex = 0
# In the following loop we dynamically read the cell numbers and create dataframes based on that. Each timepoint will have one dataframe
for (i in 3:length(args)-1){
	cell_no = args[i]
	str(i)
	if(i==2){
		str("Inside i=3")
		master_df <- data[1:cell_no,]
		endIndex = cell_no
	}else{
		#startIndex = startIndex+as.numeric(args[i-1])+1
		startIndex = as.numeric(endIndex)+1
		endIndex = as.numeric(startIndex)+as.numeric(args[i])-1
		#endIndex = endIndex+as.numeric(args[i-1])+as.numeric(args[i])
		str(startIndex)
		str(endIndex)
		master_df <- data[startIndex:endIndex,]
	}
	str(master_df)
        assign(paste0("T", i), master_df)  # dynamically (re) name df_master
    }

ls(pattern = "^T\\d+$")
final_df <- mget(ls(pattern = "^T\\d+$")) # We keep the list of dataframes at each timepoint in one final dataframe. 
save(final_df, file=op_fileName)

#do.call(rbind, mget(ls(pattern = "^df_\\d+$")))
# Get the cell nos and decide how many dataframes to split the data in based on timepoints.
#data_frame_1 <- data[1:100,]
#data_frame_2 <- data[101:600,]
#str(data_frame_1)
#str(data_frame_2)
#df_list <- list("T1" = data_frame_1,"T2" = data_frame_2) # name each dataframe as the timepoints.
#save(df_list, file=op_fileName)

