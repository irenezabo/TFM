#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Needed arguments:
# 1- inputfile: whole path to iupred results
# 2- ncolID: number of the column with the target gene or protein id
# 3- ncolDis: number of the column with disorder estimate to use 
# 4- outputfile: whole path to output file 

#Get the arguments
if (length(args)!=4) {
  stop("At least 4 arguments needes: input file, column with gene/prot id, column with dis information (num) and outputdir.\n", call.=FALSE)
} else if (length(args)==4) {
  inputfile=args[1]
  ncolID=paste0("V",args[2])
  ncolDis=paste0("V",args[3])
  outputfile=args[4]
}

# Load library
require(dplyr)

# Load file:
iupred<-read.table(inputfile, 
                   stringsAsFactors = T, header=F, na.strings="na", sep = "\t",
                   quote = "",comment.char = "", fill = TRUE)


# Colnames
iupred <- iupred %>% rename("ID" =ncolID,
                            "Dis"=ncolDis)
#print(head(iupred))


# Get table
dis_prot<-iupred %>%
  group_by(ID) %>%
  summarise(pDis=round(mean(Dis),4))

# Save table
write.table(dis_prot, outputfile, 
            append = FALSE, sep = "\t", dec = ".", 
            quote=F, row.names = F, col.names = F)


