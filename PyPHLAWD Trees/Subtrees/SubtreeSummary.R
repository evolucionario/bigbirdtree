


setwd("~/Google Drive/BigBirdTree Calibration/Joseph BBtree/Subtrees/")

clades <- list.files()

clades

# Exclude the Readme fole
clades <- clades[clades != "README"]
clades <- clades[clades != "summary.txt"]
clades <- clades[clades != "Redundant"]
clades <- clades[clades != "SubtreeSummary.R"]

library(ape)



ntaxa <- numeric()

for(i in 1:length(clades)) {
	
	tree <- read.tree(file=clades[i])
	ntaxa[i] <- Ntip(tree)
}

ntaxa

sum(ntaxa)

# 5556

