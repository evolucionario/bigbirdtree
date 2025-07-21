######################################
### Assembly of the Big Bird Tree  ###
######################################



library(ape)
library(phangorn)

#setwd('~/Google Drive/BigBirdTree Calibration/BigBirdTree Assembly/') # Home

setwd('~/GoogleDrive/BigBirdTree Calibration/BigBirdTree Assembly')


### Read BBtree (if needed) ###
# BBtree <- read.tree(file="BBtree.tre")

Sptable <- read.csv("SpeciesTable.csv")

summary(Sptable)



### Included species per family ###

Included <- table(Sptable$Family2)

length(Included)



### Calculate total species per family ###
# Using Clement's checklist

Clement <- read.csv("NEW_Clements-Checklist-v2022-October-2022.csv")

# Kee only species-level taxa

Clement <- Clement[Clement$category=="species",]

nrow(Clement)

# Delete common name from family 
ClementFamTemp <- strsplit(Clement$family, split=" ")

ClementFam <- character()
for(i in 1:length(ClementFamTemp)) {
ClementFam[i] <- ClementFamTemp[[i]][1]
}

FamilyRichness <- table(ClementFam)

length(FamilyRichness)


# Expand richness values to Sptable

FamilyRichness[Sptable$Family2]

Included[Sptable$Family2]


Sptable2 <- cbind(Sptable, FamilyTotal=FamilyRichness[Sptable$Family2], FamilyIncluded= Included[Sptable$Family2])

write.csv(Sptable2, file="SpeciesTable3.csv")