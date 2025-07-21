######################################
### Assembly of the Big Bird Tree  ###
######################################


library(ape)
library(phangorn)


setwd('BigBirdTree Assembly')


### Ancillary Functions ####

# getMRCAage returns the age of the MRCA of a pair of tips

getMRCAage <- function(phy, tip) {

	# Obtain branching times
	
	BTimes <- branching.times(phy)
	
	# Eliminate names
	names(BTimes) <- NULL
	
	# identify the MRCA node number
	
	MRCA <- getMRCA(phy, tip)

	# substract Ntip to obtain the index of the node in the branching times result

	MRCAnode <- MRCA - Ntip(phy)
	
	return(BTimes[MRCAnode])	
	
}



### Load RAG Backnone Tree ###

BackTre <- read.tree("ChronogramNeornithes.tre")

plot(BackTre, edge.width=0.5, cex=0.1)

axisPhylo(cex.axis=0.7)


### Load phyPhlawd big tree

PhlawdTree <- read.tree("PyPHLAWD Trees/RAxML_bestTree.Aves_bait_13.tre") # Home






#####################
### Palaeognathae ###
#####################

# Load PyPHLAWD Palaeognathae tree 

Palaeo <- read.tree("PyPHLAWD Trees/Subtrees/Palaeognathae_iqtree_partitioned-ML_aLRT.treefile")

# Root tree
Palaeo <- root(Palaeo, outgroup="Struthio_camelus")


### APTERYGIDAE ###

plot(Palaeo, cex=0.7)

# Extract Apterygidae and outgroups, as we need to calibrate the timne of origin
Apterygidae <- extract.clade(Palaeo, node=getMRCA(Palaeo, tip=c("Apteryx_haastii", "Casuarius_bennetti")))

plot(Apterygidae)

# Drop Subspecies
Apterygidae <- drop.tip(Apterygidae, tip=c("Apteryx_australis_australis", "Apteryx_mantelli_mantelli", "Dromaius_ater"))

plot(Apterygidae)


# Calibrate with RAG tree date

CladeAge <- getMRCAage(BackTre, tip=c("Apteryx_australis", "Casuarius_casuarius"))
#7.45

CALIB <- makeChronosCalib(Apterygidae, node="root", age.min=CladeAge)

ApterygidaeC <- chronos(Apterygidae, model="clock", calibration=CALIB)

plot(ApterygidaeC); axisPhylo()


# Delete outgroup
ApterygidaeC <- drop.tip(ApterygidaeC, tip=c("Casuarius_bennetti", "Dromaius_novaehollandiae"), trim.internal = TRUE)

plot(ApterygidaeC); axisPhylo()

# Restore root branch
ApterygidaeC$root.edge <-  getMRCAage(BackTre, tip=c("Apteryx_australis", "Eudromia_elegans")) - getMRCAage(ApterygidaeC, tip=c("Apteryx_australis", "Apteryx_haastii"))


plot(ApterygidaeC, root.edge=TRUE); axisPhylo()


# Delete clade in Backbone tree

BackTre2 <- drop.tip(BackTre, tip="Apteryx_australis", trim.internal=FALSE)

# Because the single deleted tip, there is no internal node converted to tip to which to attach the subclade.

# Graft Apterygidae tree into Backbone Tree
# Identify the edge of BackTre2 (phy1) where ApterygidaeC (ph2) will be attached to
# If Apterygidae sister to Tinamidae like in the RAG tree, this is the basal branch of Tinamidae!

# First get the node
Node <- getMRCA(BackTre2, tip=c("Eudromia_elegans", "Tinamus_guttatus"))

# Then get the branch
Edge <- which(BackTre2$edge[,2] == Node)

# get the age of the MRCA
Age0 <- getMRCAage(BackTre2, tip=c("Eudromia_elegans", "Tinamus_guttatus"))


#time of origin of Apterygidae 	
AgeAptery <- max(branching.times(ApterygidaeC)) + ApterygidaeC$root.edge

# Calculate the position below the target node where the clade will be attached
Pos <- AgeAptery - Age0

BackTre3 <- bind.tree(BackTre2, ApterygidaeC, where=Node, position=Pos)

plot(BackTre3, cex=0.1, edge.width=0.3)

# OR... use the new graft.tree function:
#BackTre4 <- graft.tree(BackTre2, ApterygidaeC, edge=c("Eudromia_elegans", "Tinamus_guttatus"))

plot(BackTre3, cex=0.1, edge.width=0.3)

#write.tree(BackTre3, file="BBtree1.tre")



### TINAMIDAE ###


# Use the Almeida et al. tree:

Tinamidae <- read.tree("Subclade Trees/Tinamidae/RAxML_mitnuc_longname.tre")

Tinamidae <- extract.clade(Tinamidae, node=getMRCA(Tinamidae, tip=c("Tinamus_tao", "Nothura_maculosa")))

# Delete subspecies 
Tinamidae <- drop.tip(Tinamidae, tip="Nothura_maculosa_chacoensis")

# Correct typo
Tinamidae$tip.label[Tinamidae$tip.label=="Crypturellys_variegatus"] <- "Crypturellus_variegatus"



plot(Tinamidae)

# Calibrate with RAG tree date

CladeAge <- getMRCAage(BackTre3, tip=c("Eudromia_elegans", "Tinamus_guttatus"))
#7.45

CALIB <- makeChronosCalib(Tinamidae, node="root", age.min=CladeAge)

TinamidaeC <- chronos(Tinamidae, model="clock", calibration=CALIB)

plot(TinamidaeC); axisPhylo()


# Delete clade in Backbone tree

BackTre4 <- drop.tip(BackTre3, tip=c("Nothocercus_nigrocapillus", "Nothoprocta_perdicaria", "Crypturellus_undulatus"), trim.internal = TRUE)

BackTre4 <- drop.tip(BackTre4, tip=c("Eudromia_elegans", "Tinamus_guttatus"), trim.internal = FALSE)

# Find the new tip number
# The new tip is allways a new tip so 
Tip <- Ntip(BackTre)

# Graft Nothurinae tree into Backbone Tree
BackTre4 <- bind.tree(BackTre4, TinamidaeC, where=Tip)

plot(BackTre4, cex=0.1, edge.width=0.3)






###################
### GALLIFORMES ###
###################

### Galliformes ###

# Load Galliformes tree from Kimball et al 2021

Galliformes <- read.nexus("Subclade Trees/Galliformes/Kimball_Galliform_Treefile.nex")

# Pick the Partitioned_PSR tree, the mitochondrial tree shows shorter deep branches that may be due to saturation
 
Galliformes <- Galliformes[[6]]

Galliformes <- root(Galliformes, outgroup="Anas_platyrhynchos")

Galliformes <- drop.tip(Galliformes, tip="Anas_platyrhynchos")

Galliformes <- ladderize(Galliformes, right = FALSE)

plot(Galliformes, cex=0.3)


# Get crown age from Backbone Tree
CladeAge <- getMRCAage(BackTre, tip=c("Megapodius_freycinet", "Gallus_gallus"))

# Build the calibration table
CALIB <- makeChronosCalib(Galliformes, node="root", age.min=CladeAge)

# Estimate chronogram
GalliformesC <- chronos(Galliformes, model="clock", calibration=CALIB)

plot(GalliformesC, cex=0.5); axisPhylo()


# Clip Galliformes from Backbone tree

BackTre5 <- drop.tip(BackTre4, tip=c("Leipoa_ocellata", "Crax_blumenbachii","Penelope_obscura","Callipepla_californica","Numida_meleagris","Guttera_pucherani","Rollulus_rouloul","Bonasa_umbellus"), trim.internal=TRUE)

BackTre5 <- drop.tip(BackTre5, tip=c("Megapodius_freycinet", "Gallus_gallus"), trim.internal=FALSE)


# Add the new megapode chronogram

BackTre6 <- bind.tree(BackTre5, GalliformesC, where=Ntip(BackTre5))

plot(BackTre6, cex=0.2)





####################
### ANSERIFORMES ###
####################

# Load Anseriformes tree

Anseriformes <- read.tree("PyPHLAWD Trees/Subtrees/Anseriformes_iqtree_partitioned-ML_aLRT.treefile")

# Delete extra subspecies 

Anseriformes <- drop.tip(Anseriformes, tip=c(
"Mareca_strepera_strepera",
"Lophonetta_specularioides_specularioides",
"Lophonetta_specularioides_alticola",
"Mergus_merganser_merganser",
"Somateria_mollissima_mollissima",
"Chloephaga_picta_picta",
"Chloephaga_picta_leucoptera",
"Chloephaga_hybrida_malvinarum",
"Oxyura_jamaicensis_jamaicensis",
"Oxyura_jamaicensis_andina"))


# Delete subspecific epitheth
Anseriformes$tip.label[which(Anseriformes$tip.label == "Chloephaga_hybrida_hybrida")] <- "Chloephaga_hybrida"

# Update genus name 

Anseriformes$tip.label[which(Anseriformes$tip.label == "Chloephaga_melanoptera")] <- "Oressochen_melanopterus"
Anseriformes$tip.label[which(Anseriformes $tip.label == "Neochen_jubata")] <- "Oressochen_jubatus"



Anseriformes <- midpoint(Anseriformes)

Anseriformes <- ladderize(Anseriformes, right = FALSE)

plot(Anseriformes, cex=0.2, edge.width=0.4)


### Anhimidae ###

# Extract Anhimidae
Anhimidae <- extract.clade(Anseriformes, node=getMRCA(Anseriformes, tip=c("Anhima_cornuta", "Chauna_torquata")))

# Get crown age from Backbone Tree
CladeAge <- getMRCAage(BackTre, tip=c("Anhima_cornuta", "Chauna_torquata"))

# Build the calibration table
CALIB <- makeChronosCalib(Anhimidae, node="root", age.min=CladeAge)

# Estimate chronogram
AnhimidaeC <- chronos(Anhimidae, model="clock", calibration=CALIB)

plot(AnhimidaeC); axisPhylo()


# Clip Anhimidae from Backbone tree

BackTre12 <- drop.tip(BackTre6, tip=c("Anhima_cornuta", "Chauna_torquata"), trim.internal=FALSE)

# Add the new megapode chronogram

BackTre14 <- bind.tree(BackTre12, AnhimidaeC, where=Ntip(BackTre12))

plot(BackTre14, cex=0.2)




### Anatidae ###

# Extract Clade

Anatidae <- extract.clade(Anseriformes, node=getMRCA(Anseriformes, tip=c("Dendrocygna_bicolor", "Anas_acuta")))

plot(Anatidae, cex=0.7); axisPhylo()

# Get crown age from Backbone Tree
CladeAge <- getMRCAage(BackTre, tip=c("Dendrocygna_arborea", "Anser_albifrons"))

# Build the calibration table
CALIB <- makeChronosCalib(Anatidae, node="root", age.min=CladeAge)

# Estimate chronogram
AnatidaeC <- chronos(Anatidae, model="discrete", calibration=CALIB, control=chronos.control(dual.iter.max=100, nb.rate.cat=3))

plot(AnatidaeC, cex=0.4); axisPhylo()


# CLIPING
# In two stepts to prevent problems with lingering internal branches

BackTre15 <- drop.tip(BackTre14, tip=c("Anas_platyrhynchos","Oxyura_jamaicensis"), trim.internal=TRUE)

plot(BackTre15, cex=0.2, edge.width=0.4, show.tip.label=FALSE);axisPhylo()


BackTre16 <- drop.tip(BackTre15, tip=c("Dendrocygna_arborea", "Anser_albifrons"), trim.internal=FALSE)

plot(BackTre16, cex=0.2, edge.width=0.4, show.tip.label=FALSE)
tiplabels(cex=0.2); axisPhylo()


# Add the new Anatidae chronogram

BackTre17 <- bind.tree(BackTre16, AnatidaeC, where=Ntip(BackTre16))

plot(BackTre17, cex=0.2)

BackTre17$node.label <- NULL

# write.tree(BackTre17, file="BackTre17.tre")




#################################
###                           ###
###          NEOAVES          ###
###                           ###
#################################


#load(file="BackTre17.R")

#plot(BackTre17, show.tip.label=FALSE)


######################
### COLUMBIMORPHAE ###
######################

### Columbidae ###

# Using Lapiedra et al 2021 comprehensive Bayesian tree of Columbidae
# The Peristerinae (Columbina, etc.) is the most basal subfamily (Johnson et al. 2010, Heupink et al 2014, Reddy et al. 2017)


# Load Lapiedra et al 2021 MCC tree generated in TreeAnnotator

Columbidae <- read.nexus(file="Subclade Trees/Columbiformes/Lapiedraetal2021MCC.tre")

plot(Columbidae, cex=0.7); axisPhylo()

cat(Columbidae$tip.label)

# Delete Subspecies

# Get crown age from Backbone Tree
CladeAge <- getMRCAage(BackTre, tip=c("Columbina_squammata", "Treron_calva"))

# Rescale tree
wrong.age <- max(branching.times(Columbidae))

Columbidae$edge.length <- Columbidae$edge.length*CladeAge/wrong.age

plot(Columbidae, cex=0.4); axisPhylo()


# Clip crown clade from Backbone tree

BackTre18 <- drop.tip(BackTre17, tip="Columba_livia", trim.internal=TRUE)

BackTre18 <- drop.tip(BackTre18, tip=c("Columbina_squammata","Treron_calva"), trim.internal=FALSE)

# Add the new Columbdiae chronogram

BackTre18 <- bind.tree(BackTre18, Columbidae, where=Ntip(BackTre18))

plot(BackTre18, cex=0.2)




### Pteroclidae ###

# Extract Clade
Pteroclidae <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Pterocles_coronatus", "Syrrhaptes_paradoxus")))

plot(Pteroclidae, cex=0.7); axisPhylo()

# Get crown age from Backbone Tree
CladeAge <- getMRCAage(BackTre, tip=c("Pterocles_personatus", "Syrrhaptes_paradoxus"))

# Build the calibration table
CALIB <- makeChronosCalib(Pteroclidae, node="root", age.min=CladeAge)

# Estimate chronogram
PteroclidaeC <- chronos(Pteroclidae, model="clock", calibration=CALIB, control=chronos.control(dual.iter.max=100))

PteroclidaeC <- ladderize(PteroclidaeC, right=FALSE)
plot(PteroclidaeC, cex=0.4); axisPhylo()


# Clip crown clade from Backbone tree

BBtree19 <- drop.tip(BackTre18, tip=c("Pterocles_personatus", "Syrrhaptes_paradoxus"), trim.internal=FALSE)

# Add the new Columbdiae chronogram

BBtree19 <- bind.tree(BBtree19, PteroclidaeC, where=Ntip(BBtree19))

plot(BBtree19, cex=0.2, show.tip.label=FALSE)



############################
### PHOENICOPTERIMORPHAE ###
############################



# Extract Clade
Phoenico <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Podiceps_cristatus", "Phoenicopterus_ruber")))

plot(Phoenico, cex=0.7); axisPhylo()

# Get crown age from Backbone Tree
CladeAge <- getMRCAage(BackTre, tip=c("Podiceps_cristatus", "Phoenicopterus_ruber"))

# Build the calibration table
CALIB <- makeChronosCalib(Phoenico, node="root", age.min=CladeAge)

# Estimate chronogram
PhoenicoC <- chronos(Phoenico, model="clock", calibration=CALIB, control=chronos.control(dual.iter.max=100))

PhoenicoC <- ladderize(PhoenicoC, right=FALSE)
plot(PhoenicoC, cex=0.4); axisPhylo()


# Clip crown clade from Backbone tree
BBtree20 <- drop.tip(BBtree19, tip="Tachybaptus_ruficollis", trim.internal=TRUE)

BBtree20 <- drop.tip(BBtree20, tip=c("Podiceps_cristatus", "Phoenicopterus_ruber"), trim.internal=FALSE)

# Add the new Columbdiae chronogram

BBtree20 <- bind.tree(BBtree20, PhoenicoC, where=Ntip(BBtree20))

plot(BBtree20, cex=0.2, show.tip.label=FALSE)



####################
### OTIDIMORPHAE ###
####################


### Otididae ###


# Extract Clade
Otididae <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Ardeotis_kori", "Otis_tarda")))

plot(Otididae, cex=0.7); axisPhylo()

# Delete subspecies

Otididae <- drop.tip(Otididae, tip=c("Chlamydotis_undulata_undulata", "Chlamydotis_undulata_fuertaventurae", "Otis_tarda_tarda", "Otis_tarda_dybowskii"))


# Get crown age from Backbone Tree
CladeAge <- getMRCAage(BackTre, tip=c("Ardeotis_kori", "Otis_tarda"))

# Build the calibration table
CALIB <- makeChronosCalib(Otididae, node="root", age.min=CladeAge)

# Estimate chronogram
OtididaeC <- chronos(Otididae, model="clock", calibration=CALIB, control=chronos.control(dual.iter.max=100))

OtididaeC <- ladderize(OtididaeC, right=FALSE)
plot(OtididaeC); axisPhylo()


# Clip crown clade from Backbone tree

BBtree21 <- drop.tip(BBtree20, tip=c("Ardeotis_kori", "Otis_tarda"), trim.internal=FALSE)

# Add the new Columbdiae chronogram

BBtree21 <- bind.tree(BBtree21, OtididaeC, where=Ntip(BBtree21))

#plot(BBtree21, cex=0.2, show.tip.labels=FALSE)



### Musophagidae ###

# Extract Clade
Musophagidae <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Corythaeola_cristata", "Tauraco_erythrolophus")))

plot(Musophagidae, cex=0.9); axisPhylo()

# Delete duplicated taxa and subspecies

Musophagidae <- drop.tip(Musophagidae, tip=c("Tauraco_persa", "Gallirex_porphyreolophus" , "Tauraco_macrorhynchus" , "Tauraco_livingstonii", "Tauraco_fischeri"  , "Tauraco_schalowi" ,  "Tauraco_leucotis", "Corythaixoides_concolor" , "Corythaixoides_personatus", "Tauraco_corythaix", "Tauraco_persa_zenkeri"))

plot(Musophagidae, cex=0.9); axisPhylo()

# Adopt the species-level taxonomy of Perktas et al. (2020), elevanting some subspecies to species
# But still conserving a genus Tauraco sensu lato.

#Perktaş, U., Groth, J. G., & Barrowclough, G. F. (2020). Phylogeography, species limits, phylogeny, and classification of the turacos (Aves: Musophagidae) based on mitochondrial and nuclear DNA sequences. American Museum Novitates, 2020(3949), 1-61.

# Transform other subspecies names into species names

Musophagidae$tip.label[which(Musophagidae$tip.label == "Gallirex_porphyreolophus_porphyreolophus")] <- "Gallirex_porphyreolophus"
Musophagidae$tip.label[which(Musophagidae$tip.label == "Gallirex_porphyreolophus_chlorochlamys")] <-
"Gallirex_chlorochlamys"

Musophagidae$tip.label[which(Musophagidae$tip.label == "Gallirex_johnstoni_kivuensis")] <- "Gallirex_kivuensis"
Musophagidae$tip.label[which(Musophagidae$tip.label == "Gallirex_johnstoni_johnstoni")] <- "Gallirex_johnstoni"
Musophagidae$tip.label[which(Musophagidae$tip.label == "Tauraco_macrorhynchus_macrorhynchus")] <- "Musophaga_macrorhyncha"
Musophagidae$tip.label[which(Musophagidae$tip.label == "Tauraco_macrorhynchus_verreauxii")] <- "Musophaga_verreauxii"
Musophagidae$tip.label[which(Musophagidae$tip.label == "Tauraco_persa_persa")] <- "Tauraco_persa"
Musophagidae$tip.label[which(Musophagidae$tip.label == "Tauraco_persa_buffoni")] <- "Tauraco_buffoni"
Musophagidae$tip.label[which(Musophagidae$tip.label == "Tauraco_schuettii_emini")] <- "Tauraco_emini"
Musophagidae$tip.label[which(Musophagidae$tip.label == "Tauraco_fischeri_fischeri")] <- "Tauraco_fischeri"
Musophagidae$tip.label[which(Musophagidae$tip.label == "Tauraco_livingstonii_reichenowi")] <- "Tauraco_reichenowi"
Musophagidae$tip.label[which(Musophagidae$tip.label == "Tauraco_corythaix_corythaix")] <- "Tauraco_corythaix"
Musophagidae$tip.label[which(Musophagidae$tip.label == "Tauraco_livingstonii_livingstonii")] <- "Tauraco_livingstonii"
Musophagidae$tip.label[which(Musophagidae$tip.label == "Tauraco_schuettii_schuettii")] <- "Tauraco_schuettii"
Musophagidae$tip.label[which(Musophagidae$tip.label == "Tauraco_schalowi_chalcolophus")] <- "Tauraco_chalcolophus"
Musophagidae$tip.label[which(Musophagidae$tip.label == "Tauraco_schalowi_loitanus")] <- "Tauraco_loitanus"
Musophagidae$tip.label[which(Musophagidae$tip.label == "Tauraco_schalowi_marungensis")] <- "Tauraco_marungensis"
Musophagidae$tip.label[which(Musophagidae$tip.label == "Tauraco_schalowi_schalowi")] <- "Tauraco_schalowi"
Musophagidae$tip.label[which(Musophagidae$tip.label == "Tauraco_leucotis_leucotis")] <- "Tauraco_leucotis"
Musophagidae$tip.label[which(Musophagidae$tip.label == "Tauraco_leucotis_donaldsoni")] <- "Tauraco_donaldsoni"
Musophagidae$tip.label[which(Musophagidae$tip.label == "Corythaixoides_concolor_bechuanae")] <- "Corythaixoides_concolor"
Musophagidae$tip.label[which(Musophagidae$tip.label == "Corythaixoides_personatus_leopoldi")] <- "Corythaixoides_personatus"


plot(Musophagidae, cex=0.9); axisPhylo()


# Get crown age from Backbone Tree
CladeAge <- getMRCAage(BackTre, tip=c("Corythaeola_cristata", "Tauraco_erythrolophus"))

# Build the calibration table
CALIB <- makeChronosCalib(Musophagidae, node="root", age.min=CladeAge)

# Estimate chronogram
MusophagidaeC <- chronos(Musophagidae, model="clock", calibration=CALIB, control=chronos.control(dual.iter.max=100))

MusophagidaeC <- ladderize(MusophagidaeC, right=FALSE)
plot(MusophagidaeC); axisPhylo()


# Clip crown clade from Backbone tree

BBtree22 <- drop.tip(BBtree21, tip=c("Corythaeola_cristata", "Tauraco_erythrolophus"), trim.internal=FALSE)

# Add the new Columbdiae chronogram

BBtree22 <- bind.tree(BBtree22, MusophagidaeC, where=Ntip(BBtree22))

#plot(BBtree22, cex=0.2, show.tip.label=FALSE)


#write.tree(BBtree22, file="BBtree21.tre")




### Cuculidae ###

# Read Cuculiformes tree

Cuculidae <- read.tree("PyPHLAWD Trees/Subtrees/Cuculiformes_iqtree_partitioned-ML_aLRT.treefile")

Cuculidae <- midpoint(Cuculidae)

Cuculidae <- ladderize(Cuculidae, right=FALSE)

plot(Cuculidae, cex=0.75, no.margin=TRUE)


# Delete subspecies

Cuculidae <- drop.tip(Cuculidae, tip=c(
"Cuculus_canorus_canorus",
"Cuculus_canorus_bakeri",
"Surniculus_lugubris_barussarum",
"Cercococcyx_montanus_patulus",
"Cercococcyx_montanus",
"Piaya_minuta",
"Coccyzus_americanus_occidentalis",
"Coccyzus_americanus_americanus",
"Dasylophus_superciliosus_cagayanensis",
"Dasylophus_superciliosus_superciliosus"
))

# Delete subspecies names

Cuculidae$tip.label[which(Cuculidae$tip.label == "Cercococcyx_montanus_montanus")] <- "Cercococcyx_montanus"

plot(Cuculidae, cex=0.6)


# Get crown age from Backbone Tree
CladeAge <- getMRCAage(BackTre, tip=c("Guira_guira", "Centropus_phasianinus"))

# Build the calibration table
CALIB <- makeChronosCalib(Cuculidae, node="root", age.min=CladeAge)

# Estimate chronogram
CuculidaeC <- chronos(Cuculidae, model="discrete", calibration=CALIB, control=chronos.control(dual.iter.max=100, nb.rate.cat = 3))

plot(CuculidaeC, cex=0.6); axisPhylo()


# Clip crown clade from Backbone tree

BBtree23 <- drop.tip(BBtree22, tip=c("Coua_cristata", "Coccyzus_americanus"), trim.internal=TRUE)

BBtree23 <- drop.tip(BBtree23, tip=c("Guira_guira", "Centropus_phasianinus"), trim.internal=FALSE)

# Add the new Columbdiae chronogram

BBtree23 <- bind.tree(BBtree23, CuculidaeC, where=Ntip(BBtree23))

#plot(BBtree23, cex=0.2, show.tip.label=FALSE)
#write.tree(BBtree23, file="BBtree23.tre")
#BBtree23 <- read.tree(file="BBtree23.tre")



########################
### CAPRIMULGIFORMES ###
########################


# Read Caprimulgiformes tree

Caprimulgiformes <- read.tree("PyPHLAWD Trees/Subtrees/Caprimulgiformes_iqtree_partitioned-ML_aLRT.treefile")


Caprimulgiformes <- midpoint(Caprimulgiformes)

Caprimulgiformes <- ladderize(Caprimulgiformes, right=FALSE)

plot(Caprimulgiformes, cex=0.6)


### BATRACOSTOMIDAE ###

Batracostomdiae <- extract.clade(Caprimulgiformes, node=getMRCA(Caprimulgiformes, tip=c("Batrachostomus_septimus", "Podargus_strigoides")))

# Get crown age from Backbone Tree
CladeAge <- getMRCAage(BackTre, tip=c("Batrachostomus_septimus", "Podargus_strigoides"))

# Build the calibration table
CALIB <- makeChronosCalib(Batracostomdiae, node="root", age.min=CladeAge)

# Estimate chronogram
BatracostomdiaeC <- chronos(Batracostomdiae, model="clock", calibration=CALIB, control=chronos.control(dual.iter.max=100))

#plot(BatracostomdiaeC); axisPhylo()


# Clip crown clade from Backbone tree

BBtree24 <- drop.tip(BBtree23, tip=c("Batrachostomus_septimus", "Podargus_strigoides"), trim.internal=FALSE)

# Add the new chronogram

BBtree24 <- bind.tree(BBtree24, BatracostomdiaeC, where=Ntip(BBtree24))

#plot(BBtree24, cex=0.2, show.tip.label=FALSE)




### Caprimulgidae ###

Caprimulgidae <- extract.clade(Caprimulgiformes, node=getMRCA(Caprimulgiformes, tip=c("Eurostopodus_argus", "Hydropsalis_brasiliana")))


# Delete subspecies

cat(Caprimulgidae$tip.label)


Caprimulgidae <- drop.tip(Caprimulgidae, tip=c(
"Setopagis_parvulus_parvulus",
"Nyctiphrynus_ocellatus_ocellatus",
"Lyncornis_macrotis_macrotis",
"Antrostomus_vociferus_vociferus",
"Uropsalis_segmentata_kalinowskii"))


# Get crown age from Backbone Tree
CladeAge <- getMRCAage(BackTre, tip=c("Eurostopodus_argus", "Caprimulgus_longirostris"))

# Build the calibration table
CALIB <- makeChronosCalib(Caprimulgidae, node="root", age.min=CladeAge)

# Estimate chronogram
CaprimulgidaeC <- chronos(Caprimulgidae, model="discrete", calibration=CALIB, control=chronos.control(dual.iter.max=100, nb.rate.cat = 3))

#plot(CaprimulgidaeC); axisPhylo()

# Clip crown clade from Backbone tree
BBtree25 <- drop.tip(BBtree24, tip=c("Lyncornis_macrotis", "Gactornis_enarratus"), trim.internal=TRUE)

BBtree25 <- drop.tip(BBtree25, tip=c("Eurostopodus_argus", "Caprimulgus_longirostris"), trim.internal=FALSE)

# Add the new chronogram

BBtree25 <- bind.tree(BBtree25, CaprimulgidaeC, where=Ntip(BBtree25))

#plot(BBtree25, cex=0.2, show.tip.label=FALSE)

#write.tree(BBtree25, file="BBtree25.tre")





### Nyctibiidae ###

Nyctibiidae <- extract.clade(Caprimulgiformes, node=getMRCA(Caprimulgiformes, tip=c("Nyctibius_grandis", "Nyctibius_griseus")))

# Note that N. grandis is basal in our phyPHLAWD tree but several previous analysis indicated N. bracteatus as the basal lineage.

# Get crown age from Backbone Tree
CladeAge <- getMRCAage(BackTre, tip=c("Nyctibius_grandis", "Nyctibius_aethereus"))

# Build the calibration table
CALIB <- makeChronosCalib(Nyctibiidae, node="root", age.min=CladeAge)

# Estimate chronogram
NyctibiidaeC <- chronos(Nyctibiidae, model="clock", calibration=CALIB, control=chronos.control(dual.iter.max=100))

#plot(NyctibiidaeC); axisPhylo()


# Clip crown clade from Backbone tree

BBtree26 <- drop.tip(BBtree25, tip=c("Nyctibius_grandis", "Nyctibius_aethereus"), trim.internal=FALSE)

# Add the new chronogram

BBtree26 <- bind.tree(BBtree26, NyctibiidaeC, where=Ntip(BBtree26))

#plot(BBtree26, cex=0.2, show.tip.label=FALSE)

#write.tree(BBtree26, file="BBtree26.tre")



### Aegothelidae ###

Aegothelidae <- extract.clade(Caprimulgiformes, node=getMRCA(Caprimulgiformes, tip=c("Aegotheles_tatei", "Steatornis_caripensis")), root.edge=1)

#plot(Aegothelidae, root.edge=TRUE)

# Delete subspecies

Aegothelidae <- drop.tip(Aegothelidae, tip=c(
"Aegotheles_cristatus_major",
"Aegotheles_bennettii_bennettii",
"Aegotheles_bennettii_wiedenfeldi",
"Aegotheles_bennettii_plumiferus",
"Aegotheles_albertisi_salvadori",
"Aegotheles_albertisi_albertisi"
))


# Elevate affinis to species status

Aegothelidae$tip.label[which(Aegothelidae$tip.label =="Aegotheles_bennettii_affinis")] <- "Aegotheles_affinis"


# Get crown age from Backbone Tree
CladeAge <- getMRCAage(BackTre, tip=c("Aegotheles_insignis", "Steatornis_caripensis"))

# Build the calibration table
CALIB <- makeChronosCalib(Aegothelidae, node="root", age.min=CladeAge)

# Estimate chronogram
AegothelidaeC <- chronos(Aegothelidae, model="clock", calibration=CALIB, control=chronos.control(dual.iter.max=100))

#plot(AegothelidaeC); axisPhylo()

Ae2 <- drop.tip(AegothelidaeC, tip="Steatornis_caripensis", trim.internal=FALSE, collapse.singles=FALSE)

#plot(Ae2)

# time of MRCA Aegothelidae-Steatornithidae
max(branching.times(Ae2))

# MRCA time of Aegothelidae-Apodi

getMRCAage(BackTre, tip=c("Aegotheles_insignis", "Chaetura_pelagica"))

# Difference

max(branching.times(Ae2)) - getMRCAage(BackTre, tip=c("Aegotheles_insignis", "Chaetura_pelagica"))

# Delete diference to Aegothelidae stem

Ae2$edge.length[1] <- Ae2$edge.length[1] - (max(branching.times(Ae2)) - getMRCAage(BackTre, tip=c("Aegotheles_insignis", "Chaetura_pelagica")))

#plot(Ae2)

# Change the name of Aegotheles_insignis in the Backbone tre

BBtree26$tip.label[which(BBtree26$tip.label == "Aegotheles_insignis")] <- "Aegotheles_insignisDELETE"

# Bind the Aegothelidae tree at the Aegotheles_insignis stem node

BBtree27 <- bind.tree(BBtree26, Ae2, where=getMRCA(BBtree26, tip=c("Aegotheles_insignisDELETE", "Chaetura_pelagica")))

# Delete Aegotheles_insignisDELETE from the tree

BBtree27 <- drop.tip(BBtree27, tip=c("Aegotheles_insignisDELETE"))


#plot(BBtree27, cex=0.2, show.tip.label=FALSE)

# write.tree(BBtree27, file="BBtree27.tre")





### APODI ###

# Double check hummingbird tree

Apodi0 <- read.tree("PyPHLAWD Trees/Subtrees/Apodiformes_iqtree_partitioned-ML_aLRT.treefile")

Apodi0 <- midpoint(Apodi0)

Apodi0 <- ladderize(Apodi0, right=FALSE)

plot(Apodi0, show.tip.label=FALSE)

Apodidae <- extract.clade(Apodi0, node=getMRCA(Apodi0, tip=c("Streptoprocne_zonaris", "Chaetura_pelagica")))

#plot(Apodidae)

Trochilidae <- extract.clade(Apodi0, node=getMRCA(Apodi0, tip=c("Topaza_pella", "Chlorostilbon_lucidus")))

#plot(Trochilidae, cex=0.3)

# Delete Subspecies

cat(Apodi0$tip.label)

Apodi <- drop.tip(Apodi0, tip=c(
"Collocalia_esculenta_cyanoptila",
"Collocalia_linchi_dedii",
"Collocalia_linchi_linchi",
"Collocalia_affinis_cyanoptila",
"Collocalia_esculenta_marginata",
"Collocalia_esculenta_bagobo",
"Collocalia_esculenta_esculenta",
"Collocalia_esculenta_spilura",
"Collocalia_vanikorensis",
"Apus_apus_apus",
"Chaetura_vauxi_vauxi",
"Chaetura_vauxi_richmondi",
"Chaetura_vauxi_aphanes",
"Chaetura_chapmani_chapmani",
"Chaetura_brachyura",
"Chaetura_andrei",
"Chaetura_brachyura_ocypetes",
"Chaetura_brachyura_cinereocauda",
"Chaetura_cinereiventris_guianensis",
"Chaetura_cinereiventris_sclateri",
"Chaetura_spinicaudus",
"Chaetura_spinicaudus_aethalea",
"Chaetura_spinicaudus_aetherodroma",
"Phaethornis_bourcieri_major",
"Phaethornis_bourcieri_bourcieri",
"Coeligena_helianthea_helianthea",
"Coeligena_bonapartei_bonapartei",
"Coeligena_helianthea_tamai",
"Coeligena_bonapartei_consita",
"Metallura_tyrianthina_tyrianthina",
"Metallura_tyrianthina_smaragdinicollis",
"Metallura_tyrianthina_septentrionalis",
"Calliphlox_evelynae_evelynae",
"Thalurania_colombica_colombica",
"Chlorostilbon_lucidus_pucherani"))



# Delete subspecif epithets

Apodi$tip.label[which(Apodi$tip.label == "Collocalia_uropygialis_albidior")] <- "Collocalia_uropygialis"
Apodi$tip.label[which(Apodi$tip.label == "Chaetura_chapmani_viridipennis")] <- "Chaetura_viridipennis" # Considered a valid species by Marin (1997) but see Chesser et al. (2018)

Apodi$tip.label[which(Apodi$tip.label == "Chaetura_andrei_meridionalis")] <- "Chaetura_meridionalis"

Apodi$tip.label[which(Apodi$tip.label == "Chaetura_andrei_andrei")] <- "Chaetura_andrei"

Apodi$tip.label[which(Apodi$tip.label == "Chaetura_brachyura_brachyura")] <- "Chaetura_brachyura"

Apodi$tip.label[which(Apodi$tip.label == "Chaetura_cinereiventris_phaeopygos")] <- "Chaetura_phaeopygos"   # in a different clade so separate

Apodi$tip.label[which(Apodi$tip.label == "Chaetura_spinicaudus_spinicaudus")] <- "Chaetura_spinicaudus"

Apodi$tip.label[which(Apodi$tip.label == "Amazilia_saucerottei_hoffmanni")] <- "Amazilia_hoffmanni" # species leve taxon
Apodi$tip.label[which(Apodi$tip.label == "Cynanthus_latirostris_lawrencei")] <- "Cynanthus_lawrencei" # species leve taxon


anyDuplicated(Apodi$tip.label)
Apodi$tip.label[anyDuplicated(Apodi$tip.label)]


# Get crown age from Backbone Tree
CladeAge <- getMRCAage(BackTre, tip=c("Chaetura_pelagica", "Colibri_coruscans"))

# Build the calibration table
CALIB <- makeChronosCalib(Apodi, node="root", age.min=CladeAge)

# Estimate chronogram
ApodiC <- chronos(Apodi, model="discrete", calibration=CALIB, control=chronos.control(dual.iter.max=100, nb.rate.cat = 3))

#plot(ApodiC, show.tip.label=FALSE); axisPhylo()



# Clip crown clade from Backbone tree

BBtree28 <- drop.tip(BBtree27, tip=c("Hemiprocne_mystacea", "Streptoprocne_zonaris", "Phaethornis_griseogularis", "Topaza_pella", "Florisuga_mellivora"), trim.internal=TRUE)

BBtree28 <- drop.tip(BBtree28, tip=c("Chaetura_pelagica", "Colibri_coruscans"), trim.internal=FALSE)

# Add the new chronogram

BBtree28 <- bind.tree(BBtree28, ApodiC, where=Ntip(BBtree28))


#plot(BBtree28, cex=0.2, show.tip.label=FALSE)

#write.tree(BBtree28, file="BBtree28.tre")







##################
### GRUIFORMES ###
##################


Grui <- read.tree("PyPHLAWD Trees/Subtrees/Gruiformes_iqtree_partitioned-ML_aLRT.treefile")

Grui <- midpoint(Grui)

Grui <- ladderize(Grui, right=FALSE)

#plot(Grui, cex=0.5, no.margin=TRUE)


### Gruoidea ###

Gruoidea <- extract.clade(Grui, node=getMRCA(Grui, tip=c("Psophia_leucoptera", "Grus_americana")))

#plot(Gruoidea, cex=0.9)

# Delete subspecies

Gruoidea <- drop.tip(Gruoidea, tip=c(
"Antigone_antigone_sharpii",
"Antigone_antigone_antigone",
"Antigone_antigone_gillae",
"Antigone_canadensis_canadensis",
"Antigone_canadensis_pulla",
"Antigone_canadensis_rowani",
"Antigone_canadensis_tabida",
"Antigone_canadensis_pratensis",
"Balearica_regulorum_gibbericeps"))

Gruoidea <- ladderize(Gruoidea, right=FALSE)


# Get crown age from Backbone Tree
CladeAge <- getMRCAage(BackTre, tip=c("Psophia_crepitans", "Grus_canadensis"))

# Build the calibration table
CALIB <- makeChronosCalib(Gruoidea, node="root", age.min=CladeAge)

# Estimate chronogram
GruoideaC <- chronos(Gruoidea, model="clock", calibration=CALIB, control=chronos.control(dual.iter.max=100))

#plot(GruoideaC); axisPhylo()


# Clip crown clade from Backbone tree

BBtree29 <- drop.tip(BBtree28, tip=c("Aramus_guarauna"))

BBtree29 <- drop.tip(BBtree29, tip=c("Psophia_crepitans", "Grus_canadensis"), trim.internal=FALSE)

# Add the new chronogram

BBtree29 <- bind.tree(BBtree29, GruoideaC, where=Ntip(BBtree29))

#plot(BBtree29, cex=0.2, show.tip.label=FALSE)

#write.tree(BBtree29, file="BBtree29.tre")
# BBtree29 <- read.tree(file="BBtree29.tre")



### Ralloidea ###
# phPhlawd tree is more complete than J.C. Garcia-R, G.C. Gibb, S.A. Trewick Deep global evolutionary radiation in birds: Diversification and trait evolution in the cosmopolitan bird family Rallidae Mol. Phylogenet. Evol., 81 (2014), pp. 96-108 

Rallo <- extract.clade(Grui, node=getMRCA(Grui, tip=c("Heliornis_fulica", "Gallinula_chloropus")))

# Delete subspecies

Rallo <- drop.tip(Rallo, tip=c(
"Gallirallus_philippensis_philippensis",
"Gallirallus_philippensis_assimilis",
"Gallirallus_philippensis_mellori",
"Gallirallus_philippensis_sethsmithi",
"Gallirallus_philippensis_yorki",
"Gallirallus_philippensis_philippensis",
"Gallirallus_torquatus_sulcirostris",
"Gallirallus_torquatus_celebensis",
"Gallirallus_australis_australis",
"Gallirallus_australis_greyi",
"Rallus_aquaticus_aquaticus",
"Rallus_elegans_elegans",
"Rallus_crepitans_saturatus",
"Rallus_crepitans_insularum",
"Rallus_crepitans_waynei",
"Rallus_crepitans_scottii",
"Rallus_crepitans_crepitans",
"Rallus_crepitans_caribaeus",
"Rallus_crepitans_coryi",
"Rallus_crepitans_leucophaeus",
"Rallus_elegans_ramsdeni",
"Rallus_obsoletus_obsoletus",
"Rallus_obsoletus_beldingi",
"Rallus_obsoletus_levipes",
"Rallus_obsoletus_yumanensis",
"Rallus_obsoletus_rhizophorae",
"Rallus_longirostris_phelpsi",
"Rallus_longirostris_cypereti",
"Rallus_longirostris_dillonripleyi",
"Rallus_longirostris_margaritae",
"Rallus_limicola_limicola",
"Gallinula_chloropus_chloropus",
"Gallinula_galeata_sandvicensis",
"Gallinula_chloropus_guami",
"Porphyrio_poliocephalus_poliocephalus",
"Porphyrio_melanotus_melanopterus",
"Porphyrio_melanotus_melanotus",
"Porphyrio_melanotus_bellus",
"Porphyrio_porphyrio_caledonicus",
"Porphyrio_melanotus_samoensis",
"Porphyrio_porphyrio_vitiensis",
"Porphyrio_porphyrio_porphyrio",
"Nesoclopeus_woodfordi_immaculatus"
))

# Delete subspecies names 

Rallo$tip.label[which(Rallo$tip.label == "Gallirallus_torquatus_torquatus")] <- "Gallirallus_torquatus"

Rallo$tip.label[which(Rallo$tip.label == "Canirallus_kioloides_kioloides")] <- "Canirallus_kioloides"


Rallo <- ladderize(Rallo, right=FALSE)

#plot(Rallo, cex=0.6, no.margin=TRUE)

# Get crown age from Backbone Tree
CladeAge <- getMRCAage(BackTre, tip=c("Rallus_limicola", "Heliornis_fulica"))

# Build the calibration table
CALIB <- makeChronosCalib(Rallo, node="root", age.min=CladeAge)

# Estimate chronogram
RalloC <- chronos(Rallo, model="discrete", calibration=CALIB, control=chronos.control(dual.iter.max=100, nb.rate.cat=3))

#plot(RalloC, cex=0.6, no.margin=TRUE); axisPhylo()


# Clip crown clade from Backbone tree

BBtree30 <- drop.tip(BBtree29, tip=c("Sarothrura_insularis"))

BBtree30 <- drop.tip(BBtree30, tip=c("Rallus_limicola", "Heliornis_fulica"), trim.internal=FALSE)

# Add the new chronogram

BBtree30 <- bind.tree(BBtree30, RalloC, where=Ntip(BBtree30))

#plot(BBtree30, cex=0.2, show.tip.label=FALSE)

#write.tree(BBtree30, file="BBtree30.tre")
#write.nexus(BBtree30, file="BBtree30.nex")
#BBtree30 <- read.tree(file="BBtree30.tre")



#######################
### CHARADRIIFORMES ###
#######################

# ML tree from 
# Černý, D., & Natale, R. (2022). Comprehensive taxon sampling and vetted fossils help clarify the time tree of shorebirds (Aves, Charadriiformes). Molecular Phylogenetics and Evolution, 177, 107620.


Charadrii <- read.tree("Subclade Trees/Charadriiformes/concatenated.raxml.bestTree.tre")


Charadrii <- drop.tip(Charadrii, tip="OUTGROUP_Balearica_regulorum")

Charadrii <- ladderize(Charadrii, right=FALSE)


# Delete Family names

names <- strsplit(Charadrii$tip.label, split="_")

species <- character()

for(i in 1:length(names)) {

species[i] <- paste(names[[i]][2], names[[i]][3], sep="_")	
	
}

Charadrii$tip.label <- species

#plot(Charadrii, cex=0.5, no.margin=TRUE)

### Calirbating the whole tree is not working for some reason, optimiations are not converging
### Try to split the clade into two subclades ###

# Get crown age from Backbone Tree
CladeAge <- getMRCAage(BackTre, tip=c("Charadrius_vociferus", "Larus_marinus"))

# Build the calibration table
CALIB <- makeChronosCalib(Charadrii, node="root", age.min=CladeAge)

# Estimate chronogram
#CharadriiC <- chronos(Charadrii, model="discrete", calibration=CALIB, control=chronos.control(dual.iter.max=200, nb.rate.cat=3))

CharadriiC <- chronos(Charadrii, model="clock", calibration=CALIB, control=chronos.control(dual.iter.max=100))

#plot(CharadriiC, cex=0.5, show.tip.label=FALSE); axisPhylo()


# Clip crown clade from Backbone tree

BBtree31 <- drop.tip(BBtree30, tip=c(
"Burhinus_grallarius",
"Chionis_minor",
"Pluvianellus_socialis",
"Pluvianus_aegyptius",
"Ibidorhyncha_struthersii",
"Haematopus_ater",
"Recurvirostra_americana",
"Scolopax_rusticola",
"Pedionomus_torquatus",
"Thinocorus_orbignyianus",
"Rostratula_benghalensis",
"Jacana_jacana",
"Turnix_sylvatica",
"Dromas_ardeola",
"Stiltia_isabella",
"Catharacta_skua",
"Alca_torda",
"Uria_aalge",
"Alle_alle",
"Gygis_alba",
"Anous_tenuirostris",
"Rynchops_niger",
"Sterna_hirundo"
))

BBtree31 <- drop.tip(BBtree31, tip=c("Charadrius_vociferus", "Larus_marinus"), trim.internal=FALSE)

# Add the new chronogram

BBtree31 <- bind.tree(BBtree31, CharadriiC, where=Ntip(BBtree31))

#plot(BBtree31, cex=0.2, show.tip.label=FALSE)

#write.tree(BBtree31, file="BBtree31.tre")
#BBtree31 <- read.tree(file="BBtree31.tre")


### Calirbating the whole tree is not working for some reason, optimiations are not converging
### Try to split the clade into two subclades ###




#########################
###                   ###
###   AQUAEORNITHES   ###
###                   ###
#########################


###############
### Gavidae ###
###############


# Extract Clade
Gaviidae <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Gavia_stellata", "Gavia_immer")))

#plot(Gaviidae)

# Get crown age from Backbone Tree
CladeAge <- getMRCAage(BackTre, tip=c("Gavia_stellata", "Gavia_immer"))

# Build the calibration table
CALIB <- makeChronosCalib(Gaviidae, node="root", age.min=CladeAge)

# Estimate chronogram
GaviidaeC <- chronos(Gaviidae, model="clock", calibration=CALIB, control=chronos.control(dual.iter.max=100))

#plot(GaviidaeC); axisPhylo()

# Clip taxa from backbone

BBtree32 <- drop.tip(BBtree31, tip=c("Gavia_stellata", "Gavia_immer"), trim.internal=FALSE)

# Add the new chronogram

BBtree32 <- bind.tree(BBtree32, GaviidaeC, where=Ntip(BBtree32))

#plot(BBtree32, cex=0.2, show.tip.label=FALSE)




####################
### Spheniscidae ###
####################

# Extract Clade
Spheniscidae <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Aptenodytes_patagonicus", "Spheniscus_magellanicus")))

#plot(Spheniscidae)

# Get crown age from Backbone Tree
CladeAge <- getMRCAage(BackTre, tip=c("Aptenodytes_forsteri", "Spheniscus_humboldti"))

# Build the calibration table
CALIB <- makeChronosCalib(Spheniscidae, node="root", age.min=CladeAge)

# Estimate chronogram
SpheniscidaeC <- chronos(Spheniscidae, model="clock", calibration=CALIB, control=chronos.control(dual.iter.max=100))

SpheniscidaeC <- ladderize(SpheniscidaeC, right=FALSE)

#plot(SpheniscidaeC); axisPhylo()


# Clip taxa from backbone

BBtree33 <- drop.tip(BBtree32, tip="Eudyptes_chrysolophus")

BBtree33 <- drop.tip(BBtree33, tip=c("Aptenodytes_forsteri", "Spheniscus_humboldti"), trim.internal=FALSE)

# Add the new chronogram

BBtree33 <- bind.tree(BBtree33, SpheniscidaeC, where=Ntip(BBtree33))

#plot(BBtree33, cex=0.2, show.tip.label=FALSE)




#########################
### PROCELLARIIFORMES ###
#########################

Proce <- read.tree("PyPHLAWD Trees/Subtrees/Procellariiformes_iqtree_partitioned-ML_aLRT.treefile")


#plot(Proce, no.margin=TRUE, cex=0.6); nodelabels()

Proce <- root(Proce, node=209)

#plot(Proce, no.margin=TRUE, cex=0.6)


# Delete subspecies and Diomedea exulans that seems wrong.

Proce <- drop.tip(Proce, tip=c(
"Diomedea_exulans",
"Thalassarche_bulleri_platei",
"Thalassarche_chlororhynchos",
"Thalassarche_bulleri", 
"Thalassarche_cauta",
"Thalassarche_melanophris",
"Hydrobates_furcatus_furcatus",              
"Hydrobates_furcatus_plumbeus",
"Hydrobates_pelagicus_melitensis"  ,         
"Hydrobates_pelagicus_pelagicus",
"Hydrobates_tethys_tethys" ,      
"Hydrobates_tethys_kelsalli",
 "Fregetta_grallaria_segethi" ,          
"Fregetta_grallaria_titan",                  
 "Fregetta_grallaria_grallaria",              
 "Fregetta_grallaria_leucogaster",            
 "Fregetta_tropica_tropica",      
"Fregetta_tropica_melanoleuca",
 "Oceanites_oceanicus_exasperatus",
 "Pelagodroma_marina_dulciae" ,          
 "Pelagodroma_marina_albiclunis",
 "Calonectris_diomedea_diomedea" ,
  "Pseudobulweria_rostrata_rostrata",          
"Pseudobulweria_rostrata_trouessarti",
"Fulmarus_glacialis_glacialis",
 "Fulmarus_glacialis_auduboni",               
 "Fulmarus_glacialis_rodgersi",
"Hydrobates_leucorhous_leucorhous"  ,        
 "Hydrobates_leucorhous_chapmani"  
))

Proce <- ladderize(Proce, right=FALSE)

#plot(Proce, cex=0.5)


# Delete subspecific epithets

Proce $tip.label[which(Proce $tip.label == "Thalassarche_melanophrys_melanophrys")] <- "Thalassarche_melanophrys"

Proce $tip.label[which(Proce $tip.label == "Thalassarche_cauta_cauta")] <- "Thalassarche_cauta"

Proce $tip.label[which(Proce $tip.label == "Thalassarche_bulleri_bulleri")] <- "Thalassarche_bulleri"

Proce $tip.label[which(Proce $tip.label == "Thalassarche_chlororhynchos_chlororhynchos")] <- "Thalassarche_chlororhynchos"

Proce $tip.label[which(Proce $tip.label == "Procellaria_aequinoctialis_conspicillata")] <- "Procellaria_conspicillata"


#plot(Proce, cex=0.6, no.margin=TRUE)


# Get crown age from Backbone Tree
CladeAge <- getMRCAage(BackTre, tip=c("Thalassarche_bulleri", "Puffinus_creatopus"))

# Build the calibration table
CALIB <- makeChronosCalib(Proce, node="root", age.min=CladeAge)

# Estimate chronogram
ProceC <- chronos(Proce, model="clock", calibration=CALIB, control=chronos.control(dual.iter.max=100))

#plot(ProceC, cex=0.5); axisPhylo()


# Clip taxa from backbone

BBtree34 <- drop.tip(BBtree33, tip=c("Oceanites_oceanicus", "Oceanodroma_castro", "Pelecanoides_magellani"))

BBtree34 <- drop.tip(BBtree34, tip=c("Thalassarche_bulleri", "Puffinus_creatopus"), trim.internal=FALSE)

# Add the new chronogram

BBtree34 <- bind.tree(BBtree34, ProceC, where=Ntip(BBtree34))

#plot(BBtree34, cex=0.2, show.tip.label=FALSE)


#write.tree(BBtree34, file="BBtree34.tre")



##################
### Ciconiidae ###
##################


Ciconi <- read.tree("PyPHLAWD Trees/Subtrees/Ciconiiformes_iqtree_partitioned-ML_aLRT.treefile")

Ciconi <- extract.clade(Ciconi, node=getMRCA(Ciconi, tip=c("Ciconia_ciconia", "Scopus_umbretta")))

#plot(Ciconi)


# Get crown age from Backbone Tree
CladeAge <- getMRCAage(BackTre, tip=c("Ciconia_abdimii", "Scopus_umbretta"))

# Build the calibration table
CALIB <- makeChronosCalib(Ciconi, node="root", age.min=CladeAge)

# Estimate chronogram
CiconiC <- chronos(Ciconi, model="clock", calibration=CALIB, control=chronos.control(dual.iter.max=100))

CiconiC <- drop.tip(CiconiC, tip="Ciconia_ciconia_ciconia")

CiconiC <- drop.tip(CiconiC, tip="Balaeniceps_rex")

#plot(CiconiC, cex=0.9); axisPhylo()

# Change name of the outgroup so it is not repreated in the backboen tree
CiconiC$tip.label[which(CiconiC$tip.label=="Scopus_umbretta")] <- "OG"

# Change the name of the strok in the backbone tree so it is not repreated in the subclade

BBtree34$tip.label[which(BBtree34$tip.label=="Ciconia_abdimii")] <- "stork"

## Graft subclade into backbone

BBtree34 <- bind.tree(BBtree34, CiconiC, where=getMRCA(BBtree34, tip=c("stork", "Scopus_umbretta")))


# Clip taxa from backbone

BBtree35 <- drop.tip(BBtree34, tip=c("stork", "OG"))

#plot(BBtree35, cex=0.2, show.tip.label=FALSE)

#write.tree(BBtree35, file="BBtree35.tre")


######################
### PELECANIFORMES ###
######################


Peleca <- read.tree("PyPHLAWD Trees/Subtrees/Pelecaniformes_iqtree_partitioned-ML_aLRT.treefile")

Peleca <- midpoint(Peleca)

#plot(Peleca, cex=0.5)


### Ardeidae ###

Ardeidae <- extract.clade(Peleca, node=getMRCA(Peleca, tip=c("Tigrisoma_lineatum", "Ardea_alba")))

#plot(Ardeidae, no.margin=TRUE)

# Delete subspecies

Ardeidae <- drop.tip(Ardeidae, tip=c(
"Ardea_herodias_fannini",
"Ardea_herodias_wardi",
"Ardea_herodias_herodias",        
"Ardea_cinerea_cinerea",
"Ardea_cinerea_jouyi",
"Ardeola_bacchus/speciosa"
))

Ardeidae <- ladderize(Ardeidae, right=FALSE)

# Get crown age from Backbone Tree
CladeAge <- getMRCAage(BackTre, tip=c("Tigrisoma_lineatum", "Egretta_tricolor"))

# Build the calibration table
CALIB <- makeChronosCalib(Ardeidae, node="root", age.min=CladeAge)

# Estimate chronogram
ArdeidaeC <- chronos(Ardeidae, model="clock", calibration=CALIB, control=chronos.control(dual.iter.max=100))

#plot(ArdeidaeC); axisPhylo()


# Clip taxa from backbone

BBtree36 <- drop.tip(BBtree35, tip=c("Tigrisoma_lineatum", "Egretta_tricolor"), trim.internal=FALSE)

# Add the new chronogram

BBtree36 <- bind.tree(BBtree36, ArdeidaeC, where=Ntip(BBtree36))

#plot(BBtree36, cex=0.2, show.tip.label=FALSE)


### Suloidea ###

Suloidea <- extract.clade(Peleca, node=getMRCA(Peleca, tip=c("Fregata_minor", "Sula_sula")))


Suloidea <- drop.tip(Suloidea, tip=c(
"Phalacrocorax_aristotelis_aristotelis",
"Phalacrocorax_carbo_sinensis",         
"Phalacrocorax_carbo_novaehollandiae",
"Phalacrocorax_carbo_carbo" ,
"Anhinga_anhinga_leucogaster",         
"Anhinga_anhinga_anhinga"
))

Suloidea <- ladderize(Suloidea, right=FALSE)

#plot(Suloidea, no.margin=TRUE, cex=0.7)

# Get crown age from Backbone Tree
CladeAge <- getMRCAage(BackTre, tip=c("Fregata_minor", "Sula_sula"))

# Build the calibration table
CALIB <- makeChronosCalib(Suloidea, node="root", age.min=CladeAge)

# Estimate chronogram
SuloideaC <- chronos(Suloidea, model="clock", calibration=CALIB, control=chronos.control(dual.iter.max=100))

#plot(SuloideaC); axisPhylo()


# Clip taxa from backbone

BBtree37 <- drop.tip(BBtree36, tip="Phalacrocorax_carbo")
BBtree37 <- drop.tip(BBtree37, tip="Anhinga_anhinga")

BBtree37 <- drop.tip(BBtree37, tip=c("Fregata_minor", "Sula_sula"), trim.internal=FALSE)

# Add the new chronogram

BBtree37 <- bind.tree(BBtree37, SuloideaC, where=Ntip(BBtree37))

#plot(BBtree37, cex=0.2, show.tip.label=FALSE)

#write.tree(BBtree37, file="BBtree37.tre")


### Threskiornithidae ###

Thresk <- keep.tip(Peleca, tip=c(
"Pelecanus_crispus",
"Eudocimus_albus",
"Eudocimus_ruber",
"Theristicus_caerulescens",
"Theristicus_caudatus",
"Theristicus_melanopis",
"Phimosus_infuscatus",
"Mesembrinibis_cayennensis",
"Plegadis_chihi",
"Plegadis_falcinellus",
"Plegadis_ridgwayi",
"Bostrychia_hagedash",
"Geronticus_eremita",
"Geronticus_calvus",
"Nipponia_nippon",
"Threskiornis_molucca",
"Platalea_ajaja",
"Platalea_flavipes",
"Platalea_regia",
"Platalea_minor",
"Platalea_leucorodia",
"Platalea_alba",
"Threskiornis_aethiopicus"))

#plot(Thresk)

# Get crown age from Backbone Tree
CladeAge <- getMRCAage(BackTre, tip=c("Theristicus_melanopis", "Scopus_umbretta"))

# Build the calibration table
CALIB <- makeChronosCalib(Thresk, node="root", age.min=CladeAge)

# Estimate chronogram
ThreskC <- chronos(Thresk, model="clock", calibration=CALIB, control=chronos.control(dual.iter.max=100))

#plot(ThreskC); axisPhylo()

BBtree37$tip.label[which(BBtree37$tip.label == "Theristicus_melanopis")] <- "ibis"

ThreskC $tip.label[which(ThreskC $tip.label == "Pelecanus_crispus")] <- "OG"


# Add the new chronogram

BBtree38 <- bind.tree(BBtree37, ThreskC, where=getMRCA(BBtree37, tip=c("ibis", "Scopus_umbretta")))


# Clip taxa from backbone

BBtree38 <- drop.tip(BBtree38, tip=c("ibis", "OG"))

#plot(BBtree38, cex=0.2, show.tip.label=FALSE)

#write.tree(BBtree38, file="BBtree38.tre")



### Pelecanidae ###

Pelecanidae <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Pelecanus_crispus", "Scopus_umbretta")))

Pelecanidae <- drop.tip(Pelecanidae, tip=c(
"Balaeniceps_rex" ,
"Pelecanus_occidentalis_carolinensis",
"Pelecanus_occidentalis_californicus"
))

Pelecanidae <- ladderize(Pelecanidae, right=FALSE)

#plot(Pelecanidae)

# Get crown age from Backbone Tree
CladeAge <- getMRCAage(BackTre, tip=c("Pelecanus_erythrorhynchos", "Scopus_umbretta"))

# Build the calibration table
CALIB <- makeChronosCalib(Pelecanidae, node="root", age.min=CladeAge)

# Estimate chronogram
PelecanidaeC <- chronos(Pelecanidae, model="clock", calibration=CALIB, control=chronos.control(dual.iter.max=100))

#plot(PelecanidaeC); axisPhylo()

BBtree38$tip.label[which(BBtree38$tip.label == "Pelecanus_erythrorhynchos")] <- "pelican"

PelecanidaeC $tip.label[which(PelecanidaeC $tip.label == "Scopus_umbretta")] <- "OG"


# Add the new chronogram

BBtree39 <- bind.tree(BBtree38, PelecanidaeC, where=getMRCA(BBtree38, tip=c("pelican", "Scopus_umbretta")))


# Clip taxa from backbone

BBtree39 <- drop.tip(BBtree39, tip=c("pelican", "OG"))

#plot(BBtree39, cex=0.2, show.tip.label=FALSE)

#write.tree(BBtree39, file="BBtree39.tre")

#BBtree39 <- read.tree(file="BBtree39.tre")




########################
### ACCIPITRIMORPHAE ###
########################



Accip <- read.tree("PyPHLAWD Trees/Subtrees/Accipitriformes_iqtree_partitioned-ML_aLRT.treefile")

#plot(Accip, cex=0.5)

cat(Accip$tip.label)

# Delete subspecies and Cathartidae

Accip <- drop.tip(Accip, tip=c(
"Circus_aeruginosus_harterti",
"Circus_aeruginosus_aeruginosus",
"Buteo_lineatus_extimus",
"Buteo_lineatus_lineatus",
"Buteo_lineatus_elegans",
"Buteo_lineatus_alleni",
"Buteo_jamaicensis_borealis",
"Buteo_japonicus_japonicus",
"Buteo_buteo_vulpinus",
"Buteo_rufinus_cirtensis",
"Buteo_buteo_buteo",
"Buteo_lagopus_sanctijohannis",
"Buteo_lagopus_kamtschatkensis",
"Buteo_jamaicensis_harlani",
"Buteo_albonotatus_albonotatus",
"Buteo_platypterus_platypterus",
"Buteo_nitidus_costaricensis",
"Buteo_nitidus_pallida",
"Leucopternis_albicollis_albicollis",
"Leucopternis_albicollis_ghiesbreghti",
"Leucopternis_albicollis_costaricensis",
"Leucopternis_albicollis_williaminae",
"Geranoaetus_melanoleucus",
"Geranoaetus_melanoleucus_australis",
"Geranoaetus_melanoleucus_melanoleucus",
"Buteo_polyosoma_polyosoma",
"Buteo_albicaudatus_hypospodius",
"Buteo_albicaudatus_colonus",
"Parabuteo_unicinctus_harrisi",
"Parabuteo_unicinctus_unicinctus",
"Leucopternis_princeps_zimmeri",
"Leucopternis_princeps_princeps",
"Buteogallus_urubitinga_ridgwayi",
"Buteogallus_urubitinga_urubitinga",
"Buteogallus_subtilis_bangsi",
"Busarellus_nigricollis_leucocephalus",
"Milvus_migrans_lineatus",
"Milvus_migrans_govinda",
"Milvus_migrans_affinis",
"Milvus_migrans_migrans",
"Milvus_migrans_parasitus",
"Milvus_milvus_fasciicauda",
"Milvus_milvus_milvus",
"Haliastur_indus_girenera",
"Haliaeetus_leucocephalus_alascanus",
"Haliaeetus_albicilla_albicilla",
"Haliaeetus_pelagicus_pelagicus",
"Gyps_indicus_indicus",
"Gyps_fulvus_fulvus",
"Gyps_fulvus_fulvescens",
"Spilornis_cheela_hoya",
"Spilornis_cheela_burmanicus",
"Spilornis_cheela_perplexus",
"Gampsonyx_swainsonii_leonae",
"Cathartes_melambrotus",
"Cathartes_burrovianus",
"Cathartes_aura",
"Sarcoramphus_papa",
"Vultur_gryphus",
"Gymnogyps_californianus",
"Coragyps_atratus",
"Pandion_haliaetus_haliaetus",
"Aviceda_subcristata_guerneyi",
"Pernis_ptilorhynchus_orientalis",
"Elanoides_forficatus_yetapa",
"Chondrohierax_uncinatus_wilsonii",
"Neophron_percnopterus_ginginianus",
"Gypaetus_barbatus_aureus",
"Clanga_pomarina_pomarina",
"Clanga_pomarina_hastata",
"Ictinaetus_malayensis_malayensis",
"Hieraaetus_morphnoides_morphnoides",
"Hieraaetus_pennatus_pennatus",
"Hieraaetus_pennatus_minisculus",
"Lophotriorchis_kienerii_formosus",
"Lophotriorchis_kienerii_kienerii",
"Spizaetus_ornatus_vicarius",
"Spizaetus_tyrannus_serus",
"Nisaetus_nipalensis_orientalis",
"Spizaetus_cirrhatus_limnaeetus",
"Accipiter_trivirgatus_palawanus",
"Accipiter_trivirgatus_extimus",
"Accipiter_badius_cenchroides",
"Accipiter_francesiae_pusillus",
"Accipiter_francesiae_brutus",
"Accipiter_francesiae_griveaudi",
"Accipiter_virgatus_quagga",
"Accipiter_virgatus_confusus",
"Accipiter_toussenelii_canescens",
"Accipiter_toussenelii_toussenelii",
"Accipiter_toussenelii_macroscelides",
"Accipiter_tachiro_sparsimfasciatus",
"Accipiter_tachiro_pembaensis",
"Accipiter_striatus_striatus",
"Accipiter_striatus_perobscurus",
"Accipiter_striatus_velox",
"Accipiter_nisus_nisosimilis",
"Accipiter_rufiventris_rufiventris",
"Accipiter_nisus_nisus",
"Accipiter_nisus_wolterstorffi",
"Accipiter_gentilis_marginatus",
"Accipiter_gentilis_albidus",
"Accipiter_gentilis_arrigonii",
"Accipiter_gentilis_buteoides",
"Accipiter_gentilis_fujiyamae",
"Accipiter_gentilis_schvedowi",
"Accipiter_gentilis_gentilis",
"Accipiter_melanoleucus_temminckii",
"Accipiter_melanoleucus_melanoleucus",
"Accipiter_gentilis_laingi",
"Accipiter_gentilis_atricapillus",
"Circus_cyaneus_cyaneus",
"Kaupifalco_monogrammicus_meridionalis"))

# Delete subspecific epithets:

Accip$tip.label

Accip$tip.label[which(Accip$tip.label == "Macheiramphus_alcinus_anderssoni")] <- "Macheiramphus_alcinus"

Accip$tip.label[which(Accip$tip.label == "Accipiter_francesiae_francesiae")] <- "Accipiter_francesiae"

Accip$tip.label[which(Accip$tip.label == "Accipiter_cirrocephalus_cirrocephalus")] <- "Accipiter_cirrocephalus"

Accip$tip.label[which(Accip$tip.label == "Accipiter_novaehollandiae_rufoschistaceus")] <- "Accipiter_novaehollandiae"


# Root with Sagitarius
Accip <- root(Accip, outgroup="Sagittarius_serpentarius")

# Delete outgroup
Accip <- drop.tip(Accip, tip="Sagittarius_serpentarius")

#plot(Accip, cex=0.4)

# Get crown age from Backbone Tree
CladeAge <- getMRCAage(BackTre, tip=c("Pandion_haliaetus", "Elanus_caeruleus"))

# Build the calibration table
CALIB <- makeChronosCalib(Accip, node="root", age.min=CladeAge)

# Estimate chronogram
AccipC <- chronos(Accip, model="discrete", calibration=CALIB, control=chronos.control(dual.iter.max=100, nb.rate.cat=3))

#plot(AccipC, cex=0.4); axisPhylo()


# Clip taxa from backbone

BBtree40 <- drop.tip(BBtree39, tip="Buteo_jamaicensis")

BBtree40 <- drop.tip(BBtree40, tip=c("Pandion_haliaetus", "Elanus_caeruleus"), trim.internal = FALSE)

## Graft subclade into backbone

BBtree40 <- bind.tree(BBtree40, AccipC, where=Ntip(BBtree40))



#plot(BBtree40, cex=0.2, show.tip.label=FALSE)

#write.tree(BBtree40, file="BBtree40.tre")



### Cathartidae ###


Cathartidae <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Cathartes_aura", "Vultur_gryphus")))

#plot(Cathartidae)

# Correct condor name in backbone tree

BBtree40$tip.label[which(BBtree40 $tip.label == "Vultur_gryhpus")] <- "Vultur_gryphus"


# Get crown age from Backbone Tree
CladeAge <- getMRCAage(BBtree40, tip=c("Cathartes_aura", "Vultur_gryphus"))

# Build the calibration table
CALIB <- makeChronosCalib(Cathartidae, node="root", age.min=CladeAge)

# Estimate chronogram
CathartidaeC <- chronos(Cathartidae, model="clock", calibration=CALIB, control=chronos.control(dual.iter.max=100))

#plot(CathartidaeC); axisPhylo()


# Clip taxa from backbone


BBtree41 <- drop.tip(BBtree40, tip=c("Cathartes_aura", "Vultur_gryphus"), trim.internal = FALSE)

## Graft subclade into backbone

BBtree41 <- bind.tree(BBtree41, CathartidaeC, where=Ntip(BBtree41))



#plot(BBtree41, cex=0.2, show.tip.label=FALSE)

#write.tree(BBtree41, file="BBtree41.tre")






### Strigidae ###

Strigidae <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Tyto_alba", "Strix_varia")))

#plot(Strigidae, cex=0.5)

cat(Strigidae$tip.label)


# Delete subspecies 

Strigidae <- drop.tip(Strigidae, tip=c(
"Athene_cunicularia_cunicularia",
"Athene_noctua_plumipes",
"Athene_noctua_indigena",
"Athene_noctua_vidalii",
"Aegolius_acadicus_brooksi",
"Aegolius_acadicus_acadicus",
"Asio_flammeus_flammeus",
"Asio_otus_otus",
"Asio_otus_canariensis",
"Asio_otus_wilsonianus",
"Ketupa_zeylonensis_orientalis",
"Strix_uralensis_japonica",
"Strix_occidentalis_caurina",
"Strix_woodfordii_nigricantior",
"Ninox_boobook_ocellata",
"Ninox_boobook_boobook",
"Ninox_boobook_halmaturina",
"Ninox_novaeseelandiae_leucopsis",
"Ninox_novaeseelandiae_undulata",
"Ninox_novaeseelandiae_novaeseelandiae",
"Ninox_connivens_connivens",
"Ninox_philippensis_centralis",
"Ninox_philippensis_philippensis",
"Tyto_novaehollandiae_calabyi",
"Tyto_javanica_sumbaensis",
"Tyto_javanica_delicatula",
"Tyto_furcata_hellmayri",
"Tyto_javanica_lulu",
"Tyto_alba_ernesti",
"Tyto_alba_affinis",
"Tyto_alba_guttata",
"Tyto_alba_alba",
"Tyto_alba_erlangeri",
"Tyto_furcata_tuidara",
"Tyto_furcata_bargei",
"Tyto_furcata_pratincola",
"Tyto_alba_pratincola",
"Tyto_tenebricosa_arfaki",
"Tyto_tenebricosa_tenebricosa"
))


# Delete subspecific epithets

Strigidae$tip.label[which(Strigidae $tip.label == "Tyto_sororcula_sororcula")] <- "Tyto_sororcula"

Strigidae$tip.label[which(Strigidae $tip.label == "Tyto_furcata_furcata")] <- "Tyto_furcata"

plot(Strigidae)
is.rooted(Strigidae)

# Get crown age from Backbone Tree
CladeAge <- getMRCAage(BBtree41, tip=c("Tyto_alba", "Athene_cunicularia"))

# Build the calibration table
CALIB <- makeChronosCalib(phy=Strigidae, node="root", age.min=CladeAge)

# Estimate chronogram
StrigidaeC <- chronos(Strigidae, model="discrete", calibration=CALIB, control=chronos.control(dual.iter.max=100, nb.rate.cat=2))

#plot(StrigidaeC); axisPhylo()


# Clip taxa from backbone


BBtree42 <- drop.tip(BBtree41, tip=c("Phodilus_badius", "Ninox_novaeseelandiae", "Strix_occidentalis"))

BBtree42 <- drop.tip(BBtree42, tip=c("Tyto_alba", "Athene_cunicularia"), trim.internal = FALSE)

## Graft subclade into backbone

BBtree42 <- bind.tree(BBtree42, StrigidaeC, where=Ntip(BBtree42))



#plot(BBtree42, cex=0.2, show.tip.label=FALSE, edge.width=0.2)

#write.tree(BBtree42, file="BBtree42.tre")

#BBtree42  <- read.tree(file="BBtree42.tre")



######################
### CORACIIMORPHAE ###
######################


### Coliidae ###

# Complete analytical


### Trogonidae ###

# New UCE tree suggests that Apaloderma is basal and all three biogeograpic regions are monophyletic but tree not availble.
#Oliveros, C. H., Andersen, M. J., Hosner, P. A., Mauck III, W. M., Sheldon, F. H., Cracraft, J., & Moyle, R. G. (2020). Rapid Laurasian diversification of a pantropical bird family during the Oligocene–Miocene transition. Ibis, 162(1), 137-152.
# In the meantime, use phyPHLAWD tree:

Trogo <- read.tree("PyPHLAWD Trees/Subtrees/Trogoniformes_iqtree_partitioned-ML_aLRT.treefile")

#plot(Trogo); nodelabels()

Trogo <- midpoint(Trogo)


# Get crown age from Backbone Tree
CladeAge <- getMRCAage(BackTre, tip=c("Trogon_personatus", "Appaloderma_vittatum"))

# Build the calibration table
CALIB <- makeChronosCalib(Trogo, node="root", age.min=CladeAge)

# Estimate chronogram
TrogoC <- chronos(Trogo, model="clock", calibration=CALIB, control=chronos.control(dual.iter.max=100))

#plot(TrogoC); axisPhylo()

# Clip taxa from backbone
BBtree43 <- drop.tip(BBtree42, tip="Apalharpactes")

BBtree43 <- drop.tip(BBtree43, tip=c("Trogon_personatus", "Appaloderma_vittatum"), trim.internal = FALSE)

## Graft subclade into backbone

BBtree43 <- bind.tree(BBtree43, TrogoC, where=Ntip(BBtree43))

#plot(BBtree43, cex=0.2, show.tip.label=FALSE)

#write.tree(BBtree43, file="BBtree43.tre")



### Bucerotiformes ###

# From PhlawdTree that is rooted and includes Upupidae and Phoneniculidae

Bucero <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Upupa_epops", "Buceros_bicornis")))

#plot(Bucero, no.margin=TRUE, cex=0.8)


# Delete semingly erroneous sequences: 
Bucero <- drop.tip(Bucero, tip="Penelopides_panini")

# Delete subspecific epithets

Bucero $tip.label[which(Bucero $tip.label == "Tockus_monteiri_majoriae")] <- "Tockus_monteiri"

Bucero $tip.label[which(Bucero $tip.label == "Tockus_pallidirostris_neumanni")] <- "Tockus_pallidirostris"

Bucero $tip.label[which(Bucero $tip.label == "Tockus_fasciatus_fasciatus" )] <- "Tockus_fasciatus" 

Bucero $tip.label[which(Bucero $tip.label == "Penelopides_panini_panini")] <- "Penelopides_panini"

Bucero $tip.label[which(Bucero $tip.label == "Penelopides_manillae_subniger")] <- "Penelopides_manillae"

Bucero $tip.label[which(Bucero $tip.label == "Penelopides_affinis_affinis")] <- "Penelopides_affinis"

Bucero $tip.label[which(Bucero $tip.label == "Penelopides_exarhatus_exarhatus")] <- "Penelopides_exarhatus"

Bucero $tip.label[which(Bucero $tip.label == "Bycanistes_fistulator_duboisi" )] <- "Bycanistes_fistulator"

Bucero $tip.label[which(Bucero $tip.label == "Bycanistes_subcylindricus_subquadratus")] <- "Bycanistes_subcylindricus"

Bucero $tip.label[which(Bucero $tip.label == "Bycanistes_subcylindricus_subquadratus")] <- "Bycanistes_subcylindricus"


# Delete subspecies
Bucero <- drop.tip(Bucero, tip=c(
"Tockus_leucomelas_parvior",
"Tockus_flavirostris_flavirostris",
"Tropicranus_albocristatus_cassini",
"Rhyticeros_plicatus_ruficollis",
"Buceros_rhinoceros_silvestris",
"Aceros_corrugatus_rugosus",
"Buceros_hydrocorax_hydrocorax",       
"Buceros_hydrocorax_mindanensis",
"Buceros_hydrocorax_semigaleatus"))


#plot(Bucero, cex=0.8)


# Get ages from Backbone Tree
RootAge <- getMRCAage(BackTre, tip=c("Upupa_epops", "Buceros_bicornis"))

BuceroAge <- getMRCAage(BackTre, tip=c("Bucorvus_leadbeateri", "Buceros_bicornis"))

UpupiAge <- getMRCAage(BackTre, tip=c("Upupa_epops", "Phoeniculus_purpureus"))


# Build the calibration table

#plot(Bucero, cex=0.8); nodelabels()

CALIB <- makeChronosCalib(Bucero, node=c(55,56, 106), age.min=c(RootAge, BuceroAge, UpupiAge))

# Estimate chronogram
BuceroC <- chronos(Bucero, model="discrete", calibration=CALIB, control=chronos.control(dual.iter.max=100, nb.rate.cat=3))

#plot(BuceroC); axisPhylo()

# Clip taxa from backbone
BBtree44 <- drop.tip(BBtree43, tip=c("Bucorvus_leadbeateri", "Phoeniculus_purpureus"))

BBtree44 <- drop.tip(BBtree44, tip=c("Upupa_epops", "Buceros_bicornis"), trim.internal = FALSE)

## Graft subclade into backbone

BBtree44 <- bind.tree(BBtree44, BuceroC, where=Ntip(BBtree44))

#plot(BBtree44, cex=0.2, show.tip.label=FALSE)

#write.tree(BBtree44, file="BBtree44.tre")


# BBtree44  <- read.tree(file="BBtree44.tre")





#####################
### CORACIIFOMRES ###
#####################

### Read the McCullough et al 2019 tree

coraci <- read.tree("Subclade Trees/Coraciiformes/figS1_200_90perc_part_100bs_bipart.tre")

#plot(coraci, cex=0.5)

# Delete extra taxa

coraci <- drop.tip(coraci, tip=c(
"merops_viridis_americanus_kunhm27326", 
"halcyon_senegalensis_kunhm19982", 
"halcyon_smyrnensis_kunhm23257",
"dacelo_leachii_anwcb56231", 
"todiramphus_tutus_mhngt3", 
"todiramphus_veneratus_youngi_mhngt4",
"todiramphus_chloris_pealei_uwbm89771",
"todiramphus_chloris_manuae_kunhm107630",
"alcedo_atthis_wam22763",
"alcedo_atthis_kunhm23444",
"alcedo_euryzona_uwbm67494",
"ceyx_erithaca_wam22719",
"ceyx_melanurus_samarensis_kunhm31654",
"ceyx_melanurus_mindanensis_kunhm18127",
"momotus_momota_microstephanus_fmnh456580",
"hornbill",
"trogon", 
"cuckoo",
"picoides_pubescens_ncbi32059"
))


# Elevate to species status

coraci$tip.label[which(coraci$tip.label == "ceyx_meeki_pallidus_kunhm32075")] <- "ceyx_pallidus"



# Delete specimen info from species name

names <- strsplit(coraci$tip.label, split="_")

species <- character()

for(i in 1:length(names)) {

species[i] <- paste(names[[i]][1], names[[i]][2], sep="_")	
	
}

species

# Capitalize first letter of genera

# Capitalizing function (from the toupper help):
# Capitalizes the first letter of each word

simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1, 1)), substring(s, 2),
          sep = "", collapse = " ")}
          
species <- sapply(species, simpleCap)

# Check for duplicated names

anyDuplicated(species)

# species[anyDuplicated(species)]

cbind(species, coraci$tip.label) 

#Substitute old labels by new labels

coraci$tip.label <- species

#plot(coraci, cex=0.5)



# Get age from Backbone Tree

RootAge <- getMRCAage(BackTre, tip=c("Merops_pusillus", "Alcedo_leucogaster"))

# Build the calibration table

CALIB <- makeChronosCalib(coraci, age.min=RootAge)

# Estimate chronogram
coraciC <- chronos(coraci, model="discrete", calibration=CALIB, control=chronos.control(dual.iter.max=100, nb.rate.cat=3))

coraciC <- ladderize(coraciC, right = FALSE)

#plot(coraciC, cex=0.5); axisPhylo()


# Clip taxa from backbone
BBtree45 <- drop.tip(BBtree44, tip=c("Nyctyornis_amicus", "Brachypteracias_leptosomus","Coracias_caudata", "Todus_angustirostris", "Momotus_momota", "Hylomanes_momotula","Chloroceryle_americana", "Halcyon_malimbica"))

BBtree45 <- drop.tip(BBtree45, tip=c("Merops_pusillus", "Alcedo_leucogaster"), trim.internal = FALSE)

## Graft subclade into backbone

BBtree45 <- bind.tree(BBtree45, coraciC, where=Ntip(BBtree45))

#plot(BBtree45, cex=0.2, show.tip.label=FALSE)

#write.tree(BBtree45, file="BBtree45.tre")

# BBtree45  <- read.tree(file="BBtree45.tre")



##################
### Piciformes ###
##################


### Galbuli ###

Galbuli <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Galbula_dea", "Bucco_tamatia")))

#plot(Galbuli)

# Brachygalba polyphyletic (but may be correct) and with suspiciously long branches.

# Delete subspecies #

Galbuli <- drop.tip(Galbuli, tip="Malacoptila_striata_minor")

#plot(Galbuli)

# Get age from Backbone Tree
RootAge <- getMRCAage(BBtree45, tip=c("Bucco_capensis", "Galbula_albirostris"))

# Build the calibration table

CALIB <- makeChronosCalib(Galbuli, age.min=RootAge)

# Estimate chronogram
GalbuliC <- chronos(Galbuli, model="discrete", calibration=CALIB, control=chronos.control(dual.iter.max=100, nb.rate.cat=3))

GalbuliC <- ladderize(GalbuliC, right = FALSE)

#plot(GalbuliC, cex=0.5); axisPhylo()


# Clip taxa from backbone

BBtree46 <- drop.tip(BBtree45, tip=c("Bucco_capensis", "Galbula_albirostris"), trim.internal = FALSE)

## Graft subclade into backbone

BBtree46 <- bind.tree(BBtree46, GalbuliC, where=Ntip(BBtree46))

#plot(BBtree46, cex=0.2, show.tip.label=FALSE)

#write.tree(BBtree46, file="BBtree46.tre")

#BBtree46  <- read.tree(file="BBtree46.tre")



### Picides ###

# Picidae and Indicatoridae
# Shakya, S. B., Fuchs, J., Pons, J. M., & Sheldon, F. H. (2017). Tapping the woodpecker tree for evolutionary insight. Molecular Phylogenetics and Evolution, 116, 182-191.

Picides <- read.tree("Subclade Trees/Piciformes/1-s2.0-S1055790317300477-mmc7.tre")

#plot(Picides)


# Delete duplicated specimens

Picides <- drop.tip(Picides, tip=c(
"Jynx_torquilla_DD33",
"Picumnus_albosquamatus_Z16",
"Picumnus_exilis_Guy0302",
"Picumnus_aurifrons_LSU46136",
"Picumnus_lafresnayi_LSU25473", # keeping FMNH 473979, closer to the type locality
"Dinopium_benghalense_Z7",
"Dinopium_javanense_M78",
"Picus_viridis_JF3718",
"Campethera_nivosa_JF3553",
"Dryocopus_javensis_W46",
"Melanerpes_formicivorus_CAS92857",
"Yungipicus_moluccensis_AMNH9620",
"Dendrocoptes_medius_DD34",
"Dendropicos_griseocephalus_DB30",
"Leuconotopicos_villosus_LSU19818",
"Veniliornis_callonotus_Z25",
"Veniliornis_passerinus_LSU26094"
))


# Delete Specimen numbers #

Nnames <- strsplit(Picides$tip.label, split="_")

Species <- character()

for(i in 1:length(Nnames)) {
	 Species[i] <- paste(Nnames[[i]][1], Nnames[[i]][2], sep="_")
	 }

Picides$tip.label <- Species

which(duplicated(Species))

#plot(Picides, cex=0.5)

# Get age from Backbone Tree
RootAge <- getMRCAage(BBtree46, tip=c("Indicator_variegatus", "Melanerpes_carolinus"))

# Build the calibration table

CALIB <- makeChronosCalib(Picides, age.min=RootAge)

# Estimate chronogram
PicidesC <- chronos(Picides, model="discrete", calibration=CALIB, control=chronos.control(dual.iter.max=100, nb.rate.cat=2))

PicidesC <- ladderize(PicidesC, right = FALSE)

#plot(PicidesC, cex=0.5); axisPhylo()


# Clip taxa from backbone

BBtree47 <- drop.tip(BBtree46, tip="Picumnus_cirratus", trim.internal = TRUE)

BBtree47 <- drop.tip(BBtree47, tip=c("Indicator_variegatus", "Melanerpes_carolinus"), trim.internal = FALSE)


## Graft subclade into the backbone

BBtree47 <- bind.tree(BBtree47, PicidesC, where=Ntip(BBtree47))

#plot(BBtree47, cex=0.2, show.tip.label=FALSE)

#write.tree(BBtree47, file="BBtree47.tre")



### Ramphastides ###

# phPhlawd tree: Piciformes_iqtree_partitioned-ML_aLRT.treefile
# Eventially change for tree of Otrow et al. 2023

Ramphastides <- read.tree("PyPHLAWD Trees/Subtrees/Piciformes_iqtree_partitioned-ML_aLRT.treefile")

Ramphastides <- midpoint(Ramphastides)

#plot(Ramphastides, cex=0.5)

Ramphastides <- extract.clade(Ramphastides, node=getMRCA(Ramphastides, tip=c("Megalaima_virens", "Ramphastos_toco")))


#plot(Ramphastides, cex=0.5)

cat(Ramphastides$tip.label) 

# Delete Subspecies #

Ramphastides <- drop.tip(Ramphastides, tip=c("Capito_auratus_transilens",
"Capito_auratus_aurantiicinctus",
"Capito_auratus_punctatus",
"Capito_auratus_nitidior",
"Capito_auratus_auratus",
"Capito_auratus_orosae",
"Capito_auratus_amazonicus",
"Capito_auratus_insperatus",
"Capito_maculicoronatus_rubrilateralis",
"Eubucco_bourcierii_tucinkae",
"Eubucco_bourcierii_aequitorialis",
"Selenidera_reinwardtii_reinwardtii",
"Andigena_hypoglauca_lateralis",
"Andigena_hypoglauca_hypoglauca",
"Pteroglossus_mariae",
"Pteroglossus_azara_flavirostris",
"Pteroglossus_azara_azara",
"Pteroglossus_bitorquatus_sturmii",
"Pteroglossus_castanotis_castanotis",
"Pteroglossus_castanotis_australis",
"Pteroglossus_aracari_atricollis",
"Pteroglossus_aracari_aracari",
"Pteroglossus_inscriptus_humboldti",
"Pteroglossus_inscriptus_inscriptus",
"Pteroglossus_bitorquatus_reichenowi",
"Ramphastos_sulfuratus_sulfuratus",
"Ramphastos_sulfuratus_brevicarinatus",
"Ramphastos_vitellinus_vitellinus",
"Ramphastos_ambiguus_ambiguus",
"Ramphastos_tucanus_tucanus",
"Ramphastos_toco_albogularis",
"Pogoniulus_bilineatus_fischeri",
"Pogoniulus_bilineatus_conciliator",
"Pogoniulus_bilineatus_bilineatus",
"Pogoniulus_bilineatus_jacksoni",
"Pogoniulus_bilineatus_mfumbiri",
"Pogoniulus_chrysoconus_extoni",
"Pogoniulus_pusillus_pusillus",
"Pogoniulus_chrysoconus_chrysoconus",
"Pogoniulus_pusillus_affinis",
"Megalaima_franklinii_franklinii",
"Megalaima_asiatica_davisoni",
"Megalaima_asiatica_asiatica",
"Megalaima_oorti_annamensis",
"Megalaima_oorti_sini",
"Megalaima_oorti_faber",
"Megalaima_oorti_oorti",
"Megalaima_haemacephala_mindanensis",
"Megalaima_haemacephala_intermedia",
"Megalaima_haemacephala_celestinoi",
"Megalaima_haemacephala_haemacephala",
"Megalaima_haemacephala_indica",
"Megalaima_haemacephala_rosea",
"Megalaima_haemacephala_delica"))

#plot(Ramphastides, cex=0.5)

# Elevate subspecies to species

Ramphastides$tip.label[which(Ramphastides$tip.label == "Capito_maculicoronatus_maculicoronatus")] <- "Capito_maculicoronatus"

Ramphastides$tip.label[which(Ramphastides$tip.label == "Aulacorhynchus_caeruleogularis_caeruleogularis")] <- "Aulacorhynchus_caeruleogularis"

Ramphastides$tip.label[which(Ramphastides$tip.label == "Aulacorhynchus_atrogularis_atrogularis")] <- "Aulacorhynchus_atrogularis"

Ramphastides$tip.label[which(Ramphastides$tip.label == "Andigena_nigrirostris_spilorhynchus")] <- "Andigena_nigrirostris"

duplicated(Ramphastides$tip.label)
which(duplicated(Ramphastides$tip.label))

# Get age from Backbone Tree
RootAge <- getMRCAage(BBtree47, tip=c("Megalaima_oorti", "Capito_niger"))

# Build the calibration table

CALIB <- makeChronosCalib(Ramphastides, age.min=RootAge)

# Estimate chronogram
RamphastidesC <- chronos(Ramphastides, model="discrete", calibration=CALIB, control=chronos.control(dual.iter.max=100, nb.rate.cat=3))

RamphastidesC <- ladderize(RamphastidesC, right = FALSE)

#plot(RamphastidesC, cex=0.5); axisPhylo()


# Clip taxa from backbone

BBtree48 <- drop.tip(BBtree47, tip=c("Lybius_hirsutus", "Semnornis_frantzii", "Pteroglossus_aracari"), trim.internal = TRUE)

BBtree48 <- drop.tip(BBtree48, tip=c("Megalaima_oorti", "Capito_niger"), trim.internal = FALSE)


## Graft subclade into the backbone

BBtree48 <- bind.tree(BBtree48, RamphastidesC, where=Ntip(BBtree48))

#plot(BBtree48, cex=0.2, show.tip.label=FALSE)

#write.tree(BBtree48, file="BBtree48.tre")


#####################
### FALCONIFORMES ###
#####################


Falcon <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Herpetotheres_cachinnans", "Falco_peregrinus")))

#plot(Falcon, cex=0.7, no.margin=TRUE)

cat(Falcon$tip.label, sep="', '")

# Delete subspecies 

Falcon <- drop.tip(Falcon, tip=c(
"Milvago_chimango_chimango",
"Caracara_plancus_plancus",
"Microhierax_erythrogenys_meridionalis",
"Microhierax_caerulescens_caerulescens",
"Polihierax_semitorquatus_semitorquatus",
"Falco_sparverius_sparverius",
"Falco_sparverius_cinnamominus",
"Falco_tinnunculus_canariensis",
"Falco_tinnunculus_tinnunculus",
"Falco_tinnunculus_interstinctus",
"Falco_tinnunculus_dacotiae",
"Falco_tinnunculus_rufescens",
"Falco_cenchroides_cenchroides",
"Falco_vespertinus_vespertinus",
"Falco_femoralis_pichinchae",
"Falco_subbuteo_subbuteo",
"Falco_chicquera_horsbrughi",
"Falco_cherrug_cyanopus",
"Falco_cherrug_milvipes",
"Falco_biarmicus_feldeggii",
"Falco_cherrug_cherrug",
"Falco_biarmicus_biarmicus",
"Falco_chicquera_chicquera",
"Falco_rusticolus_candidans",
"Falco_rusticolus_obsoletus",
"Falco_peregrinus_brookei",
"Falco_peregrinus_nesiotes",
"Falco_pelegrinoides_babylonicus",
"Falco_peregrinus_tundrius",
"Falco_peregrinus_minor",
"Falco_peregrinus_calidus",
"Falco_peregrinus_pealei",
"Falco_peregrinus_macropus",
"Falco_peregrinus_peregrinator",
"Falco_peregrinus_cassini",
"Falco_pelegrinoides_pelegrinoides",
"Falco_peregrinus_anatum",
"Falco_peregrinus_peregrinus",
"Falco_peregrinus_japonensis",
"Falco_columbarius_columbarius",
"Falco_columbarius_aesalon",
"Falco_columbarius_richardsonii"))

#plot(Falcon, cex=0.7, no.margin=TRUE)

# Get age from Backbone Tree
RootAge <- getMRCAage(BBtree48, tip=c("Micrastur_gilvicollis", "Falco_peregrinus"))

# Build the calibration table

CALIB <- makeChronosCalib(Falcon, age.min=RootAge)

# Estimate chronogram
FalconC <- chronos(Falcon, model="discrete", calibration=CALIB, control=chronos.control(dual.iter.max=100, nb.rate.cat=3))

FalconC <- ladderize(FalconC, right=FALSE)

#plot(FalconC, cex=0.5); axisPhylo()

# Clip taxa from backbone
BBtree49 <- drop.tip(BBtree48, tip=c("Daptrius_ater", "Microhierax_caerulescens"), trim.internal = TRUE)

BBtree49 <- drop.tip(BBtree49, tip=c("Micrastur_gilvicollis", "Falco_peregrinus"), trim.internal = FALSE)


## Graft subclade into the backbone

BBtree49 <- bind.tree(BBtree49, FalconC, where=Ntip(BBtree49))

#plot(BBtree49, cex=0.2, show.tip.label=FALSE)

#write.tree(BBtree49, file="BBtree49.tre")
#BBtree49 <- read.tree(file="BBtree49.tre")




######################
### PSITTACIFORMES ###
######################


# Picidae and Indicatoridae
# Smith, B. T., Merwin, J., Provost, K. L., Thom, G., Brumfield, R. T., Ferreira, M., ... & Joseph, L. (2022). Phylogenomic analysis of the parrots of the world distinguishes artifactual from biological sources of gene tree discordance. Systematic Biology.

Psitta <- read.tree("Subclade Trees/Psittaciformes/mafft-nexus-clean-75p.charsets.treefile")

# Root and delete outgroup
Psitta <- root(Psitta, outgroup=c("Calyptomena_viridis_b36377", "Icterus_cucullatus_b24601", "Caracara_cheriway_b5328"))

Psitta <- drop.tip(Psitta, tip=c("Calyptomena_viridis_b36377", "Icterus_cucullatus_b24601", "Caracara_cheriway_b5328"))

#plot(Psitta, cex=0.5)

# Delete Specimen numbers #

Nnames <- strsplit(Psitta$tip.label, split="_")

Species <- character()

for(i in 1:length(Nnames)) {
	 Species[i] <- paste(Nnames[[i]][1], Nnames[[i]][2], sep="_")
	 }

Psitta$tip.label <- Species

# Replace name with typo:

Psitta$tip.label[which(Psitta$tip.label=="Myiopsitta_fuchsi")] <- "Myiopsitta_luchsi"

Psitta$tip.label[which(Psitta$tip.label=="Psittacula_calthorpae")] <- "Psittacula_calthrapae"

Psitta$tip.label[which(Psitta$tip.label=="Eunymphicus_cornotus")] <- "Eunymphicus_cornutus"

Psitta$tip.label[which(Psitta$tip.label=="Eunymphicus_uvaensis")] <- "Eunymphicus_uvaeensis"




which(duplicated(Psitta))

#plot(Psitta, cex=0.5)

# Get age from Backbone Tree
RootAge <- getMRCAage(BBtree49, tip=c("Nestor_notabilis", "Psittacus_erithacus"))

# Build the calibration table

CALIB <- makeChronosCalib(Psitta, age.min=RootAge)

# Estimate chronogram
PsittaC <- chronos(Psitta, model="discrete", calibration=CALIB, control=chronos.control(dual.iter.max=100, nb.rate.cat=3))

PsittaC <- ladderize(PsittaC, right = FALSE)

#plot(PsittaC, cex=0.5); axisPhylo()


# Clip taxa from backbone

BBtree50 <- drop.tip(BBtree49, tip=c("Alisterus_scapularis","Myiopsitta_monachus", "Calyptorhynchus_funereus"), trim.internal = TRUE)

BBtree50 <- drop.tip(BBtree50, tip=c("Nestor_notabilis", "Psittacus_erithacus"), trim.internal = FALSE)


## Graft subclade into the backbone

BBtree50 <- bind.tree(BBtree50, PsittaC, where=Ntip(BBtree50))

#plot(BBtree50, cex=0.2, show.tip.label=FALSE)

#write.tree(BBtree50, file="BBtree50.tre")


# Check Tree #

#Phylogenetic tree with 3743 tips and 3740 internal nodes.
# Unrooted; includes branch lengths.

#Check for duplicated tips
anyDuplicated(BBtree50$tip.label)

# Check branch lengths
summary(BBtree50$edge.length)

# Check for ultrametricide
is.ultrametric(BBtree50)

# Make the tree exactly unrametric:

BBtree50 <- nnls.tree(cophenetic(BBtree50), BBtree50, method="ultrametric", rooted=TRUE)

is.ultrametric(BBtree50)

#write.tree(BBtree50, file="BBtree50.tre")




###############################
### Non-passerines Big Tree ###
###############################

BBtreeNP <- BBtree50

plot(BBtreeNP, show.tip.label=FALSE, edge.width=0.2)

plot(BBtreeNP, type="fan", show.tip.label=FALSE, edge.width=0.2)


#write.tree(BBtreeNP, file="BBtreeNP.tre")


##############################
##### Graft The Two Trees ####
##############################

# Load trees
BBtreeNP <- read.tree(file="BBtreeNP.tre")
BBtreeP <- read.tree(file="PasseriformesBigTree Michelle/PasBackTre37.tre")

plot(BBtreeP, show.tip.label=FALSE, edge.width=0.2)

# Identify grafting node: the MRCA (crown node) of Passeriformes in the main tree

MRCApass <- getMRCA(BBtreeNP, tip=c("Acanthisitta_chloris","Passer_montanus"))

# Delete Passeriformes in main tree
# First identify the species to delete (the desdendants from MRCApass)
# This gets the indices
Pass <- Descendants(BBtreeNP, node=MRCApass)[[1]]

# Now get the names
Pass <- BBtreeNP$tip.label[Pass]

# Drop all but the first two (Acanthisitta and Pitta)
todrop <- Pass[-which(Pass=="Acanthisitta_chloris")]
todrop <- todrop[-which(Pass=="Pitta_sordida")]

BBtreeNP <- drop.tip(BBtreeNP, tip=todrop)

plot(BBtreeNP, show.tip.label=FALSE, edge.width=0.2)

# Drop the remaining two without deleting internal branch
BBtreeNP <- drop.tip(BBtreeNP, tip= c("Acanthisitta_chloris", "Pitta_sordida"), trim.internal=FALSE)

plot(BBtreeNP, show.tip.label=FALSE, edge.width=0.2)

# Graft Passeriformes tree into the main tree
BBtree <- bind.tree(BBtreeNP, BBtreeP, where=Ntip(BBtreeNP))

plot(BBtree, show.tip.label=FALSE, edge.width=0.2)

# Delete nodel labels

BBtree$node.label <- NULL

# Delete Repeated species and subspecies
# Non-passerines (eventually incorporate in the code above)
BBtree <- drop.tip(BBtree, tip=c("Piaya_minuta", "Collocalia_vanikorensis","Chroicocephalus_poiocephalus", "Larus_mongolicus", "Larus_heuglini", "Larus_vegae", "Phalacrocorax_lucidus", "Eudyptula_albosignata","Puffinus_haurakiensis", "Puffinus_haurakiensis","Puffinus_kermadecensis", "Puffinus_temptator",
"Puffinus_dichrous",
"Puffinus_polynesiae",
"Puffinus_colstoni",
"Puffinus_nicolae",
"Phoebastria_nihonus",
"Diomedea_gibsoni",
"Heliactin_bilophus", # H. cornutus represents the species in the tree
"Tyto_sororcula_cayelii",
"Athene_lilith",
"Strix_omanensis", # junior synonym
"Harpactes_ardens_ardens",
"Harpactes_ardens_herberti",
"Megalaima_chrysopogon",    # duplicated
"Megalaima_monticola",     # duplicated
"Baillonius_bailloni",  # duplicated
"Hydropsalis_brasiliana",
"Polyplancta_aurescens",
"Calliphlox_mitchellii",
"Ketupa_flavipes",
"Megalaima_mystacophanos"  # duplicated
))


# Passerines
BBtree <- drop.tip(BBtree, tip=c("Euchrepomis_humeralis2","Laniarius_erlangeri","Phylloscopus_brehmi", "Cettia_robustipes", "Pnoepyga_albiventer_albiventer", "Pnoepyga_pusilla_pusilla", "Zosterops_capensis", "Zosterops_ugiensis", "Hedydipna", "Nectarinia_sperata","Hedydipna_collaris", "Motacilla_yarrellii", "Anthus_nairobi", "Hemignathus_munroi", "Hemignathus_kauaiensis", "Spizella_taverneri", "Lophorina_magnifica", "Spreo_fischeri", "Neocossyphus_fraseri", "Thryothorus_leucotis", "Cantorchilus_longirostris", "Bradypterus_seebohmi", "Anthreptes_fraseri", "Nectarinia_senegalensis", "Nectarinia_minulla", "Nectarinia_olivacea"))


### Substitute faulty names and subspecies names:

# Non-passerines (incorporate above)
BBtree$tip.label[BBtree$tip.label == "Halcyon_badial"] <- "Halcyon_badia"
BBtree$tip.label[BBtree$tip.label == "Halcyon_pileataioz18911"] <- "Halcyon_pileata"

# Passerines
BBtree$tip.label[BBtree$tip.label == "Anthreptes_griseigularis_birgitae"] <- "Anthreptes_griseigularis"

BBtree$tip.label[BBtree$tip.label == "Orthonyx_novaeguineae_victoriana"] <- "Orthonyx_novaeguineae"

BBtree$tip.label[BBtree$tip.label == "Calandrella_cinerea_williamsi"] <- "Calandrella_cinerea"

BBtree$tip.label[BBtree$tip.label == "Scotocerca_inquieta_saharae"] <- "Scotocerca_saharae"

BBtree$tip.label[BBtree$tip.label == "Arizelocichla_masukensis"] <- "Arizelocichla_masukuensis"

BBtree$tip.label[BBtree$tip.label == "Arizelocichla_kikuyensis"] <- "Arizelocichla_kikuyuensis"

BBtree$tip.label[BBtree$tip.label == "Vidua_orientalis_aucupum"] <- "Vidua_orientalis"

BBtree$tip.label[BBtree$tip.label == "Neochmia_ruficaudaU"] <- "Neochmia_ruficauda"
BBtree$tip.label[BBtree$tip.label == "Neochmia_modestaU"] <- "Neochmia_modesta"

BBtree$tip.label[BBtree$tip.label =="Glycifohia_undulata"] <- "Gliciphila_undulata"


# Ladderize 
BBtree <- ladderize(BBtree, right = FALSE)


### Save BBtree ###
# write.tree(BBtree, file="BBtree.tre")
### BBtree <- read.tree(file="BBtree.tre")


### PLOT BBtree ###

plot(BBtree, show.tip.label=FALSE, edge.width=0.15)

plot(BBtree, show.tip.label=FALSE, edge.width=0.3, direction="upward")

plot(BBtree, type="fan", show.tip.label=FALSE, edge.width=0.2)

ltt.plot(BBtree, log="y")

# write.tree(BBtree, file="BBtree.tre")
# BBtree <- read.tree("~/Library/CloudStorage/GoogleDrive-sclaramunt.uy@gmail.com/My Drive/BigBirdTree Calibration/BigBirdTree Assembly/BBtree.tre")


#######################################################
### Version 2: binarized and no zero branch lengths ###
#######################################################

is.binary(BBtree)
BBtree2 <- multi2di(BBtree)
is.binary(BBtree2)

summary(BBtree2$edge.length)
min(BBtree2$edge.length)

plot(density(BBtree2$edge.length), log="x"); rug(BBtree2$edge.length)

BBtree2$edge.length[BBtree2$edge.length < 0.01]

BBtree2$edge.length[BBtree2$edge.length < 0.01] <- 0.01

plot(BBtree2, show.tip.label=FALSE, edge.width=0.15)

# Force Stric Ultrametricity
library(phytools)
is.ultrametric(BBtree2)
BBtree2 <- force.ultrametric(BBtree2)
is.ultrametric(BBtree2)

plot(BBtree2, show.tip.label=FALSE, edge.width=0.15)

# write.tree(BBtree2, file="BBtree2.tre")


ltt.plot(BBtree2, log="y")


#################################
### Version 3: clean taxonomy ###
#################################

# As now, I cleaned the taxonomy a bit but still not using the Avonet eBirdClements column.

BBtree3 <- read.tree(file="BBtree2.tre")

# Delete duplicated species
BBtree3 <- drop.tip(BBtree3, tip=c("Collocalia_vanikorensis"))

# Load conversion table.
# This is an updated version of Hosner et al. .2021. I added about 700 unmatched names that were represented in our phylogenetic trees but not in Hosner's reconciliation table, probably because the names were already absent from GenBank.

Recon <- read.csv("~/GoogleDrive/BigBirdTree Calibration/BigBirdTree Assembly/GenBank_eBird_Clements_2022_taxonomic_reconciliationUPDATED.csv")

# Add a low dash to the spaces in the Recon table.

Recon$GenBank_name <- sapply(Recon$GenBank_name, chartr, old=" ", new="_")
Recon$eBirdClements2022_Names <- sapply(Recon$eBirdClements2022_Names, chartr, old=" ", new="_")
Recon$BigBird_Names <- sapply(Recon$BigBird_Names, chartr, old=" ", new="_")


# Substitute tree names by updated names

head(Recon)

# Find the eBird name that matches the tree name
# First find the pisition in the Recon table of the matching name in the tree

match.index <- match(BBtree3$tip.label, Recon$GenBank_name)

BBtree3$tip.label <- Recon$BigBird_Names[match.index]

# Check for duplicated names
BBtree3$tip.label[duplicated(BBtree3$tip.label)]

# Old Code
# Restore original names for those with a match but no assinged eBird names
# BBtree3$tip.label[which(BBtree3$tip.label == "")] <- BBtree2$tip.label[which(BBtree3$tip.label == "")]
# Keep the original name
# BBtree3$tip.label[anyDuplicated(BBtree3$tip.label)] <- BBtree2$tip.label[anyDuplicated(BBtree3$tip.label)]
# or delete de duplicated taxon
# BBtree3 <- drop.tip(BBtree3, tip= BBtree3$tip.label[anyDuplicated(BBtree3$tip.label)])

write.tree(BBtree3, file="BBtree3.tre")



###################################################
### Version 4 Clements/eBird 2022 species names ###
###################################################

BBtreeC2022 <- read.tree(file="BBtree2.tre")

BBtreeC2022 <- drop.tip(BBtreeC2022, tip=c("Collocalia_vanikorensis"))

# Substitute tree names with Clements/eBird 2022 names.

# Find the eBird name that matches the tree name
# First find the pisition in the Recon table of the matching name in the tree

match.index <- match(BBtreeC2022$tip.label, Recon$GenBank_name)

BBtreeC2022$tip.label <- Recon$BigBird_Names[match.index]

# Check for duplicated names
BBtreeC2022$tip.label[duplicated(BBtreeC2022$tip.label)]


write.tree(BBtreeC2022, file="BBtreeC2022.tre")
