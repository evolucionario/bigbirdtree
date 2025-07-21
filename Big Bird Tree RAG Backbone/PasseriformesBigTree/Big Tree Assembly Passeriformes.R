###############################################
### Assembly of the Big Passeriformes Tree  ###
###############################################
#install.packages("phylotools")

rm(list=ls())


library(ape)
library(phangorn)
library(phylotools)
library(stringr)


# new and cleared laptop
setwd("~/DIRECTORY")



##### FINAL tree: ####
PasBackTre37  <- read.tree(file="PasBackTre37.tre")


### Ancillary Functions ####

# getMRCAage returns the age of the MRCA of a pair of tips

getMRCAage <- function(phy, tip) {

	# Obtain branching times
	
	BTimes <- branching.times(phy)
	
	# identify the MRCA node number
	
	MRCA <- getMRCA(phy, tip)

	# substract Ntip to obtain the index of the node in the branching times result

	MRCAnode <- MRCA - Ntip(phy)
	
	return(BTimes[MRCAnode])	
	
}

#####################################################
### GraftTree: ADD A CLADE TO A PHYLOGENETIC TREE ###
#####################################################

#This function attaches a clade to a branch in another phylogenetic tree. It uses functions in the package "ape". It differs from bind.tree in that the clade is not attached to a tip or an internatl node but to a branch.

### Usage ###

#graft.tree(phy1, phy2, edge, age)

### Arguments ###

# phy1		an object of class "phylo"

# phy2		an object of class "phylo"

# edge		a vector of mode numeric or character specifying the branch to which phy2 will be attached to. The branch is specified based on terminal nodes indicated through the getMRCA function, thus a vector of tip labels.

# age		the time of origin of phy2. Default to NULL assumes that phy2 includes its complete root edge

# TO DO: make it work with trees without root edge and add root edge based on given age

### FUNCTION ###

graft.tree <- function(phy1, phy2, edge, age=NULL) {
  
  # Identify the edge of phy1 where ph2 will be attached to.
  
  if(length(edge) == 1) {
    # When a single tip is indicated, the edge is terminal 
    Node <- which(phy1$tip.label == edge) 
    mrcaage0 <- 0
  }
  else {
    # The edge is identified using the terminal node via specifyin MRCA of two tips
    Node <- getMRCA(phy1, tip= edge)
    mrcaage0 <- getMRCAage(phy1, edge)
  }
  
  # Calculate the grafting pisition along the branch
  # by substracting the terminal age of the target branch (0 if a terminal branch, MRCA age if subtending a clade) from the time of origin of phy2
  # If time of origin ("age") is not provided, the information can be extracted from the root edge of phy2, if present.
  
  if(is.null(age)) {
    if(is.null(phy2$root.edge)) {
      stop("Either a root edge for phy2 or a specified age is needed")
    } else {
      age <- max(branching.times(phy2)) + phy2$root.edge
    } 
  }	
  
  Pos <- age - mrcaage0
  
  new.tree <- bind.tree(phy1, phy2, where= Node, position=Pos)
  
  return(new.tree)
}









##### Load RAG Backbone Tree ######

BackTre <- read.tree("ChronogramNeornithes.tre")

plot(BackTre, edge.width=0.5, cex=0.1, no.margin=T, show.tip.label=F)

axisPhylo(cex.axis=0.7)

# Extract Passeriformes 

PasBackTre <- extract.clade(BackTre, node=getMRCA(BackTre, tip=c("Acanthisitta_chloris", "Passer_montanus")))


plot(PasBackTre, edge.width=0.5, cex=0.1)


##### Load phyPhlawd big tree ######

PhlawdTree <- read.tree("RAxML_bestTree.Aves_bait_13.tre") # ROM

PhlawdSpp <- as.data.frame(PhlawdTree$tip.label)

####### Load Passeroidea_iqtree_partitioned-ML_aLRT.tre ########
# help understanding the structure of this and which group to use as outgroup
# midpoint and ascending node order
Passeroidea <- read.tree("Passeroidea_iqtree_partitioned-ML_aLRT.tre")
spp <- as.data.frame(Passeroidea$tip.label)

Pas <- midpoint(Passeroidea)
PasSpp <- as.data.frame(Pas$tip.label)

plot(Pas, show.tip.label = F)
#write.tree(PasseroideaR, file="test.tre")


#### Tyranni ####
# load Emberizoidea code first!

# Load Harvey et al. 2020 tree
# Given the size, it's going to be difficult to calibrate with chronos so use a TreePL tree and rescale using a simple branch multiplication using the appropriate ratio. 
# can use the same strategy for some other trees if they are in the form of a calibrated Bayesian MCC tree

# Load the tree
Tyranni <- read.tree("T400F_AOS_Clements_sppnames.tre")

# Delete outgroup
Tyranni <- extract.clade(Tyranni, node=getMRCA(Tyranni, tip=c("Eurylaimus_javanicus","Tyrannus_tyrannus")))

plot(Tyranni, show.tip.label=FALSE)

### Determine scaling factor ###

#Obtain backbone age
CladeAge <- getMRCAage(PasBackTre, tip=c("Sapayoa_aenigma", "Tyrannus_tyrannus"))

#Obtain tree age
TreeAge <- getMRCAage(Tyranni, tip=c("Eurylaimus_javanicus","Tyrannus_tyrannus"))

# Convert branch lengths
Tyranni$edge.length <- Tyranni$edge.length*CladeAge/TreeAge

# Trim species in backbone

PasBackTre5 <- drop.tip(PasBackTre5, tip=c(
  "Pitta_sordida",
  "Smithornis_rufolateralis",
  "Psarisomus_dalhousiae",
  "Philepitta_castanea",
  "Pipra_coronata",
  "Cotinga_cayana",
  "Tityra_semifasciata",
  "Onychorhynchus_coronatus",
  "Oxyruncus_cristatus",
  "Rhynchocyclus_brevirostris",
  "Tachuris_rubrigastra",
  "Piprites_chloris",
  "Platyrinchus_coronatus",
  "Conopophaga_ardesiaca",
  "Melanopareia_torquata",
  "Terenura_sharpei",
  "Myrmornis_torquata",
  "Thamnophilus_nigrocinereus",
  "Scytalopus_magellanicus",
  "Grallaria_ruficapilla",
  "Formicarius_colma",
  "Sclerurus_mexicanus",
  "Dendrocolaptes_certhia",
  "Furnarius_rufus"), trim.internal=TRUE)

PasBackTre5 <- drop.tip(PasBackTre5, tip=c("Sapayoa_aenigma", "Tyrannus_tyrannus"),  trim.internal=FALSE)

# Graft Harvey et al tree into Backbone

PasBackTre6 <- bind.tree(PasBackTre5, Tyranni, where=Ntip(PasBackTre5))

plot(PasBackTre6, show.tip.label=FALSE)

### SAVE work ###

write.tree(PasBackTre6, file="PasBackTre6.tre")

# PasBackTre6  <- read.tree(file="PasBackTre6.tre")





#### Passeri ####


#### Climacteridae and Ptilonorhynchidae ####

### Get source TREE ###

# Use PhlawdTree tree, since there are no better phylogeny available in the literature

ClimacPtilo <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Climacteris_picumnus", "Ptilonorhynchus_violaceus")))

plot(ClimacPtilo)


### Manage TAXA and tipl labels ###

# Show tip labels to select subspecies to delete

ClimacPtilo$tip.label

# Can use cat() and BBedit to arrange lists for cut & paste

# Make sure that one sample remains to represent the species

# Atention to the A. crassirostris/melanotis complex
# delete the A. melanotis and keep the sample identified as nominate subspecies, to be sure

# Drop Subspecies

ClimacPtilo <- drop.tip(ClimacPtilo, tip=c(
"Ailuroedus_geislerorum_geislerorum",
"Ailuroedus_geislerorum_molestus",
"Ailuroedus_buccoides_buccoides",
"Ailuroedus_stonii_cinnamomeus",
"Ailuroedus_stonii_stonii",
"Ailuroedus_melanotis_facialis",
"Ailuroedus_melanotis",
"Ailuroedus_jobiensis_guttaticollis",
"Ailuroedus_jobiensis_jobiensis",
"Ailuroedus_arfakianus_misoliensis",
"Ailuroedus_crassirostris_joanae",
"Amblyornis_macgregoriae_nubicola",
"Chlamydera_guttata_carteri",
"Chlamydera_nuchalis_nuchalis",
"Chlamydera_nuchalis_orientalis",
"Ptilonorhynchus_violaceus_violaceus"))

plot(ClimacPtilo)

# Replace remaining subspecific names with species names if any

ClimacPtilo$tip.label[which(ClimacPtilo$tip.label == "Chlamydera_guttata_guttata" )] <- "Chlamydera_guttata" 

ClimacPtilo$tip.label[which(ClimacPtilo$tip.label == "Ailuroedus_arfakianus_arfakianus" )] <- "Ailuroedus_arfakianus"

ClimacPtilo$tip.label[which(ClimacPtilo$tip.label == "Ailuroedus_melanotis_melanotis" )] <- "Ailuroedus_melanotis"


plot(ClimacPtilo)

# Check for duplicated tip labels

anyDuplicated(ClimacPtilo$tip.label)



### CALIBRATE with backbone tree date ###

CladeAge <- getMRCAage(PasBackTre, tip=c("Climacteris_erythrops", "Ptilonorhynchus_violaceus"))
#7.45

CALIB <- makeChronosCalib(ClimacPtilo, node="root", age.min=CladeAge)

ClimacPtiloC <- chronos(ClimacPtilo, model="clock", calibration=CALIB)

plot(ClimacPtiloC); axisPhylo()


### GRAFT new tree into bakcbone ###

# Delete clade in Backbone tree
# Because there are only to terminals involved, it's simple.

PasBackTre2 <- drop.tip(PasBackTre, tip=c("Climacteris_erythrops", "Ptilonorhynchus_violaceus"), trim.internal=FALSE)


# Add the new ClimacPtiloC chronogram to the backbone tree

PasBackTre2 <- bind.tree(PasBackTre2, ClimacPtiloC, where=Ntip(PasBackTre2))


plot(PasBackTre2, cex=0.2, show.tip.label=FALSE); axisPhylo()


### SAVE work ###

write.tree(PasBackTre2, file="PasBackTre2.tre")

# PasBackTre2  <- read.tree(file="PasBackTre2.tre")




#### Corvides ####

# in backbone: Ptilorrhoa caerulescens to Monarcha axillaris

### Get source TREE ###

# Corvides tree from McCullough et al 2022
# extract Corvides, remove the outgroup

CorvTree <- read.nexus("corvides.728-simplified-nexus.tre")
CorvTree <- extract.clade(CorvTree, node=getMRCA(CorvTree, tip=c("ptilorrhoa_leucosticta_16514", "Pericrocotus_albifrons")))

plot(CorvTree)


### Manage TAXA and tipl labels ###

# Show tip labels to select subspecies to delete

CorvTree$tip.label

# Can use cat() and BBedit to arrange lists for cut & paste

# Make sure that one sample remains to represent the species

### No subspecies to drop

# Removing numbers from the ends of species name

CorvTree$tip.label[which(CorvTree$tip.label == "corvus_corax_30042" )] <- "Corvus_corax" 
CorvTree$tip.label[which(CorvTree$tip.label == "ptilostomus_afer_b39279" )] <- "Ptilostomus_afer" 
CorvTree$tip.label[which(CorvTree$tip.label == "cyanocorax_violaceus_21655" )] <- "Cyanocorax_violaceus" 
CorvTree$tip.label[which(CorvTree$tip.label == "cyanocitta_cristata_2271" )] <- "Cyanocitta_cristata"  
CorvTree$tip.label[which(CorvTree$tip.label == "gymnorhinus_cyanocephalus_nk165437" )] <- "Gymnorhinus_cyanocephalus" 
CorvTree$tip.label[which(CorvTree$tip.label == "cyanolyca_viridicyanus_17037" )] <- "Cyanolyca_viridicyanus" 
CorvTree$tip.label[which(CorvTree$tip.label == "perisoreus_canadensis_dot15892" )] <- "Perisoreus_canadensis" 
CorvTree$tip.label[which(CorvTree$tip.label == "cissa_hypoleuca_23241" )] <- "Cissa_hypoleuca"
CorvTree$tip.label[which(CorvTree$tip.label == "dendrocitta_cinerascens_17738" )] <- "Dendrocitta_cinerascens"
CorvTree$tip.label[which(CorvTree$tip.label == "pyrrhocorax_pyrrhocorax_28865" )] <- "Pyrrhocorax_pyrrhocorax"
CorvTree$tip.label[which(CorvTree$tip.label == "lanius_excubitor_28984" )] <- "Lanius_excubitor"
CorvTree$tip.label[which(CorvTree$tip.label == "corvinella_melanoleuca_26670" )] <- "Corvinella_melanoleuca"
CorvTree$tip.label[which(CorvTree$tip.label == "platylophus_galericulatus_24459" )] <- "Platylophus_galericulatus"
CorvTree$tip.label[which(CorvTree$tip.label == "melampitta_lugubris_16552" )] <- "Melampitta_lugubris"
CorvTree$tip.label[which(CorvTree$tip.label == "monarcha_takatsukasae_22596" )] <- "Monarcha_takatsukasae"
CorvTree$tip.label[which(CorvTree$tip.label == "symposiachrus_verticalis_27668" )] <- "Symposiachrus_verticalis"
CorvTree$tip.label[which(CorvTree$tip.label == "myiagra_azureocapilla_24306" )] <- "Myiagra_azureocapilla"
CorvTree$tip.label[which(CorvTree$tip.label == "arses_telescophthalmus_12218" )] <- "Arses_telescophthalmus"
CorvTree$tip.label[which(CorvTree$tip.label == "grallina_cyanoleuca_22946" )] <- "Grallina_cyanoleuca"
CorvTree$tip.label[which(CorvTree$tip.label == "terpsiphone_affinis_23463" )] <- "Terpsiphone_affinis"
CorvTree$tip.label[which(CorvTree$tip.label == "hypothymis_azurea_20189" )] <- "Hypothymis_azurea"
CorvTree$tip.label[which(CorvTree$tip.label == "trochocercus_nitens_15689" )] <- "Trochocercus_nitens"
CorvTree$tip.label[which(CorvTree$tip.label == "paradisaea_minor_16148" )] <- "Paradisaea_minor"
CorvTree$tip.label[which(CorvTree$tip.label == "epimachus_meyeri_16421" )] <- "Epimachus_meyeri"
CorvTree$tip.label[which(CorvTree$tip.label == "lophorina_minor_16456" )] <- "Lophorina_minor"
CorvTree$tip.label[which(CorvTree$tip.label == "parotia_lawesii_16428" )] <- "Parotia_lawesii"
CorvTree$tip.label[which(CorvTree$tip.label == "phonygammus_keraudrenii_12256" )] <- "Phonygammus_keraudrenii"
CorvTree$tip.label[which(CorvTree$tip.label == "ifrita_kowaldi_12106" )] <- "Ifrita_kowaldi"
CorvTree$tip.label[which(CorvTree$tip.label == "corcorax_melanorhamphos_10720" )] <- "Corcorax_melanorhamphos"
CorvTree$tip.label[which(CorvTree$tip.label == "struthidea_cinerea_10728" )] <- "Struthidea_cinerea"
CorvTree$tip.label[which(CorvTree$tip.label == "dicrurus_paradiseus_23288" )] <- "Dicrurus_paradiseus"
CorvTree$tip.label[which(CorvTree$tip.label == "dicrurus_aeneus_23352" )] <- "Dicrurus_aeneus"
CorvTree$tip.label[which(CorvTree$tip.label == "rhipidura_javanica_17717" )] <- "Rhipidura_javanica"
CorvTree$tip.label[which(CorvTree$tip.label == "lamprolia_victoriae_24329" )] <- "Lamprolia_victoriae"
CorvTree$tip.label[which(CorvTree$tip.label == "chaetorhynchus_papuensis_16414" )] <- "Chaetorhynchus_papuensis"
CorvTree$tip.label[which(CorvTree$tip.label == "mystacornis_crossleyi_345860" )] <- "Mystacornis_crossleyi"
CorvTree$tip.label[which(CorvTree$tip.label == "hemipus_picatus_23303" )] <- "Hemipus_picatus"
CorvTree$tip.label[which(CorvTree$tip.label == "tephrodornis_virgatus_30886" )] <- "Tephrodornis_virgatus"
CorvTree$tip.label[which(CorvTree$tip.label == "prionops_plumatus_26690" )] <- "Prionops_plumatus"
CorvTree$tip.label[which(CorvTree$tip.label == "philentoma_pyrhoptera_12323" )] <- "Philentoma_pyrhoptera"
CorvTree$tip.label[which(CorvTree$tip.label == "platysteira_castanea_29120" )] <- "Platysteira_castanea"
CorvTree$tip.label[which(CorvTree$tip.label == "batis_senegalensis_15403" )] <- "Batis_senegalensis"
CorvTree$tip.label[which(CorvTree$tip.label == "laniarius_barbarus_15451" )] <- "Laniarius_barbarus"
CorvTree$tip.label[which(CorvTree$tip.label == "chlorophoneus_sulfureopectus_15498" )] <- "Chlorophoneus_sulfureopectus"
CorvTree$tip.label[which(CorvTree$tip.label == "dryoscopus_cubla_26685" )] <- "Dryoscopus_cubla"
CorvTree$tip.label[which(CorvTree$tip.label == "tchagra_senegalus_15518" )] <- "Tchagra_senegalus"
CorvTree$tip.label[which(CorvTree$tip.label == "aegithina_lafresnayei_23213" )] <- "Aegithina_lafresnayei"
CorvTree$tip.label[which(CorvTree$tip.label == "rhagologus_leucostigma_18382" )] <- "Rhagologus_leucostigma"
CorvTree$tip.label[which(CorvTree$tip.label == "melloria_quoyi_4831" )] <- "Melloria_quoyi"
CorvTree$tip.label[which(CorvTree$tip.label == "strepera_graculina_9660" )] <- "Strepera_graculina"
CorvTree$tip.label[which(CorvTree$tip.label == "artamus_cinereus_6183" )] <- "Artamus_cinereus"
CorvTree$tip.label[which(CorvTree$tip.label == "machaerirhynchus_nigripectus_4734" )] <- "Machaerirhynchus_nigripectus"
CorvTree$tip.label[which(CorvTree$tip.label == "pachycephala_vanikorensis_19410" )] <- "Pachycephala_vanikorensis"
CorvTree$tip.label[which(CorvTree$tip.label == "melanorectes_nigrescens_18387" )] <- "Melanorectes_nigrescens"
CorvTree$tip.label[which(CorvTree$tip.label == "colluricincla_harmonica_8866" )] <- "Colluricincla_harmonica"
CorvTree$tip.label[which(CorvTree$tip.label == "pseudorectes_ferrugineus_9610" )] <- "Pseudorectes_ferrugineus"
CorvTree$tip.label[which(CorvTree$tip.label == "oriolus_chinensis_10450" )] <- "Oriolus_chinensis"
CorvTree$tip.label[which(CorvTree$tip.label == "sphecotheres_vieilloti_10752" )] <- "Sphecotheres_vieilloti"
CorvTree$tip.label[which(CorvTree$tip.label == "pitohui_dichrous_12200" )] <- "Pitohui_dichrous"
CorvTree$tip.label[which(CorvTree$tip.label == "vireo_solitarius_25186" )] <- "Vireo_solitarius"
CorvTree$tip.label[which(CorvTree$tip.label == "erpornis_zantholeuca_27942" )] <- "Erpornis_zantholeuca"
CorvTree$tip.label[which(CorvTree$tip.label == "pteruthius_flaviscapis_17806" )] <- "Pteruthius_flaviscapis"
CorvTree$tip.label[which(CorvTree$tip.label == "oreocharis_arfaki_16440" )] <- "Oreocharis_arfaki"
CorvTree$tip.label[which(CorvTree$tip.label == "aleadryas_rufinucha_16569" )] <- "Aleadryas_rufinucha"
CorvTree$tip.label[which(CorvTree$tip.label == "ornorectes_cristatus_4697" )] <- "Ornorectes_cristatus"
CorvTree$tip.label[which(CorvTree$tip.label == "oreoica_gutturalis_8884" )] <- "Oreoica_gutturalis"
CorvTree$tip.label[which(CorvTree$tip.label == "falcunculus_frontatus_57627" )] <- "Falcunculus_frontatus"
CorvTree$tip.label[which(CorvTree$tip.label == "eulacestoma_nigropectus_16397" )] <- "Eulacestoma_nigropectus"
CorvTree$tip.label[which(CorvTree$tip.label == "psophodes_occidentalis_6204" )] <- "Psophodes_occidentalis"
CorvTree$tip.label[which(CorvTree$tip.label == "mohoua_albicilla_96717" )] <- "Mohoua_albicilla"
CorvTree$tip.label[which(CorvTree$tip.label == "daphoenositta_chrysoptera_23086" )] <- "Daphoenositta_chrysoptera"
CorvTree$tip.label[which(CorvTree$tip.label == "edolisoma_tenuirostre_23644" )] <- "Edolisoma_tenuirostre"
CorvTree$tip.label[which(CorvTree$tip.label == "lalage_maculosa_26428" )] <- "Lalage_maculosa"
CorvTree$tip.label[which(CorvTree$tip.label == "coracina_papuensis_27758" )] <- "Coracina_papuensis"
CorvTree$tip.label[which(CorvTree$tip.label == "pericrocotus_divaricatus_10261" )] <- "Pericrocotus_divaricatus"
CorvTree$tip.label[which(CorvTree$tip.label == "cinclosoma_punctatum_57790" )] <- "Cinclosoma_punctatum"
CorvTree$tip.label[which(CorvTree$tip.label == "ptilorrhoa_leucosticta_16514" )] <- "Ptilorrhoa_leucosticta"

plot(CorvTree)

# Check for duplicated tip labels

anyDuplicated(CorvTree$tip.label)



### CALIBRATE with backbone tree date ###

CladeAge <- getMRCAage(PasBackTre, tip=c("Ptilorrhoa_caerulescens", "Monarcha_axillaris"))
#27.78264

CALIB <- makeChronosCalib(CorvTree, node="root", age.min=CladeAge)

CorvTreeC <- chronos(CorvTree, model="clock", calibration=CALIB)

plot(CorvTreeC); axisPhylo()


### GRAFT new tree into bakcbone ###

# Delete clade in Backbone tree
# Isolate the clade and delete all tips and internal nodes, except 2 species (to get the MRCA node)
CladeToDelete <- extract.clade(PasBackTre2, node=getMRCA(PasBackTre2, tip=c("Ptilorrhoa_caerulescens", "Monarcha_axillaris")))
TipsToDelete <- CladeToDelete$tip.label[CladeToDelete$tip.label!=c("Ptilorrhoa_caerulescens", "Monarcha_axillaris")]
PasBackTre3 <- drop.tip(PasBackTre2, tip=TipsToDelete, trim.internal=T)
plot(PasBackTre3, cex=0.3)

# Delete the 2 species, keeping their internal node
PasBackTre3 <- drop.tip(PasBackTre3, tip=c("Ptilorrhoa_caerulescens", "Monarcha_axillaris"), trim.internal=FALSE)

 
# Add the new Ccorvides chronogram to the backbone tree

PasBackTre3 <- bind.tree(PasBackTre3, CorvTreeC, where=Ntip(PasBackTre3))


plot(PasBackTre3, cex=0.2, show.tip.label=FALSE); axisPhylo()

# PasBackTre3  <- read.tree(file="PasBackTre3.tre")


### SAVE work ###

write.tree(PasBackTre3, file="PasBackTre3.tre")

# PasBackTre2  <- read.tree(file="PasBackTre2.tre")




#### Babblers of Sylvioidea ####

### Get source TREE ###

# Cai et al 2019 Near-complete phylogeny and taxonomic revision of the world’s babblers (Sylvioidea tree)
# babblers_phylogeny.tre
# sylvia to pellorneum, NO bernieria

Babblers <- read.nexus("babblers_phylogeny.tre")
Babblers <- extract.clade(Babblers, node=getMRCA(Babblers, tip=c("Sylvia_atricapilla", "Pellorneum_ruficeps")))

plot(Babblers)


### Manage TAXA and tipl labels ###

# Show tip labels to select subspecies to delete

Babblers$tip.label

# Can use cat() and BBedit to arrange lists for cut & paste

# Make sure that one sample remains to represent the species

## No Subspecies To Drop ##
Babblers$tip.label[(lengths(strsplit(Babblers$tip.label, split = "_")))!=2] #code to find names to change


# Check for duplicated tip labels

anyDuplicated(Babblers$tip.label)



### CALIBRATE with backbone tree date ###

CladeAge <- getMRCAage(PasBackTre, tip=c("Sylvia_nana", "Pellorneum_ruficeps"))
#22.6563

CALIB <- makeChronosCalib(Babblers, node="root", age.min=CladeAge)

BabblersC <- chronos(Babblers, model="clock", calibration=CALIB)

plot(BabblersC); axisPhylo()


### GRAFT new tree into backbone ###

# Delete clade in Backbone tree

# Delete 3 species from the clade of 5, removing internal nodes, leaving 2 for next step
PasBackTre4 <- drop.tip(PasBackTre3, tip=c("Zosterops_senegalensis", "Garrulax_milleti", "Timalia_pileata"), trim.internal=T)

# Delete the 2 species, keeping their internal node
PasBackTre4 <- drop.tip(PasBackTre4, tip=c("Sylvia_nana", "Pellorneum_ruficeps"), trim.internal=FALSE)

# Add the new Babblers chronogram to the backbone tree

PasBackTre4 <- bind.tree(PasBackTre4, BabblersC, where=Ntip(PasBackTre4))


plot(PasBackTre4, cex=0.2, show.tip.label=FALSE); axisPhylo()



### SAVE work ###

write.tree(PasBackTre4, file="PasBackTre4.tre")

# PasBackTre4  <- read.tree(file="PasBackTre4.tre")






#### Emberizoidea ####

#Plectophenax to Setophaga in backbone

### Get source TREE ###
# Barker et al. 2015 

## UPDATE: Santiago edits, removed Oreomystis which was in there for some reason??

Ember <- read.nexus("sptree_MCC.tre")

# Delete outgroup
Ember <- extract.clade(Ember, node=getMRCA(Ember, tip=c("Rhodinocichla_rosea","Emberiza_cia")))

plot(Ember, cex=0.5)
EmberSpp <- as.data.frame((Ember$tip.label))

### Manage TAXA and tipl labels ###

# Show tip labels to select subspecies to delete

Ember$tip.label

## No Subspecies To Drop ##
Ember$tip.label[(lengths(strsplit(Ember$tip.label, split = "_")))!=2] #code to find names to change

# Check for duplicated tip labels

anyDuplicated(Ember$tip.label)

### Determine scaling factor ###

#Obtain backbone age
#CladeAge <- getMRCAage(PasBackTre, tip=c("Plectophenax_nivalis", "Setophaga_americana"))
CladeAge <- getMRCAage(PasBackTre35, tip=c("Oreomystis_bairdii", "Emberiza_rutila"))
#17.15195

#Obtain tree age
TreeAge <- getMRCAage(Ember, tip=c("Rhodinocichla_rosea","Emberiza_cia"))

#Convert branch lengths
Ember$edge.length <- Ember$edge.length*CladeAge/TreeAge


plot(Ember, show.node.label = T); axisPhylo()


### GRAFT new tree into backbone ###

# Delete clade in Backbone tree
# Isolate the clade and delete all tips and internal nodes, except 2 species (to get the MRCA node)
# CladeToDelete <- extract.clade(PasBackTre4, node=getMRCA(PasBackTre4, tip=c("Plectophenax_nivalis", "Setophaga_americana")))
# keep <- c("Plectophenax_nivalis", "Setophaga_americana")
# TipsToDelete <- CladeToDelete$tip.label[!(CladeToDelete$tip.label) %in% keep]

CladeToDelete <- extract.clade(PasBackTre35, node=getMRCA(PasBackTre35, tip=c("Oreomystis_bairdii", "Emberiza_rutila")))
keep <- c("Oreomystis_bairdii", "Emberiza_rutila")
TipsToDelete <- CladeToDelete$tip.label[!(CladeToDelete$tip.label) %in% keep]


#PasBackTre5 <- drop.tip(PasBackTre4, tip=TipsToDelete, trim.internal=T)
PasBackTre36 <- drop.tip(PasBackTre35, tip=TipsToDelete, trim.internal=T)

plot(PasBackTre5, cex=0.3)

# Delete the 2 species, keeping their internal node
#PasBackTre5 <- drop.tip(PasBackTre5, tip=c("Plectophenax_nivalis", "Setophaga_americana"), trim.internal=FALSE)
PasBackTre36 <- drop.tip(PasBackTre36, tip=c("Oreomystis_bairdii", "Emberiza_rutila"), trim.internal=FALSE)


# Add the new Ember tree to the backbone tree

# PasBackTre5 <- bind.tree(PasBackTre5, Ember, where=Ntip(PasBackTre5))
PasBackTre36 <- bind.tree(PasBackTre36, Ember, where=Ntip(PasBackTre36))


plot(PasBackTre36, show.tip.label=FALSE); axisPhylo()


### SAVE work ###

#write.tree(PasBackTre5, file="PasBackTre5.tre")
write.tree(PasBackTre36, file="PasBackTre36.tre")







#### Muscicapidae ####
Musc <- read.tree("musc_newick")

# root tree, idk did this work?
#Musc <- root(Musc, outgroup=c("Catharus_fuscescens", "Sialia_sialis", "Turdus_rufiventris"))

Musc <- midpoint(Musc)

# remove outgroup
#Musc <- extract.clade(Musc, node=getMRCA(Musc, tip=c("Alethe_diademata", "Ficedula_hyperythra")))

plot(Musc, cex=0.3)


# Show tip labels to select subspecies to delete

Musc$tip.label

## No Subspecies To Drop ##
Musc$tip.label[(lengths(strsplit(Musc$tip.label, split = "_")))!=2] #code to find names to change

# Check for duplicated tip labels

anyDuplicated(Musc$tip.label)

### CALIBRATE with backbone tree date ###

# !!! it is only Muscicapa ferruginea

CladeAge <- getMRCAage(PasBackTre, tip=c("Muscicapa_ferruginea", "Catharus_ustulatus"))
#25.97041

CALIB <- makeChronosCalib(Musc, node="root", age.min=CladeAge)

MuscC <- chronos(Musc, model="clock", calibration=CALIB)

plot(MuscC, cex=0.3); axisPhylo()

#MuscC<- extract.clade(MuscC, node=getMRCA(MuscC, tip=c("Alethe_diademata", "Ficedula_hyperythra")))

plot(MuscC, cex=0.3); axisPhylo()


# Get the Muscicapidae root branch length
# Which is just the difference between the root age (CladeAGE) and the crown age

CrownAge <- getMRCAage(MuscC, tip=c("Muscicapa_striata", "Saxicola_ferreus"))

MuscC2 <- drop.tip(MuscC, tip=c("Catharus_fuscescens", "Sialia_sialis", "Turdus_rufiventris"))


MuscC2$root.edge <- CladeAge - CrownAge

plot(MuscC2, cex=0.3, root.edge=TRUE); axisPhylo()

# Delete Muscicapa

PasBackTre7 <- drop.tip(PasBackTre6, tip=c("Nectarinia_olivacea", "Muscicapa_ferruginea"))

#write.tree(PasBackTre7, file="PasBackTre7.tre")

# Bind tree but specifying the position of the new node
PasBackTre7 <- bind.tree(PasBackTre7, MuscC2, where=which(PasBackTre7$tip.label =="Catharus_ustulatus") , position=CladeAge)

### GRAFT new tree into backbone ###



#PasBackTre7 <- graft.tree(phy1=PasBackTre7, phy2=MuscC, edge=c("Muscicapa_ferruginea", "Catharus_ustulatus"), age=CladeAge)
#test <- graft.tree(phy1=PasBackTre, phy2=MuscC, edge=c("Muscicapa_ferruginea", "Catharus_ustulatus"), age=CladeAge)
#PasBackTre7 <- drop.tip(PasBackTre7, tip=2071)
anyDuplicated(PasBackTre7$tip.label)
#temp <- which(PasBackTre7$tip.label %in% c("Muscicapa_ferruginea"))
#which(PasBackTre7$tip.label %in% c("Muscicapa_muttui"))


### SAVE work ###

write.tree(PasBackTre7, file="PasBackTre7.tre")



### Meliphagides ####
# originally MCC_NoOutgroups_0612.tre
# Marki et al 2020
# Data from: Adaptive radiation and the evolution of nectarivory in a large songbird clade
# in backbone: Pardalotus_punctatus to Malurus_melanocephalus

Meli <- read.tree("MCC_NoOutgroups_0612_meli.tre")
#plot(Meli, cex=0.3)
#length(Meli$tip.label)


### Manage TAXA and tipl labels ###
## No Subspecies To Drop ##
Meli$tip.label[(lengths(strsplit(Meli$tip.label, split = "_")))!=2] #code to find names to change

# Check for duplicated tip labels

anyDuplicated(Meli$tip.label)


### CALIBRATE with backbone tree date ###

CladeAge <- getMRCAage(PasBackTre, tip=c("Pardalotus_punctatus", "Malurus_melanocephalus"))
#28.0921

CALIB <- makeChronosCalib(Meli, node="root", age.min=CladeAge)

MeliC <- chronos(Meli, model="clock", calibration=CALIB)

plot(MeliC); axisPhylo()


### GRAFT new tree into backbone ###
PasBackTre8 <- graft.tree(phy1=PasBackTre7, phy2=MeliC, edge=c("Pardalotus_punctatus", "Malurus_melanocephalus"),age=CladeAge)

plot(PasBackTre8, cex=0.2, show.tip.label=FALSE); axisPhylo()
plot(PasBackTre7, cex=0.2, show.tip.label=FALSE); axisPhylo()

# duplicated taxa to get rid of
anyDuplicated(PasBackTre8$tip.label)
temp <- which(PasBackTre8$tip.label %in% c("Malurus_melanocephalus"))
temp <- which(PasBackTre8$tip.label %in% c("Meliphaga_analoga"))
temp <- which(PasBackTre8$tip.label %in% c("Dasyornis_broadbenti"))
temp <- which(PasBackTre8$tip.label %in% c("Acanthiza_pusilla"))
temp <- which(PasBackTre8$tip.label %in% c("Pardalotus_punctatus"))

PasBackTre8 <- drop.tip(PasBackTre8, tip=c(1326, 1327, 1328, 1329, 1330))


### SAVE work ###

write.tree(PasBackTre8, file="PasBackTre8.tre")






###### Petroicidae #######
Petroic <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, c("Sitta_carolinensis", "Petroica_archboldi"))) # the whole big clade, ie with outgroup
#Tregellasia_leucops to extract JUST Petroica

PetroicOnly <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, c("Tregellasia_leucops", "Petroica_archboldi"))) # the whole big clade, ie with outgroup

tips.to.drop <- setdiff(Petroic$tip.label, PetroicOnly$tip.label)
tips.to.drop <- tips.to.drop[! tips.to.drop %in% c("Sitta_carolinensis")]

Petroic <- drop.tip(Petroic, tip=tips.to.drop)

# Petroic <- root(Petroic, outgroup="Sitta_carolinensis") not needed?

#Petroic <- midpoint(Petroic) # midpoint root

plot(Petroic, cex=0.5)


# Show tip labels to select subspecies to delete

Petroic$tip.label

# multicolor was nested in pusilla so i took the one outside
# Heteromyias_albispecularis was formerly in Poecilodryas so i got rid of those ones
# Pachycephalopsis_hattamensis was outside of its group? so i kept the other
Petroic <- drop.tip(Petroic, tip=c(
  "Amalocichla_incerta_brevicauda",
  "Amalocichla_sclateriana_sclateriana",
  "Pachycephalopsis_poliosoma_poliosoma",
  "Petroica_pusilla_polymorpha",
  "Petroica_pusilla_kulambangrae",
  "Petroica_pusilla_dennisi",
  "Petroica_pusilla_septentrionalis",
  "Petroica_pusilla_feminina",
  "Petroica_pusilla_cognata",
  "Petroica_pusilla_similis",  
  "Petroica_pusilla_ambrynensis",
  "Petroica_pusilla_taveunensis",
  "Petroica_pusilla_soror",
  "Petroica_pusilla_kleinschmidti",
  "Petroica_pusilla_becki",
  "Petroica_multicolor",
  "Petroica_rodinogaster_rodinogaster",
  "Petroica_rodinogaster_inexpectata",
  "Petroica_boodang_leggi",
  "Petroica_boodang_boodang",
  "Petroica_boodang_campbelli",
  "Petroica_macrocephala_chathamensis",
  "Petroica_macrocephala_toitoi",
  "Petroica_australis_australis",
  "Peneothello_sigillata_sigillata",
  "Peneothello_cyanus_subcyanea",
  "Peneothello_pulverulenta_alligator",
  "Melanodryas_cucullata_westralensis",
  "Melanodryas_cucullata_cucullata",
  "Melanodryas_cucullata_picata",
  "Melanodryas_vittata_vittata",
  "Peneothello_bimaculata_bimaculata",
  "Eopsaltria_australis_australis",
  "Eopsaltria_australis_chrysorrhos",
  "Eopsaltria_griseogularis_rosinae",      
  "Eopsaltria_griseogularis_griseogularis",
  "Tregellasia_leucops_nigriceps",
  "Tregellasia_leucops_albigularis",       
  "Tregellasia_leucops_albifacies",
  "Tregellasia_capito_capito",
  "Tregellasia_capito_nana",
  "Poecilodryas_hypoleuca_hypoleuca",
  "Poecilodryas_albispecularis_centralis",
  "Poecilodryas_albispecularis_armiti",
  "Drymodes_superciliaris_brevirostris",
  "Microeca_flavovirescens_cuicui",
  "Microeca_griseoceps_kempi",
  "Monachella_muelleriana_muelleriana",
  "Microeca_fascinans_assimilis",
  "Microeca_fascinans_pallida",
  "Microeca_fascinans_fascinans",
  "Microeca_flavigaster_flavigaster",
  "Microeca_flavigaster_laetissima",
  "Microeca_flavigaster_flavissima",
  "Eugerygone_rubra_saturatior",
  "Pachycephalopsis_hattamensis"))

## rename species to remove subspecies name
Petroic$tip.label[(lengths(strsplit(Petroic$tip.label, split = "_")))!=2] #code to find names to change

# replace name with removed subspecies
# IT WORKS OMG
for (tip in Petroic$tip.label) {
  if ((lengths(strsplit(tip, split="_"))>2)) {
    Petroic$tip.label[Petroic$tip.label==tip] <- sub("_[^_]+$", "", tip)
  }
}


# Check for duplicated tip labels

anyDuplicated(Petroic$tip.label)


### CALIBRATE with backbone tree date ###

CladeAge <- getMRCAage(PasBackTre, tip=c("Petroica_cucullata", "Regulus_calendula")) # Petroica is sister to a huge clade
#34.40419 where Petroica splits from the huge clade

CALIB <- makeChronosCalib(Petroic, node="root", age.min=CladeAge)

PetroicC <- chronos(Petroic, model="clock", calibration=CALIB)

plot(PetroicC, cex=0.5); axisPhylo()


# Get the Muscicapidae root branch length
# Which is just the difference between the root age (CladeAGE) and the crown age

CrownAge <- getMRCAage(PetroicC, tip=c("Tregellasia_leucops", "Petroica_archboldi"))

PetroicC <- drop.tip(PetroicC, tip="Sitta_carolinensis")

#MuscC2$root.edge <- CladeAge - CrownAge

plot(MuscC2, cex=0.3, root.edge=TRUE); axisPhylo()

# Delete Petroica

#PasBackTre9 <- drop.tip(PasBackTre8, tip=c("Petroica_cucullata"))

#write.tree(PasBackTre7, file="PasBackTre7.tre")

# Bind tree but specifying the position of the new node
PasBackTre9 <- bind.tree(PasBackTre8, PetroicC, where=which(PasBackTre8$tip.label =="Petroica_cucullata"), position=CrownAge)

PasBackTre9 <- drop.tip(PasBackTre9, tip="Petroica_cucullata")

anyDuplicated(PasBackTre9$tip.label)



### SAVE work ###

write.tree(PasBackTre9, file="PasBackTre9.tre")




###### Turdidae #######
Turdidae <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, c("Turdus_plebejus", "Muscicapa_striata"))) # Turdidae and sister Muscicapidae
#Tregellasia_leucops to extract JUST Petroica
#plot(Turdidae, cex=0.2)


TurdidaeOnly <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, c("Turdus_plebejus", "Grandala_coelicolor"))) # Only Turdidae
#plot(TurdidaeOnly, cex=0.4)


tips.to.drop <- setdiff(Turdidae$tip.label, TurdidaeOnly$tip.label)
tips.to.drop <- tips.to.drop[! tips.to.drop %in% c("Muscicapa_striata")]

Turdidae <- drop.tip(Turdidae, tip=tips.to.drop)
#plot(Turdidae, cex=0.4)



# Show tip labels to select subspecies to delete

Turdidae$tip.label

# remove subspecies
Turdidae <- drop.tip(Turdidae, tip=c(
  "Myadestes_coloratus_coloratus",
  "Psophocichla_litsitsirupa_litsitsirupa",
  "Psophocichla_litsitsirupa_simensis",
  "Geokichla_citrina_innotata",
  "Geokichla_sibirica_davisoni"))


## no subspecies to rename
Turdidae$tip.label[(lengths(strsplit(Turdidae$tip.label, split = "_")))!=2] #code to find names to change


# Check for duplicated tip labels
anyDuplicated(Turdidae$tip.label)


### CALIBRATE with backbone tree date ###

CladeAge <- getMRCAage(PasBackTre, tip=c("Catharus_ustulatus", "Muscicapa_ferruginea")) # Petroica is sister to a huge clade
#25.97041 is MRCA of Turdidae and Muscicapidae

CALIB <- makeChronosCalib(Turdidae, node="root", age.min=CladeAge)

TurdidaeC <- chronos(Turdidae, model="clock", calibration=CALIB)

plot(TurdidaeC, cex=0.5); axisPhylo()


# get age of crown ie only Turdidae
CrownAge <- getMRCAage(TurdidaeC, tip=c("Turdus_plebejus", "Grandala_coelicolor"))

# get rid of outgroup Muscicapa
TurdidaeC <- drop.tip(TurdidaeC, tip="Muscicapa_striata")

plot(TurdidaeC, cex=0.3, root.edge=TRUE); axisPhylo()


# Bind tree but specifying the position of the new node
PasBackTre10 <- bind.tree(PasBackTre9, TurdidaeC, where=which(PasBackTre9$tip.label =="Catharus_ustulatus"), position=CrownAge)
plot(PasBackTre10, cex=0.1)

anyDuplicated(PasBackTre10$tip.label)

PasBackTre10 <- drop.tip(PasBackTre10, tip=2394)



### SAVE work ###

write.tree(PasBackTre10, file="PasBackTre10.tre")




###### Sturnidae, Buphagidae, and Mimidae #######
# backbone tree sais Buphagidae sister to Sturnidae but probably actually the outgroup
SBM <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, c("Sturnus_vulgaris", "Cinclus_cinclus"))) # The 3 fams plus representative of sister clade
#Mimus_triurus to extract JUST SBM
plot(SBM, cex=0.2)


SBMOnly <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, c("Sturnus_vulgaris", "Buphagus_erythrorhynchus"))) # Only SBM
plot(SBMOnly, cex=0.4)


tips.to.drop <- setdiff(SBM$tip.label, SBMOnly$tip.label)
tips.to.drop <- tips.to.drop[! tips.to.drop %in% c("Cinclus_cinclus")] # keep representative of outgroup

SBM <- drop.tip(SBM, tip=tips.to.drop)
plot(SBM, cex=0.4)



# Show tip labels to select subspecies to delete

SBM$tip.label

# remove subspecies
SBM <- drop.tip(SBM, tip=c("Sarcops_calvus_melanonotus", "Sarcops_calvus_calvus"))


# find subspecies to rename
SBM$tip.label[(lengths(strsplit(SBM$tip.label, split = "_")))!=2] #code to find names to change


# Check for duplicated tip labels
anyDuplicated(SBM$tip.label)


### CALIBRATE with backbone tree date ###

CladeAge <- getMRCAage(PasBackTre, tip=c("Sturnus_vulgaris", "Cinclus_cinclus")) # node between SBM and sister clade
#28.38762

CALIB <- makeChronosCalib(SBM, node="root", age.min=CladeAge)

SBMC <- chronos(SBM, model="clock", calibration=CALIB)

plot(SBMC, cex=0.5); axisPhylo()


# get age of crown ie only SBM
CrownAge <- getMRCAage(SBMC, tip=c("Sturnus_vulgaris", "Buphagus_erythrorhynchus")) 

# get rid of outgroup 
SBMC <- drop.tip(SBMC, tip="Cinclus_cinclus")
plot(SBMC, cex=0.3, root.edge=TRUE); axisPhylo()

# Remove spp from backbone so clade of 3 becomes 1 branch
PasBackTre11 <- drop.tip(PasBackTre10, tip=c("Sturnus_vulgaris", "Buphagus_erythrorhynchus"))

# Bind tree but specifying the position of the new node
PasBackTre11 <- bind.tree(PasBackTre11, SBMC, where=which(PasBackTre10$tip.label =="Mimus_patagonicus"), position=CrownAge)
plot(PasBackTre11, cex=0.1)

# check for duplicate
anyDuplicated(PasBackTre11$tip.label)

# remove duplicate
PasBackTre11 <- drop.tip(PasBackTre11, tip=2390)



### SAVE work ###

write.tree(PasBackTre11, file="PasBackTre11.tre")



####### Locustelidae +  Donacobiidae + Bernieriidae #######
Locust <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, c("Megalurus_palustris", "Acrocephalus_newtoni")))
plot(Locust, cex=0.4)


LocustOnly <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, c("Megalurus_palustris", "Bernieria_madagascariensis"))) # Only SBM
plot(LocustOnly, cex=0.5)


tips.to.drop <- setdiff(Locust$tip.label, LocustOnly$tip.label)
tips.to.drop <- tips.to.drop[! tips.to.drop %in% c("Acrocephalus_newtoni")] # keep representative of outgroup

Locust <- drop.tip(Locust, tip=tips.to.drop)
plot(Locust, cex=0.5)



# Show tip labels to select subspecies to delete

Locust$tip.label

# remove subspecies
Locust <- drop.tip(Locust, tip=c("Donacobius_atricapilla_atricapilla",
                                 "Locustella_fasciolata_fasciolata", #outside fasciolata and amnicola?
                                 "Locustella_pryeri_pryeri", # pryeri outside this and Acrocephalus palustris sinensis?
                                 "Locustella_naevia_naevia",
                                 "Locustella_certhiola", # again, weird
                                "Locustella_luscinioides_luscinioides",
                                "Locustella_castanea_castanea", #outside a group that had castanea?
                                "Bradypterus_caudatus_unicolor",
                                "Locustella_mandelli_mandelli",
                                "Bradypterus_davidi_davidi",
                                "Bradypterus_davidi_suschkini",
                                "Bradypterus_baboecala_elgonensis", # this and centralis clustered together so should i use?
                                "Bradypterus_baboecala_centralis",
                                "Bradypterus_baboecala_tongensis/transvaalensis", # outside a group of baboecala
                                "Bradypterus_cinnamomeus_cinnamomeus",
                                "Cincloramphus_timoriensis_tweeddalei", # actually no idea 
                                "Cincloramphus_timoriensis_crex", # also no idea, so weird
                                "Poodytes_gramineus", #CHECK, either i keep Poodytes together or trust the ones without subspecies
                                "Poodytes_punctatus_vealeae",
                                "Megalurus_palustris_toklao",
                                "Megalurus_palustris_forbesi",
                                "Poodytes_punctatus", # again, CHECK
                                "Bernieria_madagascariensis_inceleber",
                                "Bernieria_madagascariensis_madagascariensis"))


# find subspecies to rename
Locust$tip.label[(lengths(strsplit(Locust$tip.label, split = "_")))!=2] #code to find names to change


# replace name with removed subspecies
# IT WORKS OMG
for (tip in Locust$tip.label) {
  if ((lengths(strsplit(tip, split="_"))>2)) {
    Locust$tip.label[Locust$tip.label==tip] <- sub("_[^_]+$", "", tip)
  }
}

# Check for duplicated tip labels
anyDuplicated(Locust$tip.label)


### CALIBRATE with backbone tree date ###

CladeAge <- getMRCAage(PasBackTre, tip=c("Megalurus_palustris", "Acrocephalus_newtoni")) # node between Locustelidae+ and sister clade
#27.60971

CALIB <- makeChronosCalib(Locust, node="root", age.min=CladeAge)

LocustC <- chronos(Locust, model="clock", calibration=CALIB)

plot(LocustC, cex=0.5); axisPhylo()


# get age of crown ie only Locust+
CrownAge <- getMRCAage(LocustC, tip=c("Megalurus_palustris", "Bernieria_madagascariensis")) 

# get rid of outgroup 
LocustC <- drop.tip(LocustC, tip="Acrocephalus_newtoni")
plot(LocustC, cex=0.3, root.edge=TRUE); axisPhylo()

# Remove spp from backbone so clade of 3 becomes 1 branch
PasBackTre12 <- drop.tip(PasBackTre11, tip=c("Donacobius_atricapillus", "Megalurus_palustris")) # i don't know why atricapillus not atricapilla?

# Bind tree but specifying the position of the new node
PasBackTre12 <- bind.tree(PasBackTre12, LocustC, where=which(PasBackTre12$tip.label =="Bernieria_madagascariensis"), position=CrownAge)
plot(PasBackTre12, cex=0.1)

# check for duplicate
anyDuplicated(PasBackTre12$tip.label)

# remove duplicate
PasBackTre12 <- drop.tip(PasBackTre12, tip=3766)


### SAVE work ###

write.tree(PasBackTre12, file="PasBackTre12.tre")



######## Pycnonotidae ########
# Shakya and Sheldon 2017
# genera polyphyletic/paraphyletic at times?
 Pycno <- read.tree("pycno.tre")
 Pycno <- Pycno[[2]]
plot(Pycno, cex=0.5)

# fix names
Pycno$tip.label <- gsub("'", '', Pycno$tip.label)
Pycno$tip.label <- gsub(" ", '_', Pycno$tip.label)


# Show tip labels to select subspecies to delete
Pycno$tip.label

# no subspp to remove

# find subspecies to rename
Pycno$tip.label[(lengths(strsplit(Pycno$tip.label, split = "_")))!=2] # no names to change

# Check for duplicated tip labels
anyDuplicated(Pycno$tip.label)

# bayesian tree but not all tips aligned on right - is it calibrated??
# ### Determine scaling factor ###
# # 
# # Obtain backbone age
# CladeAge <- getMRCAage(PasBackTre, tip=c("Pycnonotus_barbatus", "Sylvia_nana")) # clade origin
# #25.10564
# #
# #Obtain tree age
# TreeAge <- getMRCAage(Pycno, tip=c("Pycnonotus_barbatus", "Andropadus_importunus"))
# 
# #Convert branch lengths
# Pycno$edge.length <- Pycno$edge.length*CladeAge/TreeAge
# 
# 
# plot(Pycno, show.node.label = T); axisPhylo() 
# 

# ### GRAFT new tree into backbone ###
# 
# # Delete clade in Backbone tree
# # Isolate the clade and delete all tips and internal nodes, except 2 species (to get the MRCA node)
# CladeToDelete <- extract.clade(PasBackTre4, node=getMRCA(PasBackTre4, tip=c("Plectophenax_nivalis", "Setophaga_americana")))
# keep <- c("Plectophenax_nivalis", "Setophaga_americana")
# TipsToDelete <- CladeToDelete$tip.label[!(CladeToDelete$tip.label) %in% keep]
# 
# 
# PasBackTre5 <- drop.tip(PasBackTre4, tip=TipsToDelete, trim.internal=T)
# plot(PasBackTre5, cex=0.3)
# 
# # Delete the 2 species, keeping their internal node
# PasBackTre5 <- drop.tip(PasBackTre5, tip=c("Plectophenax_nivalis", "Setophaga_americana"), trim.internal=FALSE)
# 
# 
# # Add the new Ember tree to the backbone tree
# 
# PasBackTre5 <- bind.tree(PasBackTre5, Ember, where=Ntip(PasBackTre5))
# 
# 
# plot(PasBackTre5, cex=0.2, show.tip.label=FALSE); axisPhylo()
# 
# 
# ### SAVE work ###
# 
# write.tree(PasBackTre5, file="PasBackTre5.tre")

### CALIBRATE with backbone tree date ###

CladeAge <- getMRCAage(PasBackTre, tip=c("Pycnonotus_barbatus", "Sylvia_nana")) # node between Pycnonotidae and rest
#25.10564

CALIB <- makeChronosCalib(Pycno, node="root", age.min=CladeAge)

PycnoC <- chronos(Pycno, model="clock", calibration=CALIB)

plot(PycnoC, cex=0.5); axisPhylo()


# get age of crown ie only Pycno
#Pycno <- extract.clade(Pycno, node=getMRCA(Pycno, c("Pycnonotus_barbatus", "Andropadus_importunus"))) # Only Pycnonotidae

CrownAge <- getMRCAage(PycnoC, tip=c("Pycnonotus_barbatus", "Andropadus_importunus"))

# get rid of outgroup 
PycnoC <- drop.tip(PycnoC, tip=c("Acrocephalus_newtoni", "Berniera_madagascariensis", "Nicator_chloris", "Hirundo_rustica", "Malia_grata"))
plot(PycnoC, cex=0.3, root.edge=TRUE); axisPhylo()


# Bind tree but specifying the position of the new node
PasBackTre13 <- bind.tree(PasBackTre12, PycnoC, where=which(PasBackTre12$tip.label =="Pycnonotus_barbatus"), position=CrownAge)
plot(PasBackTre13, cex=0.1)

# check for duplicate
anyDuplicated(PasBackTre13$tip.label)

# remove duplicate
PasBackTre13 <- drop.tip(PasBackTre13, tip=3844)


### SAVE work ###

write.tree(PasBackTre13, file="PasBackTre13.tre")




######## Estrildidae #######
# Olsson 2020
# calibrated so only need to rescale
# sister to vidua
Est <- read.tree("estrildidae.tre")
plot(Est, cex=0.4)
Estspp <- as.data.frame(Est$tip.label)

# fix names
Est$tip.label <- gsub(" ", '_', Est$tip.label)
Est$tip.label <- gsub("'", '', Est$tip.label)
Est$tip.label<- gsub('[[:digit:]]+', '', Est$tip.label)
Est$tip.label <- gsub(" ", '_', Est$tip.label)

Est$tip.label[[69]] <- "Oreostruthus_fuliginosus"

# keep only binomial name
Est$tip.label <- stringr::str_extract(Est$tip.label, "[^_]*_[^_]*") # the NAs are Menura and Oreostruthus (okay to get rid of)

# # extract just estrildidae and viduidae
Est <- extract.clade(Est, node=getMRCA(Est, tip=c("Taenopygia_castanotis", "Vidua_macroura")))
plot(Est, cex=0.5)

# check for duplicates
anyDuplicated(Est$tip.label)

# drop duplicates
Est <- drop.tip(Est, tip=c(26, 36))


### Determine scaling factor ###

#Obtain backbone age
CladeAge <- getMRCAage(PasBackTre, tip=c("Taeniopygia_guttata", "Vidua_macroura"))
#19.44128
# 
#Obtain tree age
TreeAge <- getMRCAage(Est, tip=c("Taenopygia_castanotis","Vidua_macroura"))
# 27.79876
#Convert branch lengths
Est$edge.length <- Est$edge.length*CladeAge/TreeAge

# get Estrildidae age
CrownAge <- getMRCAage(Est, tip=c("Taenopygia_castanotis", "Erythrura_hyperythra"))
# 13.5964


plot(Est, show.node.label = T); axisPhylo()


# remove outgroup 
Est <- drop.tip(Est, tip=c("Vidua_chalybeata", "Vidua_funerea", "Vidua_regia", "Vidua_fischeri", "Vidua_macroura", "Vidua_paradisaea", "Anomalospiza_imberbis"))


# Bind tree but specifying the position of the new node
PasBackTre14 <- bind.tree(PasBackTre13, Est, where=which(PasBackTre13$tip.label =="Taeniopygia_guttata"), position=CrownAge)
plot(PasBackTre14, cex=0.1) # i think it worked!

# check for duplicate
anyDuplicated(PasBackTre14$tip.label)

# remove duplicate
PasBackTre14 <- drop.tip(PasBackTre14, tip="Taeniopygia_guttata") # remove the one used for graftin


### SAVE work ###

write.tree(PasBackTre14, file="PasBackTre14.tre")






##### Fringillidae ######
# extract from passeroidea fringillidae + 1 member of outgroup, extract from backbone the emberizoidea plus fringilla, bind together (4-branch polytomy), remove 2 outgroups
#Fringilla_montifringilla, Emberiza_affinis
Fring <- extract.clade(Pas, node=getMRCA(Pas, c("Fringilla_montifringilla", "Emberiza_affinis"))) # the whole big clade, ie with outgroup
#Tregellasia_leucops to extract JUST Petroica

FringOnly <- extract.clade(Pas, node=getMRCA(Pas, c("Fringilla_montifringilla", "Pinicola_subhimachala"))) # only fringilla

tips.to.drop <- setdiff(Fring$tip.label, FringOnly$tip.label)
tips.to.drop <- tips.to.drop[! tips.to.drop %in% c("Emberiza_affinis")]

Fring <- drop.tip(Fring, tip=tips.to.drop)
plot(Fring, cex=0.4)

# remove tips
# convert Hemignathus_stejnegeri, Carpodacus_rubescens, Carpodacus_erythrinus_roseatus,Pyrrhula_leucogenis_leucogenis,Spinus_magellanicus_magellanicus, 
Fring <- drop.tip(Fring, tip=c("Crithagra_sulphurata_sulphurata",
                                 "Fringilla_teydea_teydea", 
                               "Fringilla_coelebs_coelebs", 
                               "Fringilla_coelebs_palmae", 
                               "Euphonia_affinis_olmecorum", 
                               "Euphonia_affinis_affinis", 
                               "Carpodacus_sibiricus_lepidus", 
                               "Carpodacus_vinaceus_vinaceus", 
                               "Loxops_coccineus_coccineus", 
                               "Hemignathus_virens_virens",
                               "Bucanetes_githagineus_crassirostris", 
                               "Bucanetes_githagineus_amantum", 
                               "Pinicola_enucleator_flammula", 
                               "Pinicola_enucleator_leucura", 
                               "Pinicola_enucleator_carlottae", 
                               "Pinicola_enucleator_kamtschatkensis", 
                               "Pinicola_enucleator_enucleator", 
                               "Pyrrhula_pyrrhula_griseiventris", 
                               "Pyrrhula_pyrrhula_pileata", 
                               "Pyrrhula_pyrrhula_pyrrhula", 
                               "Pyrrhula_pyrrhula_europoea", 
                               "Pyrrhula_pyrrhula_iberiae", 
                               "Pyrrhula_pyrrhula_cineracea", 
                               "Pyrrhula_erythaca_erythaca", 
                               "Pyrrhula_leucogenis_steerei", 
                               "Pyrrhula_nipalensis_nipalensis", 
                               "Pyrrhula_nipalensis_uchidai",
                               "Pyrrhula_nipalensis_ricketti",
                               "Acanthis_flammea_islandica", 
                               "Acanthis_flammea_rostrata",
                               "Acanthis_flammea_flammea",
                               "Acanthis_hornemanni_exilipes",
                               "Acanthis_hornemanni_hornemanni",
                               "Loxia_leucoptera_bifasciata",
                               "Loxia_curvirostra_japonica",
                               "Loxia_curvirostra_guillemardi",
                               "Loxia_curvirostra_curvirostra",
                               "Loxia_curvirostra_balearica",
                               "Loxia_curvirostra_poliogyna",
                               "Loxia_curvirostra_corsicana",
                               "Spinus_pinus_perplexus",
                               "Spinus_pinus_pinus",
                               "Spinus_spinescens_spinescens",
                               "Spinus_notatus_notatus",
                               "Spinus_psaltria_hesperophilus",
                               "Spinus_psaltria_colombianus",
                               "Spinus_tristis_salicamans",
                               "Linaria_cannabina_bella",
                               "Linaria_cannabina_cannabina",
                               "Linaria_cannabina_meadewaldoi",
                               "Linaria_cannabina_harterti",
                               "Linaria_flavirostris_flavirostris",
                               "Linaria_flavirostris_rufostrigata",
                               "Chloris_ambigua_ambigua",
                               "Chloris_chloris_chloris",
                               "Chloris_chloris_aurantiiventris",
                               "Linurgus_olivaceus_olivaceus",
                               "Linurgus_olivaceus_kilimensis",
                               "Serinus_dorsostriatus_dorsostriatus",
                               "Carpodacus_erythrinus_roseatus",
                               "Spinus_magellanicus_magellanicus"))
                               
for (tip in Fring$tip.label) {
  if ((lengths(strsplit(tip, split="_"))>2)) {
    Fring$tip.label[Fring$tip.label==tip] <- sub("_[^_]+$", "", tip)
  }
}

# Check for duplicated tip labels
anyDuplicated(Fring$tip.label)

### CALIBRATE with backbone tree date ###

CladeAge <- getMRCAage(PasBackTre, tip=c("Fringilla_montifringilla", "Emberiza_schoeniclus")) # node between fringillidae and rest
#20.28843

CALIB <- makeChronosCalib(Fring, node="root", age.min=CladeAge)

FringC <- chronos(Fring, model="clock", calibration=CALIB)

plot(FringC, cex=0.5); axisPhylo()

# only crown fring age
CrownAge <- getMRCAage(FringC, tip=c("Fringilla_montifringilla", "Pinicola_subhimachala"))

# get rid of outgroup 
FringC <- drop.tip(FringC, tip=c("Emberiza_affinis"))
FringC <- drop.tip(FringC, tip=c("Carpospiza_brachydactyla")) # in Fringilla for some reason?

plot(FringC, cex=0.3, root.edge=TRUE); axisPhylo()


# Bind tree but specifying the position of the new node
PasBackTre15 <- bind.tree(PasBackTre14, FringC, where=which(PasBackTre14$tip.label =="Fringilla_montifringilla"), position=CrownAge)
plot(PasBackTre15, cex=0.1)

# check for duplicate
anyDuplicated(PasBackTre15$tip.label)

# removing duplicates
PasBackTre15 <- drop.tip(PasBackTre15, tip=3058) # fringilla montifringilla
# idk why these others are in here but as long as you run it in order, should be right
PasBackTre15 <- drop.tip(PasBackTre15, tip=3850) # Chlorophonia_flavirostris
PasBackTre15 <- drop.tip(PasBackTre15, tip=3849) # Euphonia_laniirostris
PasBackTre15 <- drop.tip(PasBackTre15, tip=3850) # Paroreomyza_montana


### SAVE work ###

write.tree(PasBackTre15, file="PasBackTre15.tre")






### Motacillidae + Passeridae#####
#Extract them, remove all others so they are sister like in backbone

MotPass <- extract.clade(Pas, node=getMRCA(Pas, c("Motacilla_cinerea", "Passer_montanus"))) # the whole big clade, both of them
plot(MotPass, cex=0.2)

# just the 2 fams
Mot <- extract.clade(Pas, node=getMRCA(Pas, c("Motacilla_cinerea", "Macronyx_croceus"))) 
plot(Mot)

Pass <- extract.clade(Pas, node=getMRCA(Pas, c("Passer_montanus", "Onychostruthus_taczanowskii"))) 
plot(Pass)

# keep only these 2 fams
# wow genius it worked 
tips.to.drop <- setdiff(MotPass$tip.label, Mot$tip.label)
tips.to.drop <- tips.to.drop[! tips.to.drop %in% Pass$tip.label]

MotPass <- drop.tip(MotPass, tip=tips.to.drop)
plot(MotPass, cex=0.4)

MotPass <- drop.tip(MotPass, tip=c("Petronia_petronia_petronia",
                                   "Montifringilla_nivalis_nivalis",
                                   "Montifringilla_nivalis_groumgrzimaili",
                                   "Passer_domesticus_domesticus",
                                   "Passer_domesticus_bactrianus",
                                   "Passer_domesticus_persicus",
                                   "Passer_domesticus_hyrcanus",
                                   "Passer_domesticus_niloticus",
                                   "Passer_domesticus_biblicus",
                                   "Passer_domesticus_indicus",
                                   "Passer_montanus_montanus",
                                   "Passer_montanus_saturatus",
                                   "Gymnoris_xanthocollis_xanthocollis",
                                   "Anthus_gutturalis_wollastoni",
                                   "Anthus_trivialis_trivialis",
                                   "Anthus_hodgsoni_hodgsoni",
                                   "Anthus_gutturalis_rhododendri",
                                   "Anthus_rubescens_japonicus",
                                   "Anthus_spinoletta_coutellii",
                                   "Anthus_petrosus_littoralis",
                                   "Anthus_spinoletta_spinoletta",
                                   "Anthus_pratensis_pratensis",
                                   "Anthus_gustavi_gustavi",
                                   "Anthus_gustavi_menzbieri",
                                   "Anthus_correndera_chilensis",
                                   "Anthus_cinnamomeus_lichenya",
                                   "Anthus_berthelotii_madeirensis",
                                   "Anthus_berthelotii_berthelotii",
                                   "Motacilla_clara_torrentium",
                                   "Motacilla_capensis_wellsi",
                                   "Motacilla_aguimp_vidua",
                                   "Motacilla_cinerea_patriciae",
                                   "Motacilla_cinerea_canariensis",
                                   "Motacilla_alba_alba",
                                   "Motacilla_alba_ocularis",
                                   "Motacilla_alba_baicalensis",
                                   "Motacilla_alba_personata",
                                   "Motacilla_alba_persica",
                                   "Motacilla_alba_alboides",
                                   "Motacilla_alba_leucopsis",
                                   "Motacilla_citreola_citreola",
                                   "Motacilla_flava_taivana",
                                   "Motacilla_flava_macronyx",
                                   "Motacilla_flava_flava",               
                                   "Motacilla_flava_flavissima",
                                   "Motacilla_flava_thunbergi",          
                                   "Motacilla_flava_pygmaea",
                                   "Motacilla_flava_cinereocapilla" ,      
                                   "Motacilla_flava_lutea",
                                   "Motacilla_flava_feldegg",             
                                   "Motacilla_flava_beema",
                                   "Motacilla_flava_iberiae",              
                                   "Motacilla_citreola_calcarata",
                                  "Motacilla_cinerea_cinerea",
                                  "Passer_melanurus_melanurus",
                                  "Anthus_rufulus_lugubris"))
                                   
for (tip in MotPass$tip.label) {
  if ((lengths(strsplit(tip, split="_"))>2)) {
    MotPass$tip.label[MotPass$tip.label==tip] <- sub("_[^_]+$", "", tip)
  }
}

# Check for duplicated tip labels
anyDuplicated(MotPass$tip.label)                             


### CALIBRATE with backbone tree date ###

CladeAge <- getMRCAage(PasBackTre, tip=c("Motacilla_cinerea", "Passer_montanus")) 
#20.40578

CALIB <- makeChronosCalib(MotPass, node="root", age.min=CladeAge)

MotPassC <- chronos(MotPass, model="clock", calibration=CALIB)

plot(MotPassC, cex=0.5); axisPhylo()

# remove 2 tips from backbone
PasBackTre16 <- drop.tip(PasBackTre15, tip=c("Motacilla_cinerea", "Passer_montanus"),  trim.internal=FALSE)

# Graft onto Backbone
PasBackTre16 <- bind.tree(PasBackTre16, MotPassC, where=Ntip(PasBackTre16))



# check for duplicate
anyDuplicated(PasBackTre16$tip.label)


plot(PasBackTre16, show.tip.label = F)


### SAVE work ###

write.tree(PasBackTre16, file="PasBackTre16.tre")
                               
               

                


###### Prunellidae ######
# extract clade + outgroup, remove all except 1 representative
# all in Prunella, and sister to Peucedramus, monophyletic group, chef's kiss
Prune <- extract.clade(Pas, node=getMRCA(Pas, tip=c("Peucedramus_taeniatus", "Prunella_collaris")))
plot(Prune, cex=0.4)

PruneOnly <- extract.clade(Pas, node=getMRCA(Pas, c("Prunella_collaris", "Prunella_rubeculoides"))) 

# keep only these 2 fams
# wow genius it worked 
tips.to.drop <- setdiff(Prune$tip.label, PruneOnly$tip.label)
tips.to.drop <- tips.to.drop[! tips.to.drop %in% c("Peucedramus_taeniatus")]

Prune <- drop.tip(Prune, tip=tips.to.drop)
plot(Prune, cex=0.6)

# drop subspp
Prune <- drop.tip(Prune, tip=c("Prunella_atrogularis_huttoni",
                               "Prunella_atrogularis_atrogularis",
                               "Prunella_fulvescens_nanschanica",
                               "Prunella_fulvescens_fulvescens",
                               "Prunella_collaris_montana",       
                               "Prunella_collaris_collaris",
                               "Prunella_collaris_erythropygia",
                               "Prunella_collaris_nipalensis" ,   
                               "Prunella_collaris_fennelli"))
                               
for (tip in Prune$tip.label) {
  if ((lengths(strsplit(tip, split="_"))>2)) {
    Prune$tip.label[Prune$tip.label==tip] <- sub("_[^_]+$", "", tip)
  }
}

# Check for duplicated tip labels
anyDuplicated(Prune$tip.label)                             


### CALIBRATE with backbone tree date ###

CladeAge <- getMRCAage(PasBackTre, tip=c("Prunella_collaris", "Peucedramus_taeniatus")) 
#21.63265

CALIB <- makeChronosCalib(Prune, node="root", age.min=CladeAge)

PruneC <- chronos(Prune, model="clock", calibration=CALIB)

plot(PruneC, cex=0.5); axisPhylo()

# remove 2 tips from backbone
PasBackTre17 <- drop.tip(PasBackTre16, tip=c("Prunella_collaris", "Peucedramus_taeniatus"),  trim.internal=FALSE)

# Graft onto Backbone
PasBackTre17 <- bind.tree(PasBackTre17, PruneC, where=Ntip(PasBackTre17))



# check for duplicate
anyDuplicated(PasBackTre17$tip.label)


plot(PasBackTre17, show.tip.label = F)

### SAVE work ###

write.tree(PasBackTre17, file="PasBackTre17.tre")                             
                                   
                                   
                                   


##### Viduidae #####
# sister in passeroidea tree, easy
# then delete taeniopygia, second step delete other 2 but trim.internal=FALSE
Vid <- extract.clade(Pas, node=getMRCA(Pas, tip=c("Vidua_macroura", "Taeniopygia_guttata")))
#Anomalospiza_imberbis
plot(Vid)

VidOnly <- extract.clade(Pas, node=getMRCA(Pas, c("Vidua_macroura", "Anomalospiza_imberbis"))) 

# keep only these 2 fams
# wow genius it worked 
tips.to.drop <- setdiff(Vid$tip.label, VidOnly$tip.label)
tips.to.drop <- tips.to.drop[! tips.to.drop %in% c("Taeniopygia_guttata")]

Vid <- drop.tip(Vid, tip=tips.to.drop)
plot(Vid, cex=0.6)

# drop subspp
Vid <- drop.tip(Vid, tip=c("Anomalospiza_imberbis_butleri",
                               "Anomalospiza_imberbis_imberbis"))
                               
for (tip in Vid$tip.label) {
  if ((lengths(strsplit(tip, split="_"))>2)) {
    Vid$tip.label[Vid$tip.label==tip] <- sub("_[^_]+$", "", tip)
  }
}

# Check for duplicated tip labels
anyDuplicated(Vid$tip.label)                             


### CALIBRATE with backbone tree date ###

CladeAge <- getMRCAage(PasBackTre, tip=c("Vidua_macroura", "Taeniopygia_guttata")) 
#19.44128

CALIB <- makeChronosCalib(Vid, node="root", age.min=CladeAge)

VidC <- chronos(Vid, model="clock", calibration=CALIB)

plot(VidC, cex=0.5); axisPhylo()

## only crown viduidae age
CrownAge <- getMRCAage(VidC, tip=c("Vidua_macroura", "Anomalospiza_imberbis"))

# get rid of outgroup 
VidC <- drop.tip(VidC, tip=c("Taeniopygia_guttata"))

plot(VidC, cex=0.8, root.edge=TRUE); axisPhylo()


# Bind tree but specifying the position of the new node
PasBackTre18 <- bind.tree(PasBackTre17, VidC, where=which(PasBackTre17$tip.label =="Vidua_macroura"), position=CrownAge)
plot(PasBackTre18, show.tip.label = F)

# check for duplicate
anyDuplicated(PasBackTre18$tip.label)

# removing duplicates
PasBackTre18 <- drop.tip(PasBackTre18, tip=3052) 

### SAVE work ###

write.tree(PasBackTre18, file="PasBackTre18.tre")




######## Phylloscopidae #######
# can i use the phlawd tree for phylloscopidae and the others or are there not enough in here?
# ie do this whole clade at once or extract the families
Phyllo <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, c("Phylloscopus_olivaceus", "Hylia_prasina")))
#Hylia_prasina
#Pnoepyga_pusilla
plot(Phyllo, cex=0.4)

PhylloOnly <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, c("Phylloscopus_olivaceus", "Phylloscopus_neglectus")))

#keep only Phyllo and outgroup
tips.to.drop <- setdiff(Phyllo$tip.label, PhylloOnly$tip.label)
tips.to.drop <- tips.to.drop[! tips.to.drop %in% c("Hylia_prasina")]

Phyllo <- drop.tip(Phyllo, tip=tips.to.drop)
plot(Phyllo, cex=0.6)

# drop subspp
Phyllo <- drop.tip(Phyllo, tip=c("Seicercus_affinis_ocularis",
                                "Seicercus_affinis_intermedius"))
                                
for (tip in Phyllo$tip.label) {
  if ((lengths(strsplit(tip, split="_"))>2)) {
    Phyllo$tip.label[Phyllo$tip.label==tip] <- sub("_[^_]+$", "", tip)
  }
}

# Check for duplicated tip labels
anyDuplicated(Phyllo$tip.label)                             


### CALIBRATE with backbone tree date ###

CladeAge <- getMRCAage(PasBackTre, tip=c("Phylloscopus_collybita", "Hylia_prasina")) 
#26.01016

CALIB <- makeChronosCalib(Phyllo, node="root", age.min=CladeAge)

PhylloC <- chronos(Phyllo, model="clock", calibration=CALIB)

plot(PhylloC, cex=0.5); axisPhylo()

## only crown viduidae age
CrownAge <- getMRCAage(PhylloC, tip=c("Phylloscopus_olivaceus", "Phylloscopus_neglectus"))

# get rid of outgroup 
PhylloC <- drop.tip(PhylloC, tip=c("Hylia_prasina"))

plot(PhylloC, cex=0.8, root.edge=TRUE); axisPhylo()


# Bind tree but specifying the position of the new node
PasBackTre19 <- bind.tree(PasBackTre18, PhylloC, where=which(PasBackTre18$tip.label =="Phylloscopus_collybita"), position=CrownAge)
plot(PasBackTre19, show.tip.label = F)

# check for duplicate
anyDuplicated(PasBackTre19$tip.label)
PasBackTre19 <- drop.tip(PasBackTre19, tip="Phylloscopus_collybita") # this isn't in fam tree but remove so no polytomy

### SAVE work ###

write.tree(PasBackTre19, file="PasBackTre19.tre")






##### Ploceidae #######
# De Silva 2019
# just rescale bc already calibrated
# It only has the clade, ie no outgroup, no clade origin age, only crown age in theory (not on backbone tho)
# 2 step using passeroidea tree to get ploceidae crown age

# extract ploc from pas/phlawd, get age of node bw ploc/taenio or vidua
# use to get age of ploc crown
# use crown age to calibrate just ploc tree
# graft ploc onto the ploceus branch, polytomy, rm OG

Ploc <- read.tree('ploceidae.tre')
plot(Ploc, cex=0.6)
summary(Ploc)

# fix names
Ploc$tip.label <- gsub(" ", '_', Ploc$tip.label)
Ploc$tip.label <- gsub("'", '', Ploc$tip.label)


# drop subspp
Phyllo <- drop.tip(Phyllo, tip=c("Seicercus_affinis_ocularis",
                                 "Seicercus_affinis_intermedius"))

# just in case but i don't think this is needed
for (tip in Ploc$tip.label) {
  if ((lengths(strsplit(tip, split="_"))>2)) {
    Ploc$tip.label[Ploc$tip.label==tip] <- sub("_[^_]+$", "", tip)
  }
}

# Check for duplicated tip labels
anyDuplicated(Ploc$tip.label)   



#Ploceus_cucullatus, Taeniopygia_guttata
PlocTaen <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Ploceus_cucullatus", "Taeniopygia_guttata")))

plot(PlocTaen, cex=0.5)

PlocOnly <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, c("Ploceus_cucullatus", "Amblyospiza_albifrons")))
plot(PlocOnly, cex=0.6)

#keep only ploc and outgroup
tips.to.drop <- setdiff(PlocTaen$tip.label, PlocOnly$tip.label)
tips.to.drop <- tips.to.drop[! tips.to.drop %in% c("Taeniopygia_guttata")]

# just ploc and one outgroup member
PlocTaen <- drop.tip(PlocTaen, tip=tips.to.drop)
plot(PlocTaen, cex=0.6)



# PlocVid <- extract.clade(Pas, node=getMRCA(Pas, tip=c("Ploceus_cucullatus", "Vidua_macroura")))
# plot(PlocVid, cex=0.4)


### CALIBRATE with backbone tree date ###

# age of node bw ploc and taenio
CladeAge <- getMRCAage(PasBackTre, tip=c("Ploceus_cucullatus", "Taeniopygia_guttata")) 
#21.95508

CALIB <- makeChronosCalib(PlocTaen, node="root", age.min=CladeAge)

PlocTaenC <- chronos(PlocTaen, model="clock", calibration=CALIB)

plot(PlocTaenC, cex=0.5); axisPhylo()

## only crown ploc age now
CrownAge <- getMRCAage(PlocTaenC, tip=c("Ploceus_cucullatus", "Amblyospiza_albifrons"))
# 20.9853

# now calibrate ploc tree w this age 
CALIB <- makeChronosCalib(Ploc, node="root", age.min=CrownAge)

PlocC <- chronos(Ploc, model="clock", calibration=CALIB)

#plot(PlocC)

# Bind tree but specifying the position of the new node
PasBackTre20 <- bind.tree(PasBackTre19, PlocC, where=which(PasBackTre19$tip.label =="Ploceus_cucullatus"), position=CrownAge)
plot(PasBackTre20, show.tip.label = F)


# check for duplicate
anyDuplicated(PasBackTre20$tip.label)
PasBackTre20 <- drop.tip(PasBackTre20, tip="Ploceus_cucullatus") # this isn't in fam tree but remove so no polytomy

### SAVE work ###

write.tree(PasBackTre20, file="PasBackTre20.tre")


# ######## just to see w other tree ###
# PlocVid <- extract.clade(Pas, node=getMRCA(Pas, tip=c("Ploceus_cucullatus", "Vidua_macroura")))
#  #plot(PlocVid, cex=0.4)
# PlocOnly <- extract.clade(Pas, node=getMRCA(Pas, c("Ploceus_cucullatus", "Amblyospiza_albifrons")))
#  #plot(PlocOnly, cex=0.6)
#  
#  #keep only ploc and outgroup
#  tips.to.drop <- setdiff(PlocVid$tip.label, PlocOnly$tip.label)
#  tips.to.drop <- tips.to.drop[! tips.to.drop %in% c("Vidua_macroura")]
#  
#  # just ploc and one outgroup member
#  PlocVid <- drop.tip(PlocVid, tip=tips.to.drop)
#  #plot(PlocTaen, cex=0.6)
# ### CALIBRATE with backbone tree date ###
# 
# # age of node bw ploc and taenio
# CladeAge <- getMRCAage(PasBackTre, tip=c("Ploceus_cucullatus", "Vidua_macroura")) 
# #21.95508
# 
# CALIB <- makeChronosCalib(PlocVid, node="root", age.min=CladeAge)
# 
# PlocVidC <- chronos(PlocVid, model="clock", calibration=CALIB)
# 
# #plot(PlocTaenC, cex=0.5); axisPhylo()
# 
# ## only crown ploc age now
# CrownAge <- getMRCAage(PlocVidC, tip=c("Ploceus_cucullatus", "Amblyospiza_albifrons"))
# # 20.80252








######## Alaudidae ######
# Alström et al. 2023
# calibrated so just rescale
# sister to Nicatoridae

# can get age of the alauda/nicator node from backbone
# use to calibrate phlawd
# get the age of crown Alaud from phlawd
# use to rescale Alaud fam tree
# graft on Alauda
Alaud <- read.tree("Fig2_BEAST_multilocus.trees")
plot(Alaud, cex=0.4)
Alaud <- drop.tip(Alaud, tip="Prinia_flavicans_U3818")

AlaudSpp <- as.data.frame(Alaud$tip.label)

# collapse spp
Alaud <- drop.tip(Alaud, tip=c("Alaemon_alaudipes_U5653",
                           "Alaemon_hamertoni_1002",
                           "Chersomanes_albofasciata_RSA100",
                           "Certhilauda_curvirostris_P215" ,               
                           "Certhilauda_curvirostris_P220",
                           "Ammomanes_cinctura_arenicolor2",
                           "Ammomanes_deserti_phoenicuroides2",
                           "Ramphocoris_clotbey_U5651",
                           "Pinarocorys_erythropygia_105748",
                           "Pinarocorys_nigricans_106597",
                           "Eremopterix_griseus_U2258",
                           "Eremopterix_melanauchen_U2260",
                           "Eremopterix_hova_FMNH449163",
                           "Alauda_arvensis_U587",
                           "Alauda_gulgula_U3267",
                           "Alauda_razae_SYSUS2265",
                           "Galerida_cristata_b",
                           "Galerida_randonii_CL2a",
                           "Galerida_malabarica_EF445430",
                           "Galerida_theklae_ruficolor",
                           "Lullula_arborea_U544",
                           "Spizocorys_fringillaris_U5659",
                           "Alaudala_cheleensis_U608",
                           "Alaudala_heinei_U585",
                           "Alaudala_raytal_U2201",
                           "Alaudala_rufescens_ad",
                           "Alaudala_somalica_AY165166",
                           "Chersophilus_duponti_U2255",
                           "Eremalauda_dunni_dunni_JL_Copete",
                           "Eremalauda_eremodites_U4578",
                           "Melanocorypha_bimaculata_U2283",
                           "Melanocorypha_calandra_U0583",
                           "Melanocorypha_yeltoniensis_U580",
                           "Melanocorypha_maxima_U588",
                           "Melanocorypha_mongolica_UWBM59839",
                           "Calandrella_c_cinerea_KX379963",
                           "Calandrella_c_rufipecta_KX379983",
                           "Calandrella_acutirostris_U577",
                           "Calandrella_dukhunensis_UWBM59838",
                           "Calandrella_brachydactyla_U5655",
                           "Calandrella_blanfordi_erlangeri",
                           "Calandrella_eremica_daaroodensis",
                           "Eremophila_alpestris_U4610",
                           "Eremophila_bilopha_U612",
                           "Eremophila_penicillata_U615",
                           "Eremophila_longirostris_U576",
                           "Calendulauda_africanoides_P175",
                           "Calendulauda_erythroclamys_U5662",
                           "Calendulauda_gilletti_U2808",
                           "Calendulauda_sabota_U2344",
                           "Mirafra_rufa_253",
                           "Heteromirafra_archeri_U2811",
                           "Heteromirafra_ruddi_U5661",
                           "Mirafra_affinis_U3268",
                           "Mirafra_assamica_U3269",
                           "Mirafra_erythroptera_U3271",
                           "Mirafra_erythrocephala_U3270",
                           "Mirafra_microptera_U3275",
                           "Mirafra_albicauda_9488",
                           "Mirafra_cheniana_U5759",
                           "Mirafra_javanica_williamsoni",
                           "Mirafra_cordofanica_1094",
                           "Mirafra_africana_5",
                           "Mirafra_tropicalis_1141",
                           "Mirafra_nyikae_1079",
                           "Mirafra_kabalii_1085",
                           "Mirafra_ashi_1144",
                           "Mirafra_somalica_1118",
                           "Mirafra_hypermetra_317",
                           "Mirafra_sharpii_1115",
                           "Mirafra_athi_U4595",
                           "Mirafra_fasciolata_U5758",
                           "Mirafra_kurrae_1109",
                           "Mirafra_rufocinnamomea_U5657",
                           "Mirafra_collaris_2634"))
                           
for (tip in Alaud$tip.label) {
  if ((lengths(strsplit(tip, split="_"))>2)) {
    Alaud$tip.label[Alaud$tip.label==tip] <- sub("_[^_]+$", "", tip)
  }
}

# don't know why it missed this one specifically
Alaud$tip.label[Alaud$tip.label=="Eremalauda_dunni_dunni_BMNH1932_8_6_268"] <- "Eremalauda_dunni"

# Check for duplicated tip labels
anyDuplicated(Alaud$tip.label)    

# remove outgroup entirely
Alaud <- extract.clade(Alaud, node=getMRCA(Alaud, c("Alauda_arvensis", "Alaemon_alaudipes")))


### get Alauda/Nicator node age
AlNic <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, c("Nicator_chloris", "Alauda_arvensis")))
plot(AlNic)

AlOnly <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, c("Alaemon_alaudipes", "Alauda_arvensis")))
plot(AlOnly)
#Alaemon_alaudipes

tips.to.drop <- setdiff(AlNic$tip.label, AlOnly$tip.label)
tips.to.drop <- tips.to.drop[! tips.to.drop %in% c("Nicator_chloris")]

AlNic <- drop.tip(AlNic, tip=tips.to.drop)

### CALIBRATE with backbone tree date ###

# age of node bw ploc and taenio
CladeAge <- getMRCAage(PasBackTre, tip=c("Alauda_arvensis", "Nicator_chloris")) 
#27.21791

CALIB <- makeChronosCalib(AlNic, node="root", age.min=CladeAge)

AlNicC <- chronos(AlNic, model="clock", calibration=CALIB)

plot(AlNicC, cex=0.5); axisPhylo()

## only crown Alaud age now
CrownAge <- getMRCAage(AlNicC, tip=c("Alauda_arvensis", "Alaemon_alaudipes"))
# 17.72489

##### now rescale Alaud tree w this age 
### Determine scaling factor ###

# use crown age

#Obtain tree age
TreeAge <- getMRCAage(Alaud, tip=c("Alauda_arvensis","Alaemon_alaudipes"))
# 18.54369


#Convert branch lengths
Alaud$edge.length <- Alaud$edge.length*CrownAge/TreeAge


plot(Alaud, cex=0.8, root.edge=TRUE); axisPhylo()


# Bind tree but specifying the position of the new node
PasBackTre21 <- bind.tree(PasBackTre20, Alaud, where=which(PasBackTre20$tip.label =="Alauda_arvensis"), position=CrownAge)
plot(PasBackTre21, show.tip.label = F)

# check for duplicate
anyDuplicated(PasBackTre21$tip.label)

# removing duplicates 
which(PasBackTre21$tip.label %in% c("Alauda_arvensis"))
PasBackTre21 <- drop.tip(PasBackTre21, tip=4130) 

### SAVE work ###

write.tree(PasBackTre21, file="PasBackTre21.tre")                           
                           


####### Stenostiridae, Remizidae, and Paridae ######
 
SteRePa <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, c("Elminia_nigromitrata", "Parus_major")))
#plot(SteRePa)

# collapse spp
SteRePa <- drop.tip(SteRePa, tip=c("Sylviparus_modestus_modestus",
                                   "Pardaliparus_elegans_mindanensis", 
                                   "Pardaliparus_elegans_elegans", 
                                   "Pardaliparus_elegans_montigenus",
                                   "Pardaliparus_elegans_gilliardi" ,
                                   "Pardaliparus_elegans_albescens",
                                   "Periparus_ater_cypriotes",
                                   "Periparus_ater_ater",
                                   "Periparus_ater_aemodius",
                                   "Periparus_rubidiventris_beavani" ,     
                                   "Periparus_rubidiventris_rubidiventris",
                                   "Sittiparus_varius_orii",
                                   "Sittiparus_varius_amamii",
                                   "Sittiparus_varius_varius" ,            
                                   "Sittiparus_varius_namiyei",
                                   "Sittiparus_varius_yakushimensis",
                                   "Sittiparus_varius_sunsunpi",
                                   "Poecile_palustris_palustris" ,         
                                   "Poecile_palustris_hensoni",
                                   "Poecile_palustris_hellmayri",          
                                   "Poecile_palustris_brevirostris" ,
                                   "Poecile_palustris_jeholicus",
                                   "Poecile_montanus_affinis",             
                                    "Poecile_songarus_stotzneri",
                                   "Poecile_montanus_salicarius",
                                  "Poecile_montanus_kamtschatkensis",
                                  "Poecile_montanus_sachalinensis",
                                  "Poecile_montanus_montanus",
                                  "Poecile_montanus_baicalensis",
                                  "Poecile_montanus_restrictus",
                                  "Poecile_montanus_rhenanus",
                                  "Poecile_gambeli_gambeli",
                                  "Poecile_gambeli_baileyi",
                                  "Poecile_hudsonicus_littoralis",
                                  "Poecile_carolinensis_extimus",
                                  "Poecile_lugubris_anatoliae",
                                  "Poecile_lugubris_lugens",
                                  "Poecile_lugubris_lugubris",
                                  "Lophophanes_cristatus_mitratus",
                                  "Melaniparus_leucomelas_insignis",
                                  "Melaniparus_niger_niger",
                                  "Melaniparus_fasciiventer_fasciiventer",
                                  "Parus_major_major",
                                  "Parus_minor_minor",
                                  "Machlolophus_spilonotus_basileus",
                                  "Machlolophus_spilonotus_subviridis",
                                  "Machlolophus_spilonotus_rex",
                                  "Machlolophus_xanthogenys_xanthogenys",
                                  "Cyanistes_cyanus_tianschanicus",
                                  "Cyanistes_cyanus_cyanus",
                                  "Cyanistes_cyanus_yenisseensis" ,
                                  "Cyanistes_caeruleus_ogliastrae",       
                                  "Cyanistes_caeruleus_obscurus",
                                  "Cyanistes_caeruleus_calamensis",
                                  "Cyanistes_caeruleus_caeruleus",
                                  "Cyanistes_caeruleus_satunini" ,
                                  "Cyanistes_caeruleus_raddei" ,
                                  "Cyanistes_cyanus_hyperriphaeus",
                                  "Cyanistes_ultramarinus_cyrenaicae",
                                  "Cyanistes_teneriffae_teneriffae",
                                  "Cyanistes_teneriffae_hedwigii",
                                  "Remiz_coronatus_stoliczkae",
                                  "Remiz_pendulinus_pendulinus",
                                  "Remiz_macronyx_ssaposhnikowi",
                                  "Remiz_pendulinus_menzbieri",
                                  "Remiz_macronyx_neglectus",
                                  "Remiz_pendulinus_caspius"
                                   ))


for (tip in SteRePa$tip.label) {
  if ((lengths(strsplit(tip, split="_"))>2)) {
    SteRePa$tip.label[SteRePa$tip.label==tip] <- sub("_[^_]+$", "", tip)
  }
}


# Check for duplicated tip labels
anyDuplicated(SteRePa$tip.label)                             


### CALIBRATE with backbone tree date ###

CladeAge <- getMRCAage(PasBackTre, tip=c("Elminia_nigromitratus", "Parus_major")) 
#31.76514

 CALIB <- makeChronosCalib(SteRePa, node="root", age.min=CladeAge) 
 
#SteRePaC <- chronos(SteRePa, model="clock", calibration=CALIB) # clock not working??
#Using model = "clock" is actually a short-cut to model = "discrete" and setting nb.rate.cat = 1 in the list passed to control.
SteRePaC <- chronos(SteRePa, model="discrete", control=chronos.control(nb.rate.cat=1), calibration=CALIB)

plot(SteRePaC, cex=0.5); axisPhylo()

# remove 1 sp so clade of 3 becomes 2
PasBackTre22 <- drop.tip(PasBackTre21, tip="Remiz_pendulinus")
# remove 2 tips from backbone
PasBackTre22 <- drop.tip(PasBackTre22, tip=c("Elminia_nigromitratus", "Parus_major"),  trim.internal=FALSE)

# Graft onto Backbone
PasBackTre22 <- bind.tree(PasBackTre22, SteRePaC, where=Ntip(PasBackTre22))

# check for duplicate
anyDuplicated(PasBackTre22$tip.label)

plot(PasBackTre22, show.tip.label = F)


### SAVE work ###

write.tree(PasBackTre22, file="PasBackTre22.tre")   





######## Chlorospeidae and Irenidae ######
ChlorIr <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Chloropsis_cochinchinensis", "Irena_cyanogastra")))
plot(ChlorIr)

ChlorIr <- drop.tip(ChlorIr, tip=c("Chloropsis_cochinchinensis_viridinucha",
                                   "Chloropsis_cyanopogon_cyanopogon",
                                   "Irena_cyanogastra_cyanogastra",
                                   "Irena_cyanogastra_hoogstraali",
                                   "Irena_puella_crinigera"
                                   ))

for (tip in ChlorIr$tip.label) {
  if ((lengths(strsplit(tip, split="_"))>2)) {
    ChlorIr$tip.label[ChlorIr$tip.label==tip] <- sub("_[^_]+$", "", tip)
  }
}


# Check for duplicated tip labels
anyDuplicated(ChlorIr$tip.label)                             


### CALIBRATE with backbone tree date ###

CladeAge <- getMRCAage(PasBackTre, tip=c("Chloropsis_cochinchinensis", "Irena_cyanogaster")) 
#27.64838

CALIB <- makeChronosCalib(ChlorIr, node="root", age.min=CladeAge) 

#Using model = "clock" is actually a short-cut to model = "discrete" and setting nb.rate.cat = 1 in the list passed to control.
ChlorIrC <- chronos(ChlorIr, model="discrete", control=chronos.control(nb.rate.cat=1), calibration=CALIB)


plot(ChlorIrC, cex=0.5); axisPhylo()


# remove 2 tips from backbone
PasBackTre23 <- drop.tip(PasBackTre22, tip=c("Chloropsis_cochinchinensis", "Irena_cyanogaster"),  trim.internal=FALSE)

# Graft onto Backbone
PasBackTre23 <- bind.tree(PasBackTre23, ChlorIrC, where=Ntip(PasBackTre23))

# check for duplicate
anyDuplicated(PasBackTre23$tip.label)

plot(PasBackTre23, show.tip.label = F)


### SAVE work ##

write.tree(PasBackTre23, file="PasBackTre23.tre")   






######## Hypocoliidae, Bombycillidae, Ptiliogonatidae, ~~Mohoidae~~, Dulidae, Hylocitreidae ########
# small, pitiful, also RIP Mohoidae	

HBPMDH <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Dulus_dominicus", "Hypocolius_ampelinus")))
plot(HBPMDH)

# collapse sp
HBPMDH <- drop.tip(HBPMDH, tip=c("Bombycilla_garrulus_garrulus"))


for (tip in HBPMDH$tip.label) {
  if ((lengths(strsplit(tip, split="_"))>2)) {
    HBPMDH$tip.label[HBPMDH$tip.label==tip] <- sub("_[^_]+$", "", tip)
  }
}


# Check for duplicated tip labels
anyDuplicated(HBPMDH$tip.label)                             


### CALIBRATE with backbone tree date ###

CladeAge <- getMRCAage(PasBackTre, tip=c("Hypocolius_ampelinus", "Hylocitrea_bonensis")) 
#.44903

CALIB <- makeChronosCalib(HBPMDH, node="root", age.min=CladeAge) 

#SteRePaC <- chronos(SteRePa, model="clock", calibration=CALIB) # clock not working??
#Using model = "clock" is actually a short-cut to model = "discrete" and setting nb.rate.cat = 1 in the list passed to control.
HBPMDHC <- chronos(HBPMDH, model="discrete", control=chronos.control(nb.rate.cat=1), calibration=CALIB)

plot(HBPMDHC, cex=0.5); axisPhylo()

# remove spp so clade becomes 2
PasBackTre24 <- drop.tip(PasBackTre23, tip=c("Bombycilla_garrulus", "Ptilogonys_cinereus", "Moho_nobilis", "Dulus_dominicus"))
# remove 2 tips from backbone
PasBackTre24 <- drop.tip(PasBackTre24, tip=c("Hypocolius_ampelinus", "Hylocitrea_bonensis"),  trim.internal=FALSE)

# Graft onto Backbone
PasBackTre24 <- bind.tree(PasBackTre24, HBPMDHC, where=Ntip(PasBackTre24))

# check for duplicate
anyDuplicated(PasBackTre24$tip.label)

plot(PasBackTre24, show.tip.label = F)


### SAVE work ###

write.tree(PasBackTre24, file="PasBackTre24.tre")   






#### Acrocephalidae #######
Acro <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Acrocephalus_newtoni", "Hirundo_rustica")))
#plot(Acro, cex=0.4)

# Acro only
AcroOnly <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Acrocephalus_newtoni", "Nesillas_aldabrana")))

tips.to.drop <- setdiff(Acro$tip.label, AcroOnly$tip.label)
tips.to.drop <- tips.to.drop[! tips.to.drop %in% c("Hirundo_rustica")]

Acro <- drop.tip(Acro, tip=tips.to.drop)
#plot(Acro, cex=0.4)


# collapse spp
Acro <- drop.tip(Acro, tip=c("Hippolais_icterina_icterina",
                             "Acrocephalus_agricola_septimus",
                             "Acrocephalus_agricola_agricola",
                             "Acrocephalus_agricola_capistrata",
                             "Acrocephalus_concinens_concinens",
                             "Acrocephalus_scirpaceus_fuscus" ,         
                             "Acrocephalus_scirpaceus_scirpaceus",
                             "Acrocephalus_baeticatus_guiersi",
                             "Acrocephalus_baeticatus_hallae",
                             "Acrocephalus_melanopogon_mimicus",
                             "Acrocephalus_melanopogon_melanopogon",
                             "Acrocephalus_familiaris_kingi",
                             "Acrocephalus_kerearako_kaoko",
                             "Acrocephalus_musae_garretti",
                             "Acrocephalus_mendanae_idae",
                             "Acrocephalus_percernis_idae",
                             "Acrocephalus_luscinius_astrolabii",
                             "Acrocephalus_percernis_aquilonis",
                             "Acrocephalus_percernis_postremus",
                            "Acrocephalus_mendanae_postremus",
                            "Acrocephalus_mendanae_aquilonis",
                           "Acrocephalus_mendanae_percernis",
                           "Acrocephalus_stentoreus_sumbae",
                           "Acrocephalus_australis_australis",
                           "Acrocephalus_australis_gouldi",
                           "Acrocephalus_aequinoctialis_pistor",
                           "Acrocephalus_luscinius_nijoi",
                           "Acrocephalus_mendanae_dido",
                           "Acrocephalus_percernis_percernis",
                           "Acrocephalus_mendanae_fatuhivae",
                           "Acrocephalus_mendanae_consobrina",
                           "Acrocephalus_stentoreus_harterti",
                           "Acrocephalus_stentoreus_amyae",
                           "Acrocephalus_stentoreus_levantinus",
                           "Acrocephalus_stentoreus_brunnescens",
                           "Acrocephalus_arundinaceus_arundinaceus",
                           "Acrocephalus_arundinaceus_zarudnyi",
                           "Acrocephalus_gracilirostris_parvus",
                           "Acrocephalus_rufescens_senegalensis",
                           "Acrocephalus_rufescens_rufescens",
                           "Chloropeta_gracilirostris_gracilirostris",
                           "Iduna_pallida_elaeica",
                           "Chloropeta_natalensis_natalensis",
                           "Chloropeta_natalensis_massaica"
                             ))




for (tip in Acro$tip.label) {
  if ((lengths(strsplit(tip, split="_"))>2)) {
    Acro$tip.label[Acro$tip.label==tip] <- sub("_[^_]+$", "", tip)
  }
}

# Check for duplicated tip labels
anyDuplicated(Acro$tip.label)

### CALIBRATE with backbone tree date ###

CladeAge <- getMRCAage(PasBackTre, tip=c("Acrocephalus_newtoni", "Hirundo_rustica")) # node between fringillidae and rest
#27.34738

CALIB <- makeChronosCalib(Acro, node="root", age.min=CladeAge)

AcroC <- chronos(Acro, model="discrete", control=chronos.control(nb.rate.cat=1), calibration=CALIB)

plot(AcroC, cex=0.5); axisPhylo()

# only crown fring age
CrownAge <- getMRCAage(AcroC, tip=c("Acrocephalus_newtoni", "Nesillas_aldabrana"))

# get rid of outgroup 
AcroC <- drop.tip(AcroC, tip=c("Hirundo_rustica"))

plot(AcroC, cex=0.3, root.edge=TRUE); axisPhylo()


# Bind tree but specifying the position of the new node
PasBackTre25 <- bind.tree(PasBackTre24, AcroC, where=which(PasBackTre24$tip.label =="Acrocephalus_newtoni"), position=CrownAge)
plot(PasBackTre25, show.tip.label = F)

# check for duplicate
anyDuplicated(PasBackTre25$tip.label)

# removing duplicates
PasBackTre25 <- drop.tip(PasBackTre25, tip=4500) # Acrocephalus palustris in Locustella for some reason?
PasBackTre25 <- drop.tip(PasBackTre25, tip=4507) # Acrocephalus newtoni, og, run these in order

### SAVE work ###
write.tree(PasBackTre25, file="PasBackTre25.tre")







###### Orthonychidae + Pomatostomidae #####
OrthPom <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Pomatostomus_ruficeps", "Orthonyx_temminckii")))

OrthPom <- drop.tip(OrthPom, tip=c("Orthonyx_spaldingii_melasmenus", "Orthonyx_spaldingii_spaldingii"))

for (tip in OrthPom$tip.label) {
  if ((lengths(strsplit(tip, split="_"))>2)) {
    OrthPom$tip.label[OrthPom$tip.label==tip] <- sub("_[^_]+$", "", tip)
  }
}


# Check for duplicated tip labels
anyDuplicated(OrthPom$tip.label)                             


### CALIBRATE with backbone tree date ###

CladeAge <- getMRCAage(PasBackTre, tip=c("Pomatostomus_isidorei", "Orthonyx_teminckii")) #temminckii
#33.86948

CALIB <- makeChronosCalib(OrthPom, node="root", age.min=CladeAge) 

#SteRePaC <- chronos(SteRePa, model="clock", calibration=CALIB) # clock not working??
#Using model = "clock" is actually a short-cut to model = "discrete" and setting nb.rate.cat = 1 in the list passed to control.
OrthPomC <- chronos(OrthPom, model="discrete", control=chronos.control(nb.rate.cat=1), calibration=CALIB)

plot(OrthPomC, cex=0.5); axisPhylo()

# remove 2 tips from backbone
PasBackTre26 <- drop.tip(PasBackTre25, tip=c("Pomatostomus_isidorei", "Orthonyx_teminckii"),  trim.internal=FALSE)

# Graft onto Backbone
PasBackTre26 <- bind.tree(PasBackTre26, OrthPomC, where=Ntip(PasBackTre26))

# check for duplicate
anyDuplicated(PasBackTre26$tip.label)

# fix that one
PasBackTre26$tip.label[which(PasBackTre26$tip.label %in% c("Orthonyx_teminckii"))] <- "Orthonyx_temminckii"

plot(PasBackTre26, show.tip.label = F)


### SAVE work###
write.tree(PasBackTre26, file="PasBackTre26.tre")












######## Melanocharitidae ######
# low priority, can use phlawd tree 
# Mila 2021
# calibrated so just rescale
# no outgroup :(
Melano <- read.tree("melanocharitidae.tre")
plot(Melano)

#Melano <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Toxorhamphus_poliopterus", "Melanocharis_nigra")))

# fix names
#Melano$tip.label<- gsub('[[:digit:]]+', '', Melano$tip.label)
#Melano$tip.label <- gsub("'", '', Melano$tip.label)
#Melano$tip.label <- str_squish(Melano$tip.label)
#Melano$tip.label <- gsub(" ", '_', Melano$tip.label)
Melano$tip.label

# manually replace names
Melano$tip.label[which(Melano$tip.label %in% c("MelanochariscitreolaLGR167"))] <- "Melanocharis_citreola"
Melano$tip.label[which(Melano$tip.label %in% c("MelanocharisstriativentrisB01008"))] <- "Melanocharis_striativentris"
Melano$tip.label[which(Melano$tip.label %in% c("MelanocharislongicaudaB26041"))] <- "Melanocharis_longicauda"
Melano$tip.label[which(Melano$tip.label %in% c("Melanocharisversteri10MV"))] <- "Melanocharis_versteri"
Melano$tip.label[which(Melano$tip.label %in% c("MelanocharisarfakianaMELARF657308"))] <- "Melanocharis_arfakiana"
Melano$tip.label[which(Melano$tip.label %in% c("Rhamphochariscrassirostris150300"))] <- "Melanocharis_crassirostris"
Melano$tip.label[which(Melano$tip.label %in% c("Melanocharisnigra150402"))] <- "Melanocharis_nigra"
Melano$tip.label[which(Melano$tip.label %in% c("Oedistomailiolophus150392"))] <- "Oedistoma_iliolophus"
Melano$tip.label[which(Melano$tip.label %in% c("OedistomapygmaeumZ44439"))] <- "Oedistoma_pygmaeum"
Melano$tip.label[which(Melano$tip.label %in% c("ToxorhamphusnovaeguineaeLGR144"))] <- "Toxorhamphus_novaeguineae"

# drop other tips
# i am proud of myself
Melano <- drop.tip(Melano, tip=Melano$tip.label[lengths(strsplit(Melano$tip.label, split="_"))<2])

### now get melano age
MelanoPhlawd <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Melanocharis_nigra", "Ptilorrhoa_caerulescens")))

MelanoOnlyPhlawd <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, c("Melanocharis_nigra", "Toxorhamphus_novaeguineae")))

#keep only ploc and outgroup
tips.to.drop <- setdiff(MelanoPhlawd$tip.label, MelanoOnlyPhlawd$tip.label)
tips.to.drop <- tips.to.drop[! tips.to.drop %in% c("Ptilorrhoa_caerulescens")]

# just ploc and one outgroup member
MelanoPhlawd <- drop.tip(MelanoPhlawd, tip=tips.to.drop)
plot(MelanoPhlawd, cex=0.6)

### CALIBRATE with backbone tree date ###

# age of node bw melano and outgroup
CladeAge <- getMRCAage(PasBackTre, tip=c("Melanocharis_nigra", "Ptilorrhoa_caerulescens")) 
#30.27485

CALIB <- makeChronosCalib(MelanoPhlawd, node="root", age.min=CladeAge)

MelanoPhlawdC <- chronos(MelanoPhlawd,model="discrete", control=chronos.control(nb.rate.cat=1), calibration=CALIB)

plot(MelanoPhlawdC, cex=0.5); axisPhylo()

## only crown melano age now
CrownAge <- getMRCAage(MelanoPhlawdC, tip=c("Melanocharis_nigra", "Toxorhamphus_novaeguineae"))
# 16.64464

##### now rescale Melano tree w this age 
### Determine scaling factor ###

# use crown age

#Obtain tree age
TreeAge <- getMRCAage(Melano, tip=c("Melanocharis_nigra","Toxorhamphus_novaeguineae"))
# 18.86064


#Convert branch lengths
Melano$edge.length <- Melano$edge.length*CrownAge/TreeAge

# Bind tree but specifying the position of the new node
PasBackTre27 <- bind.tree(PasBackTre26, Melano, where=which(PasBackTre26$tip.label =="Melanocharis_nigra"), position=CrownAge)
plot(PasBackTre27, show.tip.label = F)


# check for duplicate
anyDuplicated(PasBackTre27$tip.label)
PasBackTre27 <- drop.tip(PasBackTre27, tip=1615) 

### SAVE work ###

write.tree(PasBackTre27, file="PasBackTre27.tre")







###### Nectariniidae and Dicaeidae ######

NectDi <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Nectarinia_olivacea", "Dicaeum_hypoleucum")))

#NectOnly <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Nectarinia_olivacea", "Arachnothera_juliae")))

#keep only Phyllo and outgroup
#tips.to.drop <- setdiff(Nect$tip.label, NectOnly$tip.label)
#tips.to.drop <- tips.to.drop[! tips.to.drop %in% c("Muscicapa_ferruginea")]

#Nect <- drop.tip(Nect, tip=tips.to.drop)
#plot(Nect, cex=0.6)

# drop subspp
NectDi <- drop.tip(NectDi, tip=c("Kurochkinegramma_hypogrammicum_hypogrammicum",
                                 "Arachnothera_chrysogenys_chrysogenys",
                             "Arachnothera_longirostra_buettikoferi",
                             "Cinnyris_reichenowi_preussi",
                             "Cinnyris_whytei_skye",
                             "Cinnyris_whytei_whytei",
                             "Cinnyris_notatus_moebii",
                             "Cinnyris_notatus_notatus",
                             "Cinnyris_notatus_voeltzkowi",
                             "Cinnyris_jugularis_aurora",
                             "Cinnyris_jugularis_jugularis",
                             "Cinnyris_jugularis_obscurior",
                             "Cinnyris_sovimanga_sovimanga",
                             "Cinnyris_abbotti_buchenorum",
                             "Cinnyris_sovimanga_aldabrensis",
                             "Anthreptes_malacensis_malacensis",
                             "Anthreptes_malacensis_paraguae",
                             "Leptocoma_brasiliana_brasiliana",
                             "Aethopyga_siparaja_siparaja",
                             "Dicaeum_hypoleucum_cagayanense" ,
                             "Dicaeum_hypoleucum_pontifex",
                             "Prionochilus_maculatus_maculatus",
                             "Prionochilus_olivaceus_olivaceus" ,
                             "Prionochilus_olivaceus_parsoni",
                             "Prionochilus_olivaceus_samarensis",
                             "Dicaeum_trigonostigma_xanthopygium",          
                             "Dicaeum_trigonostigma_cinereigulare",
                             "Dicaeum_trigonostigma_dorsale",
                             "Dicaeum_cruentatum_nigrimentum",
                             "Dicaeum_minullum_borneanum"
                             
                             ))

for (tip in NectDi$tip.label) {
  if ((lengths(strsplit(tip, split="_"))>2)) {
    NectDi$tip.label[NectDi$tip.label==tip] <- sub("_[^_]+$", "", tip)
  }
}

# Check for duplicated tip labels
anyDuplicated(NectDi$tip.label)                             


### CALIBRATE with backbone tree date ###

CladeAge <- getMRCAage(PasBackTre, tip=c("Leptocoma_aspasia", "Dicaeum_aeneum")) 
#24.71527

CALIB <- makeChronosCalib(NectDi, node="root", age.min=CladeAge)

NectDiC <- chronos(NectDi, model="discrete", control=chronos.control(nb.rate.cat=1), calibration=CALIB)

plot(NectDiC, cex=0.5); axisPhylo()


# remove 2 tips from backbone
PasBackTre28 <- drop.tip(PasBackTre27, tip=c("Leptocoma_aspasia", "Dicaeum_aeneum"),  trim.internal=FALSE)

# Graft onto Backbone
PasBackTre28 <- bind.tree(PasBackTre28, NectDiC, where=Ntip(PasBackTre28))

# check for duplicate
anyDuplicated(PasBackTre28$tip.label)

# in for some reason?
PasBackTre28 <- drop.tip(PasBackTre28, tip=5165) #Cyornis_rufigastra

plot(PasBackTre28, show.tip.label = F)


### SAVE work###
write.tree(PasBackTre28, file="PasBackTre28.tre")







####### Aegithalidae, Scotocercidae, Hyliidae ######
AeScoHy <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Hylia_prasina", "Aegithalos_iouschensis")))
plot(AeScoHy, cex=0.5)

# phylloscopidae nested in this clade in phlawd tree but i took it out according to backbone
AeScoHy <- drop.tip(AeScoHy, tip=(extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Phylloscopus_neglectus", "Phylloscopus_olivaceus"))))$tip.label)

# collapse spp
AeScoHy <- drop.tip(AeScoHy, tip=c("Leptopoecile_sophiae_obscura",
                                   "Psaltriparus_minimus_plumbeus",
                                   "Psaltriparus_minimus_minimus",
                                   "Cettia_cetti_orientalis",
                                   "Cettia_cetti_cetti",
                                   "Cettia_cetti_albiventris",
                                   "Oligura_castaneocoronata", # synonym w the Cettia sp?
                                   "Phyllergates_cucullatus_hedymeles" ,    
                                   "Phyllergates_cucullatus_cucullatus",
                                   "Phyllergates_cucullatus_cinereicollis",
                                   "Phyllergates_cucullatus_philippinus"   ,
                                   "Phyllergates_cucullatus_coronatus",
                                   "Horornis_fortipes_fortipes",
                                   "Cettia_fortipes_fortipes",
                                   "Horornis_fortipes_davidianus",
                                   "Cettia_fortipes_davidiana",
                                   "Horornis_fortipes_pallidus",            
                                   "Cettia_vulcania_flaviventris" ,
                                   "Cettia_vulcania_oreophila",
                                   "Cettia_vulcania_vulcania",              
                                   "Cettia_vulcania_intricata",
                                   "Cettia_vulcania_oblita",
                                   "Horornis_flavolivaceus_flavolivaceus",  
                                   "Horornis_flavolivaceus_weberi",
                                   "Horornis_acanthizoides_acanthizoides",
                                   "Horornis_acanthizoides_concolor",
                                   "Abroscopus_superciliaris_superciliaris",
                                   "Abroscopus_albogularis_fulvifacies"
                                   ))

for (tip in AeScoHy$tip.label) {
  if ((lengths(strsplit(tip, split="_"))>2)) {
    AeScoHy$tip.label[AeScoHy$tip.label==tip] <- sub("_[^_]+$", "", tip)
  }
}

# Check for duplicated tip labels
anyDuplicated(AeScoHy$tip.label)                             


### CALIBRATE with backbone tree date ###

CladeAge <- getMRCAage(PasBackTre, tip=c("Hylia_prasina", "Aegithalos_iouschensis")) 
#21.84884

CALIB <- makeChronosCalib(AeScoHy, node="root", age.min=CladeAge)

AeScoHyC <- chronos(AeScoHy, model="discrete", control=chronos.control(nb.rate.cat=1), calibration=CALIB)

plot(AeScoHyC, cex=0.5); axisPhylo()

# make clade of 3 into 2 tip
PasBackTre29 <- drop.tip(PasBackTre28, tip=c("Cettia_brunnifrons"))
# remove 2 tips from backbone
PasBackTre29 <- drop.tip(PasBackTre29, tip=c("Aegithalos_iouschensis", "Hylia_prasina"),  trim.internal=FALSE)

# Graft onto Backbone
PasBackTre29 <- bind.tree(PasBackTre29, AeScoHyC, where=Ntip(PasBackTre29))

# check for duplicate
anyDuplicated(PasBackTre29$tip.label)


plot(PasBackTre29, show.tip.label = F)


### SAVE work###
write.tree(PasBackTre29, file="PasBackTre29.tre")










###### Hirundinidae #####
Hirun <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Hirundo_rustica", "Aegithalos_iouschensis")))

HirunOnly <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Hirundo_rustica", "Pseudochelidon_eurystomina")))

tips.to.drop <- setdiff(Hirun$tip.label, HirunOnly$tip.label)
tips.to.drop <- tips.to.drop[! tips.to.drop %in% c("Aegithalos_iouschensis")]

Hirun <- drop.tip(Hirun, tip=tips.to.drop)
#plot(Hirun, cex=0.4)


# collapse spp
Hirun <- drop.tip(Hirun, tip=c("Psalidoprocne_holomelas_ruwenzori",
                               "Psalidoprocne_orientalis_reichenowi",
                               "Hirundo_tahitica_javanica",
                               "Hirundo_rustica_rustica",
                               "Hirundo_rustica_transitiva",
                               "Hirundo_rustica_savignii" ,
                               "Hirundo_rustica_gutturalis",
                               "Hirundo_rustica_tytleri",
                               "Hirundo_rustica_erythrogaster",
                               "Cecropis_daurica_rufula",
                               "Delichon_urbicum_urbicum",
                               "Progne_subis_arboricola",
                               "Progne_subis_subis",
                               "Progne_chalybea_macrorhamphus",
                               "Progne_chalybea_chalybea",
                               "Pygochelidon_cyanoleuca_patagonica",
                               "Riparia_riparia_shelleyi",
                               "Riparia_riparia_riparia",
                               "Riparia_diluta_tibetana",
                               "Riparia_diluta_diluta",
                               "Riparia_diluta_fohkienensis" ,
                               "Riparia_diluta_indica" ))




for (tip in Hirun$tip.label) {
  if ((lengths(strsplit(tip, split="_"))>2)) {
    Hirun$tip.label[Hirun$tip.label==tip] <- sub("_[^_]+$", "", tip)
  }
}

# Check for duplicated tip labels
anyDuplicated(Hirun$tip.label)

### CALIBRATE with backbone tree date ###

CladeAge <- getMRCAage(PasBackTre, tip=c("Aegithalos_iouschensis", "Hirundo_rustica")) # 
#25.12135

CALIB <- makeChronosCalib(Hirun, node="root", age.min=CladeAge)

HirunC <- chronos(Hirun, model="discrete", control=chronos.control(nb.rate.cat=1), calibration=CALIB)

plot(HirunC, cex=0.5); axisPhylo()

# only crown  age
CrownAge <- getMRCAage(HirunC, tip=c("Hirundo_rustica", "Pseudochelidon_eurystomina"))

# get rid of outgroup 
HirunC <- drop.tip(HirunC, tip=c("Aegithalos_iouschensis"))

plot(HirunC, cex=0.3, root.edge=TRUE); axisPhylo()


# Bind tree but specifying the position of the new node
PasBackTre30 <- bind.tree(PasBackTre29, HirunC, where=which(PasBackTre29$tip.label =="Hirundo_rustica"), position=CrownAge)
plot(PasBackTre30, show.tip.label = F)

# check for duplicate
anyDuplicated(PasBackTre30$tip.label)

# removing duplicates
PasBackTre30 <- drop.tip(PasBackTre30, tip=4549) # Hirundo rustica

### SAVE work ###
write.tree(PasBackTre30, file="PasBackTre30.tre")





##### Macrosphenidae ####

Macro <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Sphenoeacus_afer", "Megalurus_palustris")))
plot(Macro)


MacroOnly <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Sphenoeacus_afer", "Sylvietta_ruficapilla")))
plot(MacroOnly) # Macrosphenus really far away for some reason, so excluded

tips.to.drop <- setdiff(Macro$tip.label, MacroOnly$tip.label)
tips.to.drop <- tips.to.drop[! tips.to.drop %in% c("Megalurus_palustris")]

Macro <- drop.tip(Macro, tip=tips.to.drop)
plot(Macro)


# collapse spp
Macro <- drop.tip(Macro, tip=c("Anthoscopus_musculus", # idk why in here
                               "Bradypterus_victorini", # idk why in here
                               "Melocichla_mentalis_mentalis",
                               "Sylvietta_whytii_whytii"))




for (tip in Macro$tip.label) {
  if ((lengths(strsplit(tip, split="_"))>2)) {
    Macro$tip.label[Macro$tip.label==tip] <- sub("_[^_]+$", "", tip)
  }
}

# Check for duplicated tip labels
anyDuplicated(Macro$tip.label)

### CALIBRATE with backbone tree date ###

CladeAge <- getMRCAage(PasBackTre, tip=c("Sphenoeacus_afer", "Megalurus_palustris")) # 
#28.11144

CALIB <- makeChronosCalib(Macro, node="root", age.min=CladeAge)

MacroC <- chronos(Macro, model="clock", calibration=CALIB)

plot(MacroC, cex=0.5); axisPhylo()

# only crown  age
CrownAge <- getMRCAage(MacroC, tip=c("Sphenoeacus_afer", "Sylvietta_ruficapilla"))

# get rid of outgroup 
MacroC <- drop.tip(MacroC, tip=c("Megalurus_palustris"))

plot(MacroC, cex=0.3, root.edge=TRUE); axisPhylo()


# Bind tree but specifying the position of the new node
PasBackTre35 <- bind.tree(PasBackTre34, MacroC, where=which(PasBackTre34$tip.label =="Sphenoeacus_afer"), position=CrownAge)
plot(PasBackTre35, show.tip.label = F)

# check for duplicate
anyDuplicated(PasBackTre35$tip.label)

# removing duplicates
PasBackTre35 <- drop.tip(PasBackTre35, tip=4705) # Sphenoeacus_afer

### SAVE work ###
write.tree(PasBackTre35, file="PasBackTre35.tre")













###### Sittidae and Certhiidae, also Tichodroma, plus polioptilidae and troglodytidae? #####
# just keep polytomy, use the topology within phlawd tree

SitCert <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Sitta_carolinensis", "Polioptila_guianensis")))
plot(SitCert, cex=0.4)

SitCert<-(root(SitCert, outgroup="Tichodroma_muraria"))
plot(SitCert, cex=0.4)


# certhia, salpornis

# collapse spp
SitCert <- drop.tip(SitCert, tip=c("Certhia_himalayana_himalayana",                   
                                   "Certhia_himalayana_yunnanensis", 
                                   "Certhia_brachydactyla_mauritanica",       
                                   "Certhia_familiaris_bianchii",
                                   "Certhia_familiaris_britannica",           
                                   "Certhia_familiaris_macrodactyla",
                                   "Certhia_familiaris_familiaris",           
                                   "Certhia_familiaris_tianschanica",
                                   "Certhia_familiaris_daurica",
                                   "Certhia_familiaris_corsa",                
                                   "Certhia_familiaris_caucasica",           
                                   "Certhia_brachydactyla_megarhynchos",
                                   "Thryothorus_maculipectus_canobrunneus" ,  
                                   "Thryothorus_rutilus_hyperythrus" ,
                                   "Pheugopedius_genibarbis_bolivianus" ,     
                                   "Pheugopedius_genibarbis_genibarbis" ,         
                                    "Pheugopedius_coraya" ,                    
                                   "Thryothorus_rutilus" ,
                                  "Thryothorus_modestus_zeledoni",           
                                   "Thryothorus_modestus_modestus",
                                   "Cantorchilus_superciliaris_superciliaris",
                                   "Thryothorus_longirostris_bahiae" ,    
                                   "Thryothorus_nigricapillus", # removed bc not grouped w the others
                                    "Thryothorus_nigricapillus_schotti",
                                   "Thryothorus_nigricapillus_costaricensis" ,
                                   "Thryothorus_nigricapillus_connectens",
                                    "Uropsila_leucogastra_leucogastra"  ,      
                                   "Henicorhina_leucophrys_leucophrys",              
                                    "Henicorhina_leucosticta" ,              
                                   "Thryophilus_rufalbus_castanonotus",                 
                                    "Thryophilus_pleurostictus_nisorius",
                                   "Thryothorus_albinucha_albinucha" ,        
                                   "Thryothorus_ludovicianus" ,          
                                   "Campylorhynchus_rufinucha_capistratus",
                                  "Campylorhynchus_rufinucha_humilis",          
                                   "Campylorhynchus_turdinus_unicolor",
                                  "Campylorhynchus_turdinus_hypostictus",
                                  "Campylorhynchus_zonatus_zonatus",         
                                   "Campylorhynchus_zonatus_brevirostris", 
                                   "Polioptila_dumicola_saturata",
                                   "Polioptila_dumicola_berlepschi",
                                  "Polioptila_dumicola_dumicola",            
                                   "Polioptila_plumbea_parvirostris" ,
                                  "Polioptila_plumbea_atricapilla",
                                  "Polioptila_plumbea_plumbea",              
                                    "Polioptila_caerulea_amoenissima",         
                                   "Polioptila_caerulea_caerulea",
                                  "Polioptila_caerulea_deppei",              
                                    "Polioptila_caerulea_mexicana" ,
                                  "Polioptila_caerulea_cozumelae" ,          
                                   "Polioptila_californica_margaritae",
                                  "Polioptila_caerulea_nelsoni",             
                                   "Polioptila_caerulea_obscura" ,
                                  "Polioptila_bilineata_superciliaris",      
                                  "Polioptila_bilineata_brodkorbi" ,         
                                   "Polioptila_bilineata_cinericia",
                                 "Polioptila_plumbea_plumbiceps",
                                 "Polioptila_plumbea_innotata",             
                                   "Polioptila_nigriceps_restricta",
                                 "Polioptila_nigriceps_nigriceps",          
                                   "Polioptila_melanura_curtata" ,
                                   "Polioptila_melanura_melanura" ,
                                   "Polioptila_melanura_lucida",
                                   "Polioptila_californica_californica" ,     
                                   "Polioptila_californica_pontilis", 
                                    "Polioptila_albiloris_vanrossemi" ,        
                                   "Polioptila_albiloris_albiloris",
                                 "Polioptila_albiloris_albiventris",        
                                    "Ramphocaenus_melanurus_trinitatis" ,
                                 "Ramphocaenus_melanurus_albiventris",      
                                   "Ramphocaenus_rufiventris_ardeleo" ,       
                                   "Ramphocaenus_melanurus_melanurus" ,              
                                   "Ramphocaenus_melanurus_amazonum" ,
                                 "Ramphocaenus_melanurus_badius" ,          
                                   "Ramphocaenus_melanurus_obscurus",
                                 "Ramphocaenus_melanurus_sticturus" ,       
                                   "Microbates_cinereiventris_cinereiventris",
                                 "Microbates_cinereiventris_semitorquatus" ,
                                   "Microbates_collaris_collaris" ,
                                 "Microbates_collaris_perlatus" ,           
                                   "Microbates_cinereiventris_peruvianus",
                                  "Sitta_oenochlamys_apo",                   
                                  "Sitta_oenochlamys_isarog",
                                 "Sitta_oenochlamys_mesoleuca",            
                                 "Sitta_tephronota_iranica" ,
                                 "Sitta_europaea_sinensis",                 
                                    "Sitta_europaea_amurensis",
                                 "Campylorhynchus_zonatus_vulcanius"))

for (tip in SitCert$tip.label) {
  if ((lengths(strsplit(tip, split="_"))>2)) {
    SitCert$tip.label[SitCert$tip.label==tip] <- sub("_[^_]+$", "", tip)
  }
}

# Check for duplicated tip labels
anyDuplicated(SitCert$tip.label)                             


### CALIBRATE with backbone tree date ###

CladeAge <- getMRCAage(PasBackTre, tip=c("Certhia_familiaris", "Polioptila_caerulea")) 
#29.89933

CALIB <- makeChronosCalib(SitCert, node="root", age.min=CladeAge)

SitCertC <- chronos(SitCert, model="clock", calibration=CALIB) # clock works again now...

plot(SitCertC, cex=0.5); axisPhylo()

# make clade of 5 into 2 tip
PasBackTre31 <- drop.tip(PasBackTre30, tip=c("Tichodroma_muraria", "Sitta_carolinensis", "Troglodytes_aedon"))
# remove 2 tips from backbone
PasBackTre31 <- drop.tip(PasBackTre31, tip=c("Certhia_familiaris", "Polioptila_caerulea"),  trim.internal=FALSE)

# Graft onto Backbone
PasBackTre31 <- bind.tree(PasBackTre31, SitCertC, where=Ntip(PasBackTre31))

# check for duplicate
anyDuplicated(PasBackTre31$tip.label)


plot(PasBackTre31, show.tip.label = F)


### SAVE work###
write.tree(PasBackTre31, file="PasBackTre31.tre")








##### Cisticolidae #####
Cisti <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Neomixis_tenella", "Pycnonotus_barbatus")))
plot(Cisti)

#Pycnonotus_barbatus

CistiOnly <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Neomixis_tenella", "Eremomela_pusilla")))

tips.to.drop <- setdiff(Cisti$tip.label, CistiOnly$tip.label)
tips.to.drop <- tips.to.drop[! tips.to.drop %in% c("Pycnonotus_barbatus")]

Cisti <- drop.tip(Cisti, tip=tips.to.drop)


# collapse spp
Cisti <- drop.tip(Cisti, tip=c("Poliolais_lopesi_manengubae",
                               "Camaroptera_brachyura_brachyura",
                               "Camaroptera_brevicaudata_brevicaudata",
                               "Camaroptera_brevicaudata_tincta",
                               "Camaroptera_brachyura_pileata",
                               "Camaroptera_chloronota_chloronota" ,
                               "Apalis_melanocephala_melanocephala",
                               "Apalis_porphyrolaema_porphyrolaema",
                               "Apalis_cinerea_cinerea",
                               "Apalis_flavida_caniceps",
                               "Apalis_flavida_florisuga",
                               "Apalis_thoracica_griseiceps",
                               "Apalis_thoracica_flaviventris",
                               "Artisornis_metopias_metopias",
                               "Oreolais_pulcher_pulcher",
                               "Schistolais_leucopogon_reichenowi",
                               "Bathmocercus_rufus_vulpinus",
                               "Cisticola_juncidis_brunniceps",
                               "Cisticola_juncidis_tinnabulans",
                               "Cisticola_juncidis_juncidis",
                               "Cisticola_tinniens_perpullus" ,
                               "Cisticola_brachypterus_brachypterus",
                               "Cisticola_cantans_belli",
                               "Prinia_crinigera_parumstriata",
                               "Prinia_polychroa_polychroa",
                               "Prinia_crinigera_crinigera",
                               "Prinia_crinigera_striatula",
                               "Prinia_crinigera_catharia" ,
                               "Prinia_polychroa_bangsi",
                               "Prinia_erythroptera_erythroptera",
                               "Prinia_subflava_subflava",             
                                "Prinia_subflava_mutatrix",
                               "Prinia_bairdii_obscura",
                               "Prinia_bairdii_bairdii",
                               "Orthotomus_sutorius_maculicollis" ,    
                               "Orthotomus_sutorius_edela",
                               "Orthotomus_sutorius_longicauda",
                               "Orthotomus_sutorius_inexpectatus",
                               "Orthotomus_sericeus_sericeus",
                               "Orthotomus_sericeus_hesperius",
                               "Orthotomus_derbianus_nilesi",
                               "Orthotomus_castaneiceps_rabori",
                               "Orthotomus_cinereiceps_obscurior",
                               "Orthotomus_atrogularis_nitidus",
                               "Orthotomus_sepium_sepium",             
                               "Orthotomus_ruficeps_borneoensis" ,
                               "Orthotomus_ruficeps_cineraceus",
                               "Orthotomus_atrogularis_atrogularis"
                               ))




for (tip in Cisti$tip.label) {
  if ((lengths(strsplit(tip, split="_"))>2)) {
    Cisti$tip.label[Cisti$tip.label==tip] <- sub("_[^_]+$", "", tip)
  }
}

# Check for duplicated tip labels
anyDuplicated(Cisti$tip.label)

### CALIBRATE with backbone tree date ###

CladeAge <- getMRCAage(PasBackTre, tip=c("Cisticola_anonymus", "Pycnonotus_barbatus")) # 
#26.00996

CALIB <- makeChronosCalib(Cisti, node="root", age.min=CladeAge)

CistiC <- chronos(Cisti, model="clock", calibration=CALIB)

plot(CistiC, cex=0.5); axisPhylo()

# only crown  age
CrownAge <- getMRCAage(CistiC, tip=c("Neomixis_tenella", "Eremomela_pusilla"))

# get rid of outgroup 
CistiC <- drop.tip(CistiC, tip=c("Pycnonotus_barbatus"))

plot(CistiC, cex=0.3, root.edge=TRUE); axisPhylo()


# Bind tree but specifying the position of the new node
PasBackTre32 <- bind.tree(PasBackTre31, CistiC, where=which(PasBackTre31$tip.label =="Cisticola_anonymus"), position=CrownAge)
plot(PasBackTre32, show.tip.label = F)

# check for duplicate
anyDuplicated(PasBackTre32$tip.label)

# removing duplicates
PasBackTre32 <- drop.tip(PasBackTre32, tip=4881) # Cisticola anonymus

### SAVE work ###
write.tree(PasBackTre32, file="PasBackTre32.tre")






##### Cinclidae #####
Cincl <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Cinclus_cinclus", "Catharus_ustulatus")))
plot(Cincl)

#Pycnonotus_barbatus

CinclOnly <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Cinclus_cinclus", "Cinclus_mexicanus")))

tips.to.drop <- setdiff(Cincl$tip.label, CinclOnly$tip.label)
tips.to.drop <- tips.to.drop[! tips.to.drop %in% c("Catharus_ustulatus")]

Cincl <- drop.tip(Cincl, tip=tips.to.drop)



### CALIBRATE with backbone tree date ###

CladeAge <- getMRCAage(PasBackTre, tip=c("Cinclus_cinclus", "Catharus_ustulatus")) # 
#27.61293

CALIB <- makeChronosCalib(Cincl, node="root", age.min=CladeAge)

CinclC <- chronos(Cincl, model="clock", calibration=CALIB)


# only crown  age
CrownAge <- getMRCAage(CinclC, tip=c("Cinclus_cinclus", "Cinclus_mexicanus"))

# get rid of outgroup 
CinclC <- drop.tip(CinclC, tip=c("Catharus_ustulatus"))


# Bind tree but specifying the position of the new node
PasBackTre33 <- bind.tree(PasBackTre32, CinclC, where=which(PasBackTre32$tip.label =="Cinclus_cinclus"), position=CrownAge)
plot(PasBackTre33, show.tip.label = F)

# check for duplicate
anyDuplicated(PasBackTre33$tip.label)

# removing duplicates
PasBackTre33 <- drop.tip(PasBackTre33, tip=2544) # Cinclus cinclus

### SAVE work ###
write.tree(PasBackTre33, file="PasBackTre33.tre")







#### Pnoepygidae #####
Pnoe <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Pnoepyga_pusilla", "Cisticola_anonymus")))
plot(Pnoe)

#Pycnonotus_barbatus

PnoeOnly <- extract.clade(PhlawdTree, node=getMRCA(PhlawdTree, tip=c("Pnoepyga_pusilla", "Pnoepyga_formosana")))

tips.to.drop <- setdiff(Pnoe$tip.label, PnoeOnly$tip.label)
tips.to.drop <- tips.to.drop[! tips.to.drop %in% c("Cisticola_anonymus")]

Pnoe <- drop.tip(Pnoe, tip=tips.to.drop)



### CALIBRATE with backbone tree date ###

CladeAge <- getMRCAage(PasBackTre, tip=c("Pnoepyga_pusilla", "Cisticola_anonymus")) # 
#26.58905

CALIB <- makeChronosCalib(Pnoe, node="root", age.min=CladeAge)

PnoeC <- chronos(Pnoe, model="clock", calibration=CALIB)


# only crown  age
CrownAge <- getMRCAage(PnoeC, tip=c("Pnoepyga_pusilla", "Pnoepyga_formosana"))

# get rid of outgroup 
PnoeC <- drop.tip(PnoeC, tip=c("Cisticola_anonymus"))


# Bind tree but specifying the position of the new node
PasBackTre34 <- bind.tree(PasBackTre33, PnoeC, where=which(PasBackTre33$tip.label =="Pnoepyga_pusilla"), position=CrownAge)
plot(PasBackTre34, show.tip.label = F)

# check for duplicate
anyDuplicated(PasBackTre34$tip.label)

# removing duplicates
PasBackTre34 <- drop.tip(PasBackTre34, tip=5004) # Pnoepyga pusilla

### SAVE work ###
write.tree(PasBackTre34, file="PasBackTre34.tre")






####### Last edits for Santiago ####
# remove duplicated hylia
PasBackTre35 <- drop.tip(PasBackTre34, tip="Hylia_prasina_prasina")

#Carpodacus sibiricus duplicated as Uranus sibiricus: delete Uranus sibiricus
PasBackTre35 <- drop.tip(PasBackTre35, tip="Uragus_sibiricus")

FinalSpp<-as.data.frame((PasBackTre35$tip.label))

### SAVE work ###
write.tree(PasBackTre35, file="PasBackTre35.tre")


##### LAST LAST edits for Santiago

# Emberizoidea fixed, see above

# Chloropsis_sonnerati_zosterops
PasBackTre37 <- drop.tip(PasBackTre36, tip="Chloropsis_sonnerati_zosterops")

# delete one of the Robsonius using drop.tip
PasBackTre37 <- drop.tip(PasBackTre37, tip="Robsonius_rabori") # rip

write.tree(PasBackTre37, file="PasBackTre37.tre")



# Locustellidae
# root is splitting the genus Robsnius, that should be monophyletic and not that divergent. Revise that rooting. 
# checked, it is a polytomy in the original phlawd tree

# Emberizoidea
#Check attachments of Emberizoidea, maybe because the inclusion of Oreomystis, which is outgroup, not crown lineage, the crown became much younger than expected.
# i checked and it seems to be fine and not include Oreomystis




# ###### Cnemophilidae, Callaeidae, Notiomystidae #####
# don't bother





