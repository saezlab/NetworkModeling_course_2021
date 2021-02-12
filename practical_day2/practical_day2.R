#Load the libraries
library(viper)
library(fgsea)
library(cosmos)

###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################

#Set working directory
setwd("~/Dropbox/Mobi_course_2021/practical/")

###############################################################################################################################################
###############################################################################################################################################
#Prepare input of kinase activity estimation

#load post-translational modification database
load("data/omnipath_ptm.RData")

#load phospho differential anaylsis
load("data/phospho_ttop.RData")

#filter and process the ptm database
omnipath_ptm <- omnipath_ptm[omnipath_ptm$modification %in% c("dephosphorylation","phosphorylation"),]
KSN <- omnipath_ptm[,c(4,3)]
KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol,omnipath_ptm$residue_type, sep ="_")
KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol,omnipath_ptm$residue_offset, sep = "")
KSN$sign <- ifelse(omnipath_ptm$modification == "phosphorylation", 1, -1)

###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################
#Estimate kianse activities

#load the ptm dabase in viper input format
load("data/omnipath_ptm_viper.RData")

#load the phospho differential anaylsis result in viper input format
load("data/viper_phospho.RData")

###Compute TF activity scores
kin_activity <- as.data.frame(viper(eset = viper_phospho, regulon = KSN_viper, minsize = 5, adaptive.size = F, eset.filter = F))

names(kin_activity) <- "NES"
kin_activity$kinase <- row.names(kin_activity)

###############################################################################################################################################
###############################################################################################################################################
##GSEA on kinactivity

#load the kianse activity estiamtion result
load("data/kin_activity.RData")

#load a pathway ontology
load("data/pathway_list.RData")

#Filter for a given type of pathways
pathways_list <- pathways_list[grepl("KEGG",names(pathways_list))]

#prepare kianse activities as gsea input
kin_activity_vector <- kin_activity$NES
names(kin_activity_vector) <- kin_activity$kinase

#run gsea on kinase activities
fgsea_res <- fgsea(pathways = pathways_list, stats = kin_activity_vector, minSize = 5, maxSize = 100, gseaParam = 1)

#load the fgsea result
load("data/fgsea_res.RData")

###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################

#In order to adapt options to users specification we can load them into a variable 
#that will then be passed to preprocess_COSMOS_signaling_to_metabolism CARNIVAL_options parameter
my_options <- default_CARNIVAL_options()

#Here the user should provide a path to its CPLEX executable (only cplex at the moment, other solvers will be documented soon !)
my_options$solverPath <- "~/Documents/cplex"

#### FORWARD run of COSMOS, to connect signaling to metabolism

#The signaling inputs are the result of footprint based TF and kinase activity estiamtion
#For more info on TF activity estiamtion from transcriptomic data, see:https://github.com/saezlab/transcriptutorial (Especially chapter 4)

#Here we use of toy PKN, to see the full meta PKN, you can load it with load_meta_pkn()

#The metabolites in the prior knowledge network are identified as XMetab__PUBCHEMid___compartment____ or XMetab__BIGGid___compartment____
#for example “XMetab__6804___m____”. The compartment code is the BIGG model standard (r, c, e, x, m, l, n, g). 
#Thus we will first need to map whatever identifer for metabolite the data has to the one of the network.
#Genes are identified as XENTREZid (in the signaling part of network) or XGene####__ENTREZid (in the reaction network part of network)

load("data/toy_sif.RData")

#prepare the input for comos forward run
test_for <- preprocess_COSMOS_signaling_to_metabolism(meta_network = toy_sif,
                                                      signaling_data = toy_signaling_input_carnival_vec,
                                                      metabolic_data = toy_metab_input_carnival_vec,
                                                      diff_expression_data = toy_RNA,
                                                      maximum_network_depth = 15,
                                                      remove_unexpressed_nodes = T,
                                                      CARNIVAL_options = my_options
                                                      
)

load("data/cosmos_test_for.RData")

#run cosmos forward
test_result_for <- run_COSMOS_signaling_to_metabolism(data = test_for,
                                                      CARNIVAL_options = my_options)

load("data/cosmos_test_result_for.RData")

metab_to_pubchem_vec <- metab_to_pubchem$name
names(metab_to_pubchem_vec) <- metab_to_pubchem$pubchem

test_result_for <- format_COSMOS_res(test_result_for,
                                     metab_mapping = metab_to_pubchem_vec,
                                     measured_nodes = unique(c(names(toy_metab_input_carnival_vec),names(toy_signaling_input_carnival_vec))),
                                     omnipath_ptm = omnipath_ptm)

load("data/cosmos_test_result_for_processed.RData")

View(test_result_for[[1]]) #SIF
View(test_result_for[[2]]) #ATTRIBUTES

#### BACKWARD run of COSMOS, to connect metabolism to signaling
test_back <- preprocess_COSMOS_metabolism_to_signaling(meta_network = toy_sif,
                                                       signaling_data = toy_signaling_input_carnival_vec,
                                                       metabolic_data = toy_metab_input_carnival_vec,
                                                       diff_expression_data = toy_RNA,
                                                       maximum_network_depth = 15,
                                                       remove_unexpressed_nodes = F,
                                                       CARNIVAL_options = my_options
                                                       
)

load("data/cosmos_test_back.RData")

test_result_back <- run_COSMOS_metabolism_to_signaling(data = test_back,
                                                       CARNIVAL_options = my_options)

load("data/cosmos_test_result_back.RData")

test_result_back <- format_COSMOS_res(test_result_back,
                                      metab_mapping = metab_to_pubchem_vec,
                                      measured_nodes = unique(c(names(toy_metab_input_carnival_vec),names(toy_signaling_input_carnival_vec))),
                                      omnipath_ptm = omnipath_ptm)

load("data/cosmos_test_result_back_processed.RData")

View(test_result_back[[1]]) #SIF
View(test_result_back[[2]]) #ATTRIBUTES

###Merge forward and backward networks

full_sif <- as.data.frame(rbind(test_result_for[[1]], test_result_back[[1]]))
full_attributes <- as.data.frame(rbind(test_result_for[[2]], test_result_back[[2]]))

full_sif <- unique(full_sif)
full_attributes <- unique(full_attributes)

###Visualise network

#you can give here a signle nodes name or a vector or names. 
#This functions will allow to display the neighboorhood of given nodes

network_plot <- display_node_neighboorhood(central_node = 'PRKACA', sif = full_sif, att = full_attributes, n = 5)

load("data/cosmos_network_PRKACA.RData")

network_plot
