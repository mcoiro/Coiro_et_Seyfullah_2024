##Loading the package
library(Claddis)
library(dispRity)

## The tree
tree <- read.nexus("data/cond.tre")

## The morphological matrix
morpho_matrix <- read_nexus_matrix("data/Matrix_20230214", equalize_weights = TRUE)

## cleaning the tree to only include taxa present in the matrix and to remove polytomies and zero-length branches
tree <- drop.tip(tree, setdiff(tree$tip.label, rownames(morpho_matrix$matrix_1$matrix)))
tree=multi2di(tree)
tree$edge.length[tree$edge.length==0]<-0.001

## Adding the tree root time
tree$root.time <- max(tree.age(tree)$age)

## Time slices
time_slices <- c(330.36,299,252,201,145,100,66,23,0)

#Rate analysis comparing time bins
rates = test_rates(
  tree,
  morpho_matrix,
  time_slices,
  branch_partitions = NULL,
  character_partitions = NULL,
  clade_partitions = NULL,
  time_partitions = partition_time_bins(8, partition_sizes_to_include = "all"),
  change_times = "random",
  test_type = "aic",
  alpha = 0.01,
  multiple_comparison_correction = "benjaminihochberg",
  polymorphism_state = "missing",
  uncertainty_state = "missing",
  inapplicable_state = "missing",
  time_binning_approach = "lloyd",
  all_weights_integers = FALSE,
  estimate_all_nodes = FALSE,
  estimate_tip_values = FALSE,
  inapplicables_as_missing = FALSE,
  polymorphism_behaviour = "equalp",
  uncertainty_behaviour = "equalp",
  threshold = 0.01,
  all_missing_allowed = FALSE
)

## obtaining the AICc values from the rates object
aics = list()
for (i in 1:120){
  aics[[i]] =  rates$time_test_results[[i]]$aicc
}
aics_l =unlist(aics)

## Finding the best model
which(aics_l == min(aics_l))

## Plotting the best model
plot_rates_time(
  test_rates_output = rates,
  model_number = 81
)

#####Nitrogen fixation rate analysis

#Set up the taxon set for Zamiaceae
zam_node = getMRCA(tree, c("Dioon_edule", "Zamia_erosa"))
zam_tax = na.omit(tree$tip.label[getDescendants(tree,zam_node )])

#Set up the taxon sets for Cycadaceae
cyc_node = getMRCA(tree, c("Cycas_revoluta", "Cycas_sp_Sakhalin"))
cyc_tax = na.omit(tree$tip.label[getDescendants(tree,cyc_node )])

nonex_tax = setdiff(tree$tip.label, c(zam_tax,cyc_tax))

#separating the taxa set between nitrogen fixating and not nitrogen fixating
nitro= c(zam_tax, cyc_tax)
nonnitro= setdiff(tree$tip.label, c(zam_tax,cyc_tax))

 
 rates_nitros = test_rates(
   tree,
   morpho_matrix,
   time_slices,
   branch_partitions = NULL,
   character_partitions = NULL,
   clade_partitions = c(zam_node,cyc_node, 170),
   time_partitions = NULL,
   change_times = "random",
   test_type = "aic",
   alpha = 0.01,
   multiple_comparison_correction = "benjaminihochberg",
   polymorphism_state = "missing",
   uncertainty_state = "missing",
   inapplicable_state = "missing",
   time_binning_approach = "lloyd",
   all_weights_integers = FALSE,
   estimate_all_nodes = FALSE,
   estimate_tip_values = FALSE,
   inapplicables_as_missing = FALSE,
   polymorphism_behaviour = "equalp",
   uncertainty_behaviour = "equalp",
   threshold = 0.01,
   all_missing_allowed = FALSE
 )
 
 
 ## obtaining the AICc values from the rates object
 aicsn = list()
 for (i in 1:3){
   aicsn[[i]] =  rates_nitros$clade_test_results[[i]]$aicc
 }
 aics_ln =unlist(aicsn)
 
 ## Finding the best model
 which(aics_ln == min(aics_ln))
 
 #############
 # Figure X 
 pdf(file = "FigX.pdf",   # The directory you want to save the file in
     width = 12, # The width of the plot in inches
     height = 12) # The height of the plot in inches
 ##Graphical parameters
 op <- par(bty = "n", mfrow = c(1,3))
 
 ## plotting model 3
 plot_rates_tree(rates_nitros, model_type = "clade", model_number = 3)
  ## plotting model 1
 plot_rates_tree(rates_nitros, model_type = "clade", model_number = 1)
  ## plotting model 2
 plot_rates_tree(rates_nitros, model_type = "clade", model_number = 2)
  ## Resetting graphical parameters
 par(op)
 # Run dev.off() to create the file!
 dev.off()
 