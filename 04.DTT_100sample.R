library(phytools)
library(Claddis)
#devtools::install_github("TGuillerme/dispRity", force=TRUE)
library(dispRity)


## The tree
trees <- read.nexus("data/sample100")



## The morphological matrix
morpho_matrix <- read_nexus_matrix("data/Matrix_20230214", equalize_weights = TRUE)

## Cleaning the trees 
treesrem = clean.data(morpho_matrix$matrix_1$matrix,trees)
trees= treesrem$tree
trees = remove.zero.brlen(trees, 0.001)

## adding node labels and root times
for (i in 1:length(trees)){
  trees[[i]]$node.label <- paste0(1:Nnode(trees[[i]]))
  trees[[i]]$root.time <- max(tree.age(trees[[i]])$age)
}

## Ancestral reconstruction across the trees
multimorph = multi.ace(morpho_matrix$matrix_1$matrix,trees ,threshold=0.99, parallel=TRUE,output = "combined.matrix")



## transforming in distances
distances <- lapply(multimorph, char.diff, method = "mord", by.col = FALSE)


## Running PCoA on the distance matrices
#pcoas = lapply(distances, cmdscale, k=100, add=TRUE)
pcoas = lapply(distances, ape::pcoa, correction="cailliez")
pcoas_list= list()

for (i in 1:length(pcoas)){
pcoas_list[[i]] = pcoas[[i]]$vectors.cor 
  
}
## 


#trees[[i]]$tip.label <- paste0(trees[[i]]$tip.label,i)


## Adding rownames to the pcoas matrices ???
for (i in 1:length(pcoas_list)){
  rownames(pcoas_list[[i]]) = c(trees[[i]]$tip.label, trees[[i]]$node.label)
}


time_slices <- seq(from = 0, to = 300, by = 10)

## Gradual splits model
gradual_model <- chrono.subsets(pcoas_list, trees, time = time_slices, method = "continuous",
                                model = "gradual.split", inc.nodes=TRUE, bind.data = TRUE)

## Bootstrapping both models (500 times)
gradual_model_nr <- boot.matrix(gradual_model, 14500)


gradual_sum_var <- dispRity(gradual_model_nr, metric = c(sum, variances))
gradual_centr0_dist <- dispRity(gradual_model_nr, metric = c(mean, centroids),centroid=0,  method = "manhattan")





## Graphical parameters
#############
# Figure S4 
pdf(file = paste0("output/sample100.pdf"),   # The directory you want to save the file in
    width = 11.69, # The width of the plot in inches
    height = 8.27) # The height of the plot in inches
op <- par(bty = "n", mfrow = c(1,2))
plot(gradual_sum_var, xlab = "")
plot(gradual_centr0_dist, xlab = "")
par(op)
# Run dev.off() to create the file!
dev.off()

