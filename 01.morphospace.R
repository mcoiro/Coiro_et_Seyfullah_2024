##Â Loading the package
library(Claddis)
library(dispRity)

## The tree
tree <- read.nexus("data/cond.tre")

## The morphological matrix
morpho_matrix <- read_nexus_matrix("data/Matrix_20230214", equalize_weights = TRUE)

## cleaning the tree to only include taxa present in the matrix and to remove polytomies and zero-length branches
tree <- drop.tip(tree, setdiff(tree$tip.label, rownames(morpho_matrix$matrix_1$matrix)))

#tree$edge.length[tree$edge.length==0]<-0.001

## Adding the tree root time
tree$root.time <- max(tree.age(tree)$age)


## Adding node labels
tree <- makeNodeLabel(tree, method = "number", prefix="")

#generate the morphospace
morphospace = ordinate_cladistic_matrix(morpho_matrix, time_tree=tree)
morphospace_ext = ordinate_cladistic_matrix(morpho_matrix, time_tree=NULL)

#Extract the vector matrix from the Claddis morphospace object to feed it to dispRity
morphosp = morphospace$vectors
rownames(morphosp) = c(tree$tip.label,tree$node.label)

morphosp_tip = morphospace_ext$vectors
rownames(morphosp_tip) = tree$tip.label

#saving the morphospaces to file
write.csv(morphosp, "output/morphospace")
write.csv(morphosp_tip, "output/morphospace_ext")


# Before we plot this lets get the data we need to make a scree plot:
scree.data <- apply(morphospace$vectors, 2, var) / sum(apply(morphospace$vectors, 2, var)) * 100
scree.data_ext <- apply(morphospace_ext$vectors, 2, var) / sum(apply(morphospace_ext$vectors, 2, var)) * 100


#############
# Figure S1 
pdf(file = "FigS1.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 12) # The height of the plot in inches
##Graphical parameters
op <- par(bty = "n", mfrow = c(2,1))
# We can make a simple plot of this:
plot(scree.data, type="l", xlab="Ordination axis", ylab="Percentage variance")
plot(scree.data_ext, type="l", xlab="Ordination axis", ylab="Percentage variance")
## Resetting graphical parameters
par(op)
# Run dev.off() to create the file!
dev.off()

##Setting up taxon groups
zam_node = getMRCA(tree, c("Dioon_edule", "Zamia_erosa"))
zam_tax = na.omit(tree$tip.label[getDescendants(tree,zam_node )])

cyc_node = getMRCA(tree, c("Cycas_revoluta", "Paracycas_raripinnata"))
cyc_tax = na.omit(tree$tip.label[getDescendants(tree,cyc_node )])

nonex_tax = setdiff(tree$tip.label, c(zam_tax,cyc_tax))
taxongroups = list(Zamiaceae = zam_tax, Cycadaceae=cyc_tax, stem=nonex_tax)


#############
# Figure S2 
pdf(file = "FigS2.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 12) # The height of the plot in inches
##Graphical parameters
op <- par(bty = "n", mfrow = c(2,1))
plot_multi_morphospace(pcoa_input= morphospace, taxon_groups = taxongroups, n_axes=8)
plot_multi_morphospace(pcoa_input= morphospace_ext, taxon_groups = taxongroups, n_axes=8)
## Resetting graphical parameters
par(op)
# Run dev.off() to create the file!
dev.off()

