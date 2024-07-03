library(phytools)
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

#the morphospaces 

morphosp =read.csv("output/morphospace", row.names = 1)
morphosp_tip = read.csv("output/morphospace_ext", row.names = 1 )


######TIP ONLY MORPHOSPACE
##########################

FADLAD = read.table('data/FADLAD.csv',  sep=',',row.names = 1, header=TRUE)

## Time slices
time_slices <- c(252,201.4,145,66,23,0)

## Gradual splits model
gradual_model_tip <- chrono.subsets(morphosp_tip,  time = time_slices,method = "discrete",FADLAD=FADLAD)

##CB Bootstrapping both models (500 times)
gradual_model_nr_tip <- boot.matrix(gradual_model_tip, 500)


##CB Bootstrapping both models (500 times) with rarefaction to 6  
gradual_model_6_tip <- boot.matrix(gradual_model_tip, 500, rarefaction =6)


##CB Bootstrapping both models (500 times) with rarefaction to 3
gradual_model_3_tip <- boot.matrix(gradual_model_tip, 500, rarefaction =3)



gradual_sum_var_tip <- dispRity(gradual_model_nr_tip, metric = c(sum, variances))
gradual_centr0_dist_tip <- dispRity(gradual_model_nr_tip, metric = c(mean, centroids),centroid=0,  method = "manhattan")


gradual_sum_var6_tip <- dispRity(gradual_model_6_tip, metric = c(sum, variances))
gradual_centr0_dist6_tip <- dispRity(gradual_model_6_tip, metric = c(mean, centroids),centroid=0,  method = "manhattan")


gradual_sum_var3_tip <- dispRity(gradual_model_3_tip, metric = c(sum, variances))
gradual_centr0_dist3_tip <- dispRity(gradual_model_3_tip, metric = c(mean, centroids),centroid=0,  method = "manhattan")




## Graphical parameters
#############
# Figure S3 
pdf(file = "output/FigS3.pdf",   # The directory you want to save the file in
    width =  8.27, # The width of the plot in inches
    height = 11.69) # The height of the plot in inches
op <- par(bty = "n", mfrow = c(3,2))
## Plotting the sum of variances
plot(gradual_sum_var_tip, xlab = "", type="continuous")
plot(gradual_centr0_dist_tip, xlab = "", type="continuous")
plot(gradual_sum_var6_tip, xlab = "", type="continuous")
plot(gradual_centr0_dist6_tip, xlab = "", type="continuous")
plot(gradual_sum_var3_tip, xlab = "", type="continuous")
plot(gradual_centr0_dist3_tip, xlab = "", type="continuous")
par(op)
# Run dev.off() to create the file!
dev.off()


######POARS analysis
####################

## Time slices
time_slices <- seq(from = 0, to = 300, by = 10)

## Gradual splits model
gradual_model <- chrono.subsets(morphosp, tree, time = time_slices, method = "continuous",
                                model = "gradual.split", inc.nodes=TRUE)
## Punctuated model
punctuated_model <- chrono.subsets(morphosp, tree, time = time_slices, method = "continuous",
                                   model = "proximity")

## Bootstrapping both models (500 times)
gradual_model_nr <- boot.matrix(gradual_model, 500)
punctuated_model_nr <- boot.matrix(punctuated_model, 500)


## Bootstrapping both models (500 times) with rarefaction to 6 
gradual_model_6 <- boot.matrix(gradual_model, 500, rarefaction =6)
punctuated_model_6 <- boot.matrix(punctuated_model, 500, rarefaction =6)


## Bootstrapping both models (500 times) with rarefaction to 3
gradual_model_3 <- boot.matrix(gradual_model, 500, rarefaction =3)
punctuated_model_3 <- boot.matrix(punctuated_model, 500, rarefaction =3)


## 
gradual_sum_var <- dispRity(gradual_model_nr, metric = c(sum, variances))
gradual_centr0_dist <- dispRity(gradual_model_nr, metric = c(mean, centroids),centroid=0,  method = "manhattan")


gradual_sum_var6 <- dispRity(gradual_model_6, metric = c(sum, variances))
gradual_centr0_dist6 <- dispRity(gradual_model_6, metric = c(mean, centroids),centroid=0,  method = "manhattan")


gradual_sum_var3 <- dispRity(gradual_model_3, metric = c(sum, variances))
gradual_centr0_dist3 <- dispRity(gradual_model_3, metric = c(mean, centroids),centroid=0,  method = "manhattan")


punctuated_sum_var_nr <- dispRity(punctuated_model_nr, metric =  c(sum, variances))
punctuated_centr0_dist_nr <- dispRity(punctuated_model_nr, metric = c(mean, centroids),centroid=0,  method = "manhattan")

punctuated_sum_var_6 <- dispRity(punctuated_model_3, metric =  c(sum, variances))
punctuated_centr0_dist_6 <- dispRity(punctuated_model_3, metric = c(mean, centroids),centroid=0,  method = "manhattan")

punctuated_sum_var_3 <- dispRity(punctuated_model_3, metric =  c(sum, variances))
punctuated_centr0_dist_3 <- dispRity(punctuated_model_3, metric = c(mean, centroids),centroid=0,  method = "manhattan")



## Graphical parameters
#############
# Figure S4 
pdf(file = "output/FigS4.pdf",   # The directory you want to save the file in
    width = 11.69, # The width of the plot in inches
    height = 8.27) # The height of the plot in inches
op <- par(bty = "n", mfrow = c(3,4))
plot(gradual_sum_var, xlab = "")
plot(punctuated_sum_var_nr, xlab = "")
plot(gradual_centr0_dist, xlab = "")
plot(punctuated_centr0_dist_nr, xlab = "")
plot(gradual_sum_var6, xlab = "")
plot(punctuated_sum_var_6, xlab = "")
plot(gradual_centr0_dist6, xlab = "")
plot(punctuated_centr0_dist_6, xlab = "")
plot(gradual_sum_var3, xlab = "")
plot(punctuated_sum_var_3, xlab = "")
plot(gradual_centr0_dist3, xlab = "")
plot(punctuated_centr0_dist_3, xlab = "")
par(op)
# Run dev.off() to create the file!
dev.off()



## Graphical parameters
#############
# Figure 1 
pdf(file = "output/Fig1.pdf",   # The directory you want to save the file in
    width = 11.69, # The width of the plot in inches
    height =  8.27) # The height of the plot in inches
op <- par(bty = "n", mfrow = c(2,2))
plot(gradual_sum_var_tip, xlab = "", type="continuous")
plot(gradual_centr0_dist_tip, xlab = "", type="continuous")
plot(gradual_sum_var, xlab = "")
plot(gradual_centr0_dist, xlab = "")
par(op)
# Run dev.off() to create the file!
dev.off()



#####ZAMIACEAE 

#Set up the taxon set for Zamiaceae
zam_node = getMRCA(tree, c("Dioon_edule", "Zamia_erosa"))
zam_tax = na.omit(tree$tip.label[getDescendants(tree,zam_node )])

#Set up the taxon sets for Cycadaceae
cyc_node = getMRCA(tree, c("Cycas_revoluta", "Paracycas_raripinnata"))
cyc_tax = na.omit(tree$tip.label[getDescendants(tree,cyc_node )])

#Set up the taxon sets for Extinct Taxa
nonex_tax = setdiff(tree$tip.label, c(zam_tax,cyc_tax))

#Trimming the tree to only include Zamiaceae
zam_tree= drop.tip(tree, c(cyc_tax,nonex_tax))

#Setting the age of the root of the Zamiaceae tree
zam_tree$root.time = 182.4865

## Time slices
time_slices <- seq(from = 0, to = 180, by = 10)

## cleaning the data
cleaned_data <- clean.data(morphosp, zam_tree, inc.nodes = TRUE)

## Gradual splits model
gradual_model_zam <- chrono.subsets(cleaned_data$data, cleaned_data$tree, time = time_slices, method = "continuous",
                                    model = "gradual.split", inc.nodes=TRUE)
## Punctuated model
punctuated_model_zam <- chrono.subsets(cleaned_data$data, cleaned_data$tre, time = time_slices, method = "continuous",
                                       model = "proximity")

##CB Bootstrapping both models (500 times)
gradual_model_zam <- boot.matrix(gradual_model_zam, 500)
punctuated_model_zam <- boot.matrix(punctuated_model_zam, 500)


gradual_sum_var_zam <- dispRity(gradual_model_zam, metric = c(sum, variances))
gradual_centr0_dist_zam <- dispRity(gradual_model_zam, metric = c(mean, centroids),centroid=0,  method = "manhattan")


punctuated_sum_var_zam <- dispRity(punctuated_model_zam, metric = c(sum, variances))
punctuated_centr0_dist_zam <- dispRity(punctuated_model_zam, metric = c(mean, centroids),centroid=0,  method = "manhattan")



## Graphical parameters
#############
# Figure S5 
pdf(file = "output/FigS5.pdf",   # The directory you want to save the file in
    width =  8.27, # The width of the plot in inches
    height = 11.69) # The height of the plot in inches
##CB Graphical parameters
op <- par(bty = "n", mfrow = c(2,2))
## Plotting the sum of variances
plot(gradual_sum_var_zam, xlab = "")
plot(gradual_centr0_dist_zam, xlab = "")
## Plotting the sum of variances
plot(punctuated_sum_var_zam, xlab = "")
plot(punctuated_centr0_dist_zam, xlab = "")
## Resetting graphical parameters
par(op)
# Run dev.off() to create the file!
dev.off()




###EXTINCT CLADE ANALYSIS


nozam_tree= drop.tip(tree,c(zam_tax,cyc_tax))

## Time slices
time_slices <- seq(from = 20, to = 300, by = 10)

## cleaning the data
cleaned_data <- clean.data(morphosp, nozam_tree, inc.nodes = TRUE)

## Gradual splits model
gradual_model_extinct <- chrono.subsets(cleaned_data$data, cleaned_data$tree, time = time_slices, method = "continuous",
                                        model = "gradual.split", inc.nodes=TRUE)
## Punctuated model
punctuated_model_extinct <- chrono.subsets(cleaned_data$data, cleaned_data$tree, time = time_slices, method = "continuous",
                                           model = "proximity")


##CB Bootstrapping both models (500 times)
gradual_model_extinct <- boot.matrix(gradual_model_extinct, 500)
punctuated_model_extinct <- boot.matrix(punctuated_model_extinct, 500)


gradual_sum_var_extinct <- dispRity(gradual_model_extinct, metric = c(sum, variances))
gradual_centr0_dist_extinct <- dispRity(gradual_model_extinct, metric = c(mean, centroids),centroid=0,  method = "manhattan")


punctuated_sum_var_extinct <- dispRity(punctuated_model_extinct, metric = c(sum, variances))
punctuated_centr0_dist_extinct <- dispRity(punctuated_model_extinct, metric = c(mean, centroids),centroid=0,  method = "manhattan")



## Graphical parameters
#############
# Figure S6 
pdf(file = "output/FigS6.pdf",   # The directory you want to save the file in
    width =  8.27, # The width of the plot in inches
    height = 11.69) # The height of the plot in inches
##CB Graphical parameters
op <- par(bty = "n", mfrow = c(2,2))
## Plotting
plot(gradual_sum_var_extinct, xlab = "")
plot(gradual_centr0_dist_extinct, xlab = "")
plot(punctuated_sum_var_extinct, xlab = "")
plot(punctuated_centr0_dist_extinct, xlab = "")
## Resetting graphical parameters
par(op)
# Run dev.off() to create the file!
dev.off()



## Graphical parameters
#############
# Figure 2
pdf(file = "output/Fig2.pdf",   # The directory you want to save the file in
    width =  8.27, # The width of the plot in inches
    height = 11.69) # The height of the plot in inches
##CB Graphical parameters
op <- par(bty = "n", mfrow = c(2,1))
## Plotting the sum of variances
plot(gradual_sum_var_zam, xlab = "")
plot(gradual_sum_var_extinct, xlab = "")
par(op)
# Run dev.off() to create the file!
dev.off()

