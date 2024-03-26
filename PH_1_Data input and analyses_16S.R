#############################################
# International core Microbiome - Genotypes #
#############################################

###############
# Philippines #
###############

### Contact: Paul Dennis (p.dennis@uq.edu.au) and Belle Clarke, 2024.

#################################
# Data input and analyses - 16S #
#################################

source('Functions.R')

### Load the OTU tables into the environment
otu.PH.tmp <- read.table("../../Data/16S/otu.16S.2400.PH.csv", header=TRUE, sep=',', row.names = 1)
otu.PH <- otu.PH.tmp[, apply(otu.PH.tmp, 2, sum) > 0]/2400 # remove OTUs with no seqs and divide by 2400 to get relative abundance (0-1)

### Load the taxonomy information into the environment
tax.PH <- read.table("../../Data/16S/taxonomy.16S.PH.csv", header=TRUE, sep=',', row.names = 1)

### Load the weighted UniFrac distance matrices into the environment
w.unifrac.PH.tmp <- read.table("../../Data/16S/w.unifrac.16S.PH.csv", header=TRUE, sep=',', row.names = 1)
w.unifrac.PH <- as.matrix(w.unifrac.PH.tmp)

### Load the environmental metadata tables with alpha diversity into the environment and set factors
env.PH <- read.table("../../Data/16S/env.16S.2400.PH.csv", header=TRUE, sep=',', row.names = 1)
env.PH$Compartment <- factor(env.PH$Compartment)
env.PH$Banana <- factor(env.PH$Banana)

### Check samples are in the same order
row.names(otu.PH) == row.names(env.PH)
row.names(w.unifrac.PH) == row.names(env.PH)

# Color palettes to represent the compartments and bananas
mypal.compartment <- colorRampPalette(c('#8c510a','#f6e8c3'))
mypal.banana <- colorRampPalette(c('#c7e9c0','#74c476','#238b45'))
mysymbol.compartment <- c(21,23)
mysymbol.banana <- c(21,23,22)

#################################
## Analyses of alpha diversity ##
#################################

alpha.indicators <- c('Sobs', 'Chao1', 'Shan', 'PD')

# ANOVA with posthoc
for(i in alpha.indicators) {
  print(i)
  print(summary(aov(env.PH[,i] ~ Compartment * Banana, data = env.PH)))
  
  print(cld(lsmeans(aov(env.PH[,i] ~ Compartment * Banana, data = env.PH), 
                    ~ pairwise ~ Compartment * Banana), 
            Letters = letters, 
            adjust = "bh"))
}

### Visualise the data using a barplot with Standard Errors of the means as error bars   
bargraph.CI(Compartment, Sobs, Banana, data = env.PH, col = mypal.banana(3), legend = T) 
bargraph.CI(Compartment, Chao1, Banana, data = env.PH, col = mypal.banana(3), legend = T) 
bargraph.CI(Compartment, Shan, Banana, data = env.PH, col = mypal.banana(3), legend = T) 
bargraph.CI(Compartment, PD, Banana, data = env.PH, col = mypal.banana(3), legend = T) 

### Visualise the data using a barplot with Standard Deviation as error bars   
#bargraph.CI(Treatment.full, Shan, data = env, col = mypal.treatment(7),
#            ci.fun= function(x) c(mean(x)-sd(x), mean(x) + sd(x))) 

##################################
### Analyses of Beta diversity ###
##################################

# Using square root transformed OTU relative abundances (aka Hellinger transformation)

###############################################################################
# Determine whether the composition of bacterial communities is significantly #
#       influenced by Compartment and/or banana using PERMANOVA               #
###############################################################################
set.seed(1)
# Hellinger transformed OTU abundances
adonis2(sqrt(otu.PH) ~ Compartment * Banana, data = env.PH, method = 'euc')
# Weighted UniFrac distances
#adonis2(w.unifrac.PH ~ Compartment * Banana, data = env.PH, method = 'euc')

############################################################
# View community similar/dissimilarity in ordination space # 
############################################################

### First we'll make a distance-based RDA ordination object.
otu.PH.pca <- rda(sqrt(otu.PH)) # PCA gives a euclidean projection of the input, which in dbPCA is transformed (e.g. Hellinger)
axis.percent(otu.PH.pca) # Returns the proportion (%) of total inertia (compositional dissimilarity) explained on each axis 
#otu.PH.rda <- rda(sqrt(otu.PH) ~ Compartment, data = env.PH) # RDA gives a euclidean projection of the input constrained by Compartment (PERMANOVA showed it was significant), which in dbPCA is transformed (e.g. Hellinger)
#axis.percent(otu.PH.rda) 

# A very basic plot can be obtained as follows:
plot(otu.PH.pca, type = 'n', main = "BS (Circles) ER (Diamonds)", scaling = 3)
points(otu.PH.pca, dis='sites', 
       pch = mysymbol.compartment[env.PH$Compartment], 
       bg = mypal.banana(3)[env.PH$Banana], scaling = 3, cex = 2)
legend("bottomright", legend=unique(env.PH$Banana), pch=19,
       col=mypal.banana(3)[unique((env.PH$Banana))])

### Here is a much more complex but informative plot: 

### Start 
# set the constant variables used below then plot the axes with %variation explained
ord = otu.PH.pca
scaling.val = 3

plot(ord, main = "BS (Circles) ER (Diamonds)",
     type='n', scaling=scaling.val, 
     xlab=paste("db-PC1 (",axis.percent(
       ord)[[1]],"%)",sep=""),
     ylab=paste("db-PC2 (",axis.percent(
       ord)[[2]],"%)",sep=""))

# Add ellipses representing the SD for each treatment combo
ordiellipse(ord, env.PH$Compartment:env.PH$Banana, 
            col=mypal.banana(3), 
            draw = 'polygon', alpha = 0.5, lty=3,  kind ='sd',scaling=scaling.val) # can also use kind = 'se' etc...

# Plot the OTUs in grey - we can't read most of them anyway (this can also be omitted)
points(ord, dis='sp', pch=4, col='grey', cex=0.6, scaling=scaling.val)

# Plot the Samples colored by Banana and shaped by Compartment 
points(ord, dis='sites', pch=mysymbol.compartment[env.PH$Compartment], scaling = scaling.val,
       bg = addTrans(mypal.banana(3)[factor(env.PH$Banana)], 220), cex = 2.5)

# The next block of code simply highlights which OTUs are most discriminating in red and adds their OTU_ID codes      
sd.val = 8 # This is the number of standard deviations from the mean position of OTUs on each axis
points(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,1] > sd.val * sd(scores(ord)$sp[,1])),],
       pch=4,col='darkred',cex=0.6)	
points(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,1] < 0 - (sd.val * sd(scores(ord)$sp[,1]))),],
       pch=4,col='darkred',cex=0.6) 
points(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,2] > sd.val * sd(scores(ord)$sp[,2])),],
       pch=4,col='darkred',cex=0.6)
points(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,2] < 0 - (sd.val * sd(scores(ord)$sp[,2]))),],
       pch=4,col='darkred',cex=0.6) 
orditorp(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,1] > sd.val * sd(scores(ord)$sp[,1])),], 
         "sp", pch="+", col="black", air = 0.8, scaling = scaling.val)
orditorp(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,1] < 0 - (sd.val * sd(scores(ord)$sp[,1]))),], 
         "sp", pch="+", col="black", air = 0.8, scaling = scaling.val)
orditorp(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,2] > sd.val * sd(scores(ord)$sp[,2])),], 
         "sp", pch="+", col="black", air = 0.8, scaling = scaling.val)
orditorp(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,2] < 0 - (sd.val * sd(scores(ord)$sp[,2]))),], 
         "sp", pch="+", col="black", air = 0.8, scaling = scaling.val)

# This adds the legend
legend("bottomright",legend=unique(env.PH$Banana), pch=19,
       col=mypal.banana(3)[unique(factor(env.PH$Banana))])

### END

#############################################################
# View the relative abundance of dominant OTUs in a heatmap #
#############################################################

# Select OTUs based on a relative abundance threshold 
hm.tmp <- otu.PH[,which(apply(otu.PH,2,max)>= 0.01)]
dim(hm.tmp) # too many OTUs? Could select those that are on average above threshold within treatment combos

# limit to OTUs that have a mean relative abundance of >=1% within treatment
otu.PH.mean <- aggregate(otu.PH,list(env.PH$Compartment:env.PH$Banana), mean)
row.names(otu.PH.mean) <- otu.PH.mean[,1]
otu.PH.mean <- otu.PH.mean[,-1]

hm.tmp <- otu.PH.mean[,which(apply(otu.PH.mean,2,max)>= 0.01)]
dim(hm.tmp) # Too few? Could drop the RA threshold 

# Get the taxonomy information for the OTUs in the heatmap
otus.for.heatmap.tmp <- tax.PH[colnames(hm.tmp),]
otus.for.heatmap <- otus.for.heatmap.tmp[order(otus.for.heatmap.tmp$Silva138_taxonomy),]

hm.mean <- hm.tmp[,row.names(otus.for.heatmap)]

hm.all <- otu.PH[,row.names(otus.for.heatmap)]

# Make a full.id label by concatenated the OTU_ID (in square brackets) with the taxonomy
otus.for.heatmap$full.id <- paste("[", 
                                  otus.for.heatmap$OTU, 
                                  "] ", 
                                  otus.for.heatmap$Silva138_taxonomy, 
                                  sep='')

# Make a heatmap
mypal.bw <- colorRampPalette(c("White","Black"))
heatmap(t(sqrt(hm.all)), revC = TRUE, 
        Colv = NA, Rowv = NA, col=mypal.bw(180),
        labRow = otus.for.heatmap$full.id, 
        labCol = env.PH$Compartment:env.PH$Banana, 
        scale = 'none', margins = c(5,15))


for(i in otus.for.heatmap$OTU){
  print(i)
  print(anova(lm(otu.PH[,i] ~ env.PH$Compartment * env.PH$Banana)))
  bargraph.CI(x.factor = env.PH$Compartment, response = otu.PH[,i], 
              group = env.PH$Banana, col = mypal.banana(3), legend = T, main = i)
}



