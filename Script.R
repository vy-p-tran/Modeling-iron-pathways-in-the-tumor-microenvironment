#Compiled by Vy Tran, email: vtran21@jhmi.edu

#Most codes were modifed from ARCHS4 https://amp.pharm.mssm.edu/archs4/help.html (Ma'ayan et al.), and from the Weight Gene 
#Correlation Network Analysis (WGCNA) tutorial https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/ 
#(Peter Langelder and Steve Horvath)
#################################################################################################
#Check working directory:
getwd()

#Load required packages:
if(!require("BiocManager", "plotrix")){
  install.packages("BiocManager")
}

BiocManager::install("WGCNA")
install.packages("tidyverse")
library(tidyverse)
library(WGCNA)
library(plotrix)
library(tidyverse)

# Important setting:
options(stringsAsFactors=FALSE) 

#Load RNA-seq data:
load(file = "Data Inputs/LGG_RNAseq.Rfile") #Loaded as LGG_RNA.seq
#################################################################################################

## 4. WGCNA Analysis for glioma (LGG) data set

#Remove uninformative data:
LGG = LGG_RNA.seq
row.names(LGG) = LGG_RNA.seq[,1]; colnames(LGG) = LGG_RNA.seq[1,]
LGG = LGG[-c(1,2),-1]

#Convert data frame to numeric matrix:
LGG1 = as.matrix(sapply(LGG, as.numeric))  
colnames(LGG1) = colnames(LGG)
row.names(LGG1) = row.names(LGG)
class(LGG1)
is.numeric(LGG1)

#Transpose matrix and filter for 10,000 most variant genes:
LGGdata = t(LGG1[order(apply(LGG1,1,mad), decreasing = T)[1:10000],])

#Check if LGGdata have many missing values:
gsg = goodSamplesGenes(LGGdata, verbose = 3)
gsg$allOK

#The command returns "TRUE", so all genes have passed the cuts.
#################################################################################################

# Re-cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers:
sampleTree = hclust(dist(LGGdata), method = "average")
# Plot the sample tree:
sizeGrWindow(12,9)
pdf(file = "Plots/Sample clustering to detect outliers_LGG.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()
#################################################################################################

#It appears there is 1 outlier. Outlier removal: 

# Plot a line to show the cut
sizeGrWindow(12,9)
pdf(file = "Plots/Sample clustering to detect outliers_LGG_outlier removal.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h = 2500000, col = "red")

dev.off()

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 2500000, minSize = 10)
table(clust)
# Clust 1 contains the samples we want to keep
keepSamples = (clust==1)
LGGdata1 = LGGdata[keepSamples, ]

# Now the LGGdata1 matrix is ready for WGCNA analysis.
#################################################################################################

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(LGGdata1, powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
pdf(file = "Plots/Soft-thresholding pwoer for LGG.pdf")
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()

# Based on the scale-free topology graph, the soft-thresholding power of 6 was chosen.
#################################################################################################

# Constructing the gene network and identifying modules:
net_LGG = blockwiseModules(LGGdata1, power = 6,
                              TOMType = "unsigned", minModuleSize = 30,
                              reassignThreshold = 0, mergeCutHeight = 0.25,
                              numericLabels = TRUE, pamRespectsDendro = FALSE,
                              saveTOMs = TRUE,
                              saveTOMFileBase = "LGGTOM",
                              verbose = 3)

# To see how many modules were identified and what the module sizes are, one can use table(net$colors).
table(net_LGG$colors)

#Now we can visualize the modules.

# Convert labels to colors for plotting
mergedColors_LGG = labels2colors(net_LGG$colors)
# Plot the dendrogram and the module colors underneath
pdf(file = "Plots/Cluster dendrogram for LGG.pdf")
plotDendroAndColors(net_LGG$dendrograms[[1]], mergedColors_LGG[net_LGG$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main = "Cluster dendrogram for LGG")
dev.off()
#################################################################################################

# We now save the module assignment and module eigengene information necessary for subsequent analysis:
moduleLabels_LGG = net_LGG$colors
moduleColors_LGG = labels2colors(net_LGG$colors)
MEs_LGG = net_LGG$MEs;
geneTree_LGG = net_LGG$dendrograms[[1]];
save(MEs_LGG, moduleLabels_LGG, moduleColors_LGG, geneTree_LGG,
     file = "Data Outputs/LGG_network_modulecolor_and_label.RData")

# Define numbers of genes and samples
nGenes = ncol(LGGdata1)
nSamples = nrow(LGGdata1)

# Recalculate MEs with color labels
MEs0_LGG = moduleEigengenes(LGGdata1, moduleColors_LGG)$eigengenes
MEs_LGG = orderMEs(MEs0_LGG)

#Calculate module membership to identify important genes. 
# names (colors) of the modules
modNames_LGG = substring(names(MEs_LGG), 3)

geneModuleMembership_LGG = as.data.frame(cor(LGGdata1, MEs_LGG, use = "p"));
MMPvalue_LGG = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_LGG), nSamples));

names(geneModuleMembership_LGG) = paste("MM", modNames_LGG, sep="");
names(MMPvalue_LGG) = paste("p.MM", modNames_LGG, sep="");

#We now create a data frame holding the following information for all genes: gene names, 
#module color,and module membership and p-values in all modules: 

# Create the starting data frame
geneInfoLGG = data.frame(Genename = colnames(LGGdata1),
                            moduleColor = moduleColors_LGG,
                            geneModuleMembership_LGG,
                            MMPvalue_LGG)

#Order the genes in the geneInfo variable by module color:
geneOrder_LGG = order(geneInfoLGG$moduleColor)
geneInfoLGG_1 = geneInfoLGG[geneOrder_LGG, ]

# Save the data frame into a text-format spreadsheet:
write.csv(geneInfoLGG_1, file = "Data Outputs/LGG_geneMM.csv")
#################################################################################################

# Now we calculate scaled connectivity of genes in the LGG network:

#Create a TOM matrix:
tom_LGG = TOMsimilarityFromExpr(LGGdata1)

#################################################################################################

# c. Extract LGG "yellow" module containing tyrobp:
     
     #Select modules
     modules = c("yellow")
     
     # Select module probes
     Genes = colnames(LGGdata1)
     inModule = is.finite(match(moduleColors_LGG, modules))
     modGenes =Genes[inModule]
     
     # Select the corresponding Topological Overlap
     modTOM = tom_LGG[inModule, inModule]
     dimnames(modTOM) = list(modGenes, modGenes)
     
     # Export the network into edge and node list files Cytoscape can read
     cyt = exportNetworkToCytoscape(modTOM,
                                    edgeFile = paste("LGG_edges_", paste(modules, collapse="_"), ".txt", sep=""),
                                    nodeFile = paste("LGG_nodes_", paste(modules, collapse="_"), ".txt", sep=""),
                                    weighted = TRUE,
                                    threshold = 0.02,
                                    nodeNames = modGenes,
                                    nodeAttr = moduleColors_LGG[inModule])
     