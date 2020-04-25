
# Load the package
library(WGCNA);
library(flashClust)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

#Read in the data
femData = read.table(file="for_WGCNA_gene.txt",sep="\t",header=T);
order_list = read.table("For_gene_order.txt",sep="\t",header = F)[,1]
femData<-femData[,order_list]

Sample_inf=femData[,1:2]
row.names(femData)=femData[,1]

allowWGCNAThreads()
#enableWGCNAThreads() # for RStudio

femData=data.frame(t(femData[,-1:-2]))

# Take a quick look at what is in the data set:
#dim(femData);
#names(femData);
#head(femData)

boxplot(femData[,c(3:25)],col="skyblue")

datExpr0 = as.data.frame(t(femData));
#head(datExpr0)


###check data###
gsg = goodSamplesGenes(datExpr0, verbose = 3);
#gsg$allOK

#If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the o
#ending genes and samples
#from the data:
if (!gsg$allOK)
{
# Optionally, print the gene and sample names that were removed:
if (sum(!gsg$goodGenes)>0)
printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
if (sum(!gsg$goodSamples)>0)
printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
# Remove the offending genes and samples from the data:
datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = flashClust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
#sizeGrWindow(12,9)
pdf(file = "sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
dev.off()
###filter data###
# Plot a line to show the cut
abline(h = 135, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 135, minSize = 3)
#table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
#head(datExpr)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#head(cor(datExpr))


###load content###
traitData = read.table(file="Clinic.txt",sep="\t",header = T);
# remove columns that hold information we do not need.
#allTraits = traitData[, -c(31, 16)];
#dim(allTraits)
allTraits = data.frame(traitData)

# Form a data frame analogous to expression data that will hold the clinical traits.
femaleSamples = rownames(datExpr);
traitRows = match(femaleSamples, allTraits$Sample);
datTraits = allTraits[traitRows, ];
rownames(datTraits) = allTraits[traitRows, 1];

collectGarbage();

# Re-cluster samples
sampleTree2 = flashClust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits[,3], signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
groupLabels = names(datTraits),
main = "Sample dendrogram and trait heatmap")

save(datExpr, datTraits, file = "FemaleLiver-01-dataInput.RData")


###Automatic network construction and module detection###
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.

# Load the data saved in the first part
lnames = load(file = "FemaleLiver-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
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


###one step network###
net = blockwiseModules(datExpr, power = sft$powerEstimate,
TOMType = "unsigned", minModuleSize = 30,
reassignThreshold = 0, mergeCutHeight = 0.35,
numericLabels = TRUE, pamRespectsDendro = FALSE,
saveTOMs = TRUE,
saveTOMFileBase = "femaleMouseTOM",
verbose = 3)

#table(net$colors)

# open a graphics window
sizeGrWindow(12, 9)

beta = 7
s = abs(bicor(datExpr0))
a = s^beta
#dissimilarity measure
w = 1-a

#create gene tree by average linkage hierarchical clustering 
geneTree = hclust(as.dist(w), method = 'average')
mergedColors = labels2colors(net$colors)

plotDendroAndColors(geneTree, mergedColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[3]], mergedColors[net$blockGenes[[3]]],
"Module colors",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
#geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
file = "FemaleLiver-02-networkConstruction-auto.RData")


###Relating modules to external information and identifying important genes###
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "FemaleLiver-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
#lnames
# Load network data saved in the second part.
lnames = load(file = "FemaleLiver-02-networkConstruction-auto.RData");
#lnames

###Quantifying module-trait associations###
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits$Point2, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

write.csv(MEs0,"MEs0.csv")

sizeGrWindow(10,10)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(0, 10, 0, 1));
# Display the correlation values within a heatmap plot

labeledHeatmap(Matrix = moduleTraitCor,
xLabels = 'Sample', # names(datTraits)
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.6,
zlim = c(-1,1),
main = paste("Module-trait relationships"))


# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$Point2);
names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));#xiugaiguo !!!!!!!

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");


All_mod<-data.frame(table(moduleColors),stringsAsFactors = F)[,1]
pdf("all_point.pdf",width =15,height = 10)

for (M in All_mod){
module = M
column = match(module, modNames);
moduleGenes = moduleColors==module;
#sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
abs(geneTraitSignificance[moduleGenes, 1]),
xlab = paste("Module Membership in", module, "module"),
ylab = "Gene significance for body weight",
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
}

dev.off()

###Summary output of network analysis results###
#names(datExpr)
#names(datExpr)[moduleColors=="red"]
library(readr)
annot <- read_delim("GPL4133_annot.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

probes = names(datExpr)
probes2annot = match(probes, annot$Gene_symbols)
# The following is the number or probes without annotation:
#sum(is.na(probes2annot))
# Should return 0.

# Create the starting data frame
geneInfo0 = data.frame(substanceBXH = probes,
geneSymbol = annot$Gene_symbols[probes2annot],
LocusLinkID = annot$Gene[probes2annot],
moduleColor = moduleColors,
geneTraitSignificance,
GSPvalue)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, weight, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
oldNames = names(geneInfo0)
geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
MMPvalue[, modOrder[mod]]);
names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight));
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "geneInfo.csv")


###Interfacing network analysis with other data such as functional annotation and gene ontology###
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "FemaleLiver-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
#lnames
# Load network data saved in the second part.
lnames = load(file = "FemaleLiver-02-networkConstruction-auto.RData");
#lnames

# Read in the probe annotation

# Match probes in the data set to the probe IDs in the annotation file
probes = names(datExpr)
probes2annot = match(probes, annot$Gene_symbols)
# Get the corresponding Locuis Link IDs
allLLIDs = annot$Gene[probes2annot];
# $ Choose interesting modules
intModules = All_mod

for (module in intModules)
{
# Select module probes
modGenes = (moduleColors==module)
# Get their entrez ID codes
modLLIDs = allLLIDs[modGenes];
# Write them into a file
fileName = paste("LocusLinkIDs-", module, ".txt", sep="");
write.table(as.data.frame(modLLIDs), file = fileName,
row.names = FALSE, col.names = FALSE)
}
# As background in the enrichment analysis, we will use all probes in the analysis.
fileName = paste("LocusLinkIDs-all.txt", sep="");
write.table(as.data.frame(allLLIDs), file = fileName,
row.names = FALSE, col.names = FALSE)


###Network visualization using WGCNA functions###
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "FemaleLiver-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
#lnames
# Load network data saved in the second part.
lnames = load(file = "FemaleLiver-02-networkConstruction-auto.RData");
#lnames
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

###Visualizing the gene network###
# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = sft$powerEstimate);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^sft$powerEstimate;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function

w=plotTOM


sizeGrWindow(9,9)
#geneTree = hclust(as.dist(w),method="average")
#TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
geneTree2 = net$dendrograms[[3]];
moduleColors2=mergedColors[net$blockGenes[[3]]]
plotTOM2=plotTOM[net$blockGenes[[3]],net$blockGenes[[3]]]
TOMplot(plotTOM2, geneTree2, moduleColors2, main = "Network heatmap plot, all genes")

nSelect = 1000
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = flashClust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^sft$powerEstimate;
diag(plotDiss) = NA;

TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")


###Visualizing the network of eigengenes###
# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes

# Isolate weight from the clinical traits
#weight = as.data.frame(datTraits$Point2);
#names(weight) = "weight"
# Add the weight to existing module eigengenes
#MET = orderMEs(cbind(MEs, weight))

MET <- orderMEs(MEs)

# Plot the relationships among the eigengenes and the trait
pdf("heat.map.pdf",width = 10,height = 6)
par(cex = 0.9)

plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90, excludeGrey = FALSE, greyLabel = "MEgrey")

# Plot the dendrogram

par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE, excludeGrey = FALSE, greyLabel = "MEgrey")
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90, excludeGrey = FALSE, greyLabel = "MEgrey")
dev.off()
###Exporting a gene network to external visualization software###
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "FemaleLiver-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
#lnames
# Load network data saved in the second part.
lnames = load(file = "FemaleLiver-02-networkConstruction-auto.RData");
#lnames



###Exporting to Cytoscape###
# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = sft$powerEstimate);
TOM=1-dissTOM


# Select modules

for (modules in intModules)
{
#modules = "yellow";
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modGenes = annot$Gene_symbols[match(modProbes, annot$Gene_symbols)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
weighted = TRUE,
threshold = 0.02,
nodeNames = modProbes,
altNodeNames = modGenes,
nodeAttr = moduleColors[inModule]);
}
###Exporting to VisANT###
# Recalculate topological overlap
#TOM = TOMsimilarityFromExpr(datExpr, power = 3);

# Read in the annotation file
#annot = read.csv(file = "GeneAnnotation_g123.csv");
# Select module

for (module in intModules)
{
#module = "yellow";
# Select module probes
probes = names(datExpr)
inModule = (moduleColors==module);
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
file = paste("VisANTInput-", module, ".txt", sep=""),
weighted = TRUE,
threshold = 0.15
#probeToGene = data.frame(annot$Gene_symbols, annot$Gene_symbols) 
)

nTop = 100;
IMConn = softConnectivity(datExpr[, modProbes]);
top = (rank(-IMConn) <= nTop)
vis = exportNetworkToVisANT(modTOM[top, top],
file = paste("VisANTInput-", module, "-top100.txt", sep=""),
weighted = TRUE,
threshold = 0.15
#probeToGene = data.frame(annot$Gene_symbols, annot$Gene_symbols) 
)

}

sessionInfo()
