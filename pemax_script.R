# Load packages
library(vegan)
library(DESeq2)
library(qiime2R)
library(phyloseq)
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(scales)
library(data.table)
library(fBasics)
library(forcats)
library(dplyr)
library(ggsci)

# Load data
physeq<-qza_to_phyloseq("rarefied-filtered-merged-table.qza",
                        "rooted-tree.qza",
                        "taxonomy.qza", 
                        "metadata.txt")

# Remove taxa with 0 abundance
otu.table = subset_taxa(physeq, rowSums(otu_table(physeq)) != 0)

# To normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
otu.relative.table = transformSampleCounts(otu.table, normalizeSample)

# Check
rownames(sample_data(otu.relative.table))
colnames(sample_data(otu.relative.table))

# Subset
cohortA.table = subset_samples(otu.relative.table, Cohort %in% c("A"))
cohortB.with.bkg.table = subset_samples(otu.relative.table, Cohort %in% c("B"))
cohortB.table = subset_samples(cohortB.with.bkg.table, Description %in% c("sputum"))
cohortA.pma.table = subset_samples(cohortA.table, Treatment %in% c("pma"))
cohortA.pemax.table = subset_samples(cohortA.table, Treatment %in% c("pemax"))
cohortA.dnase.table = subset_samples(cohortA.table, Treatment %in% c("dnase"))
cohortA.untreated.table = subset_samples(cohortA.table, Treatment %in% c("untreated"))

# Save
save.image(file="pemax.RData")

# Load
load(file="pemax.RData")

## ALPHA DIVERSITY

alpha <- estimate_richness(physeq)
write.csv(alpha, "alpha.diversity.csv")

## BETA DIVERSITY

## Figure 4B
# Create a distance matrix 
vegdist = distance(cohortA.untreated.table, "bray")

# Formulate principal component co-ordinates for PCoA plot, k is the choice of PCs
CmdScale <- cmdscale(vegdist, k = 10)

# Apply a function (variance) to the matrix
vars <- apply(CmdScale, 2, var)

# Create variable with the percent variance for each axis
percentVar <- round(100 * (vars/sum(vars)))

# Merge principal component data with metadata
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(cohortA.untreated.table), by = "row.names", all.x = TRUE)

# Rename variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

# Calculate the centroid value
centroids <- aggregate(cbind(PC1,PC2)~ Group,data= newResults, mean)

# Merge the centroid data into the PCOA data
newResults <- merge(newResults,centroids,by="Group",suffixes=c("",".centroid"))

# Plot
pdf("Figure4B.pdf", height = 10, width = 15)
ggplot(newResults, aes(PC1, PC2, color= Group)) +
  geom_point(size=3,alpha=0.7) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values=c("#4DDDDE", "#FF6F6F")) +  
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= Group), size=0) +
  theme_minimal() +
  theme(plot.title=element_text( hjust=1, vjust=0.5, face='bold', size=20)) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Group, alpha=0.2))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=Group, size=10),size=8, fontface = "bold") +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Statistics
bray = distance(cohortA.untreated.table, "bray")
adonis2(bray ~ sample_data(cohortA.untreated.table)$Group)

## Figure 5B
# Create a distance matrix 
vegdist = distance(cohortB.table, "bray")

# Formulate principal component co-ordinates for PCoA plot, k is the choice of PCs
CmdScale <- cmdscale(vegdist, k = 10)

# Apply a function (variance) to the matrix
vars <- apply(CmdScale, 2, var)

# Create variable with the percent variance for each axis
percentVar <- round(100 * (vars/sum(vars)))

# Merge principal component data with metadata
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(cohortB.table), by = "row.names", all.x = TRUE)

# Rename variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

# Calculate the centroid value
centroids <- aggregate(cbind(PC1,PC2)~ Timepoint_treatment,data= newResults, mean)

# Merge the centroid data into the PCOA data
newResults <- merge(newResults,centroids,by="Timepoint_treatment",suffixes=c("",".centroid"))

# Plot
pdf("Figure5B.pdf", height = 10, width = 15)
ggplot(newResults, aes(PC1, PC2, color= Timepoint_treatment)) +
  geom_point(size=3,alpha=0.7) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values=c("#59BE86", "#F253A9", "#FF964D", "#F5D27C")) +  
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= Timepoint_treatment), size=0) +
  #ggtitle("Stool") +
  theme_minimal() +
  theme(plot.title=element_text( hjust=1, vjust=0.5, face='bold', size=20)) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Timepoint_treatment, alpha=0.2))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=Timepoint_treatment, size=10),size=8, fontface = "bold") +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Statistics
bray = distance(cohortB.table, "bray")
vegan::adonis2(bray ~ sample_data(cohortB.table)$Timepoint_treatment)

# Subset
week0.vs.week0.pemax = subset_samples(cohortB.table, Timepoint_treatment %in% c("week0", "week0.PEMAX"))
week24.vs.week24.pemax = subset_samples(cohortB.table, Timepoint_treatment %in% c("week24", "week24.PEMAX"))
week0.vs.week24 = subset_samples(cohortB.table, Timepoint_treatment %in% c("week0", "week24"))
week0.pemax.vs.week24.pemax = subset_samples(cohortB.table, Timepoint_treatment %in% c("week0.PEMAX", "week24.PEMAX"))

# Check
rownames(sample_data(week0.vs.week0.pemax))
rownames(sample_data(week24.vs.week24.pemax))
rownames(sample_data(week0.vs.week24))
rownames(sample_data(week0.pemax.vs.week24.pemax))

bray = distance(week0.vs.week0.pemax, "bray")
adonis2(bray ~ sample_data(week0.vs.week0.pemax)$Timepoint_treatment)

bray = distance(week24.vs.week24.pemax, "bray")
adonis2(bray ~ sample_data(week24.vs.week24.pemax)$Timepoint_treatment)

bray = distance(week0.vs.week24, "bray")
adonis2(bray ~ sample_data(week0.vs.week24)$Timepoint_treatment)

bray = distance(week0.pemax.vs.week24.pemax, "bray")
adonis2(bray ~ sample_data(week0.pemax.vs.week24.pemax)$Timepoint_treatment)

## Supplementary Figure 9D
# Create a distance matrix 
vegdist = distance(cohortA.pma.table, "bray")

# Formulate principal component co-ordinates for PCoA plot, k is the choice of PCs
CmdScale <- cmdscale(vegdist, k = 10)

# Apply a function (variance) to the matrix
vars <- apply(CmdScale, 2, var)

# Create variable with the percent variance for each axis
percentVar <- round(100 * (vars/sum(vars)))

# Merge principal component data with metadata
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(cohortA.pma.table), by = "row.names", all.x = TRUE)

# Rename variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

# Calculate the centroid value
centroids <- aggregate(cbind(PC1,PC2)~ Group,data= newResults, mean)

# Merge the centroid data into the PCOA data
newResults <- merge(newResults,centroids,by="Group",suffixes=c("",".centroid"))

# Plot
pdf("SuppFigure9D.pdf", height = 10, width = 15)
ggplot(newResults, aes(PC1, PC2, color= Group)) +
  geom_point(size=3,alpha=0.7) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values=c("#4DDDDE", "#FF6F6F")) +  
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= Group), size=0) +
  theme_minimal() +
  theme(plot.title=element_text( hjust=1, vjust=0.5, face='bold', size=20)) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Group, alpha=0.2))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=Group, size=10),size=8, fontface = "bold") +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Statistics
bray = distance(cohortA.pma.table, "bray")
adonis2(bray ~ sample_data(cohortA.pma.table)$Group)

## Supplementary Figure 9E
# Create a distance matrix 
vegdist = distance(cohortA.pemax.table, "bray")

# Formulate principal component co-ordinates for PCoA plot, k is the choice of PCs
CmdScale <- cmdscale(vegdist, k = 10)

# Apply a function (variance) to the matrix
vars <- apply(CmdScale, 2, var)

# Create variable with the percent variance for each axis
percentVar <- round(100 * (vars/sum(vars)))

# Merge principal component data with metadata
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(cohortA.pemax.table), by = "row.names", all.x = TRUE)

# Rename variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

# Calculate the centroid value
centroids <- aggregate(cbind(PC1,PC2)~ Group,data= newResults, mean)

# Merge the centroid data into the PCOA data
newResults <- merge(newResults,centroids,by="Group",suffixes=c("",".centroid"))

# Plot
pdf("SuppFigure9E.pdf", height = 10, width = 15)
ggplot(newResults, aes(PC1, PC2, color= Group)) +
  geom_point(size=3,alpha=0.7) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values=c("#4DDDDE", "#FF6F6F")) +  
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= Group), size=0) +
  theme_minimal() +
  theme(plot.title=element_text( hjust=1, vjust=0.5, face='bold', size=20)) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Group, alpha=0.2))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=Group, size=10),size=8, fontface = "bold") +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Statistics
bray = distance(cohortA.pemax.table, "bray")
adonis2(bray ~ sample_data(cohortA.pemax.table)$Group)

## Supplementary Figure 9F
# Create a distance matrix 
vegdist = distance(cohortA.dnase.table, "bray")

# Formulate principal component co-ordinates for PCoA plot, k is the choice of PCs
CmdScale <- cmdscale(vegdist, k = 10)

# Apply a function (variance) to the matrix
vars <- apply(CmdScale, 2, var)

# Create variable with the percent variance for each axis
percentVar <- round(100 * (vars/sum(vars)))

# Merge principal component data with metadata
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(cohortA.dnase.table), by = "row.names", all.x = TRUE)

# Rename variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

# Calculate the centroid value
centroids <- aggregate(cbind(PC1,PC2)~ Group,data= newResults, mean)

# Merge the centroid data into the PCOA data
newResults <- merge(newResults,centroids,by="Group",suffixes=c("",".centroid"))

# Plot
pdf("SuppFigure9F.pdf", height = 10, width = 15)
ggplot(newResults, aes(PC1, PC2, color= Group)) +
  geom_point(size=3,alpha=0.7) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values=c("#4DDDDE", "#FF6F6F")) +  
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= Group), size=0) +
  theme_minimal() +
  theme(plot.title=element_text( hjust=1, vjust=0.5, face='bold', size=20)) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Group, alpha=0.2))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=Group, size=10),size=8, fontface = "bold") +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Statistics
bray = distance(cohortA.dnase.table, "bray")
adonis2(bray ~ sample_data(cohortA.dnase.table)$Group)

## DIFFERENTIAL ABUNDANCE ANALYSIS

# Load required libraries
library(vegan)
library(DESeq2)
library(qiime2R)
library(phyloseq)
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(scales)
library(data.table)
library(fBasics)
library(forcats)
library(dplyr)
library(ggsci)

# Create new phyloseq object using unrarefied table
physeq<-qza_to_phyloseq("filtered-merged-table.qza",
                        "rooted-tree.qza",
                        "taxonomy.qza", 
                        "metadata.txt")

# Define functions
normalizeSample <- function(x) {
  x / sum(x)
}

gm_mean <- function(x, na.rm = TRUE) {
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

# Remove taxa with 0 abundance
otu.table = subset_taxa(physeq, rowSums(otu_table(physeq)) != 0)

# Set Theme For Figures
theme <- theme(
  panel.background = element_blank(),
  panel.border = element_rect(fill = NA),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  axis.text.x = element_text(colour = "black"),
  axis.text.y = element_text(colour = "black"),
  axis.ticks = element_line(colour = "black"),
  plot.margin = unit(c(1, 1, 1, 1), "line"),
  legend.position = "none"
)

# Choose Alpha/FDR
alpha = 0.2

# Subset to genus level
genus.table = tax_glom(otu.table, taxrank = "Genus")

# Plot
plot_bar(genus.table, fill="Genus") + theme(legend.position = "none")

# Check number of taxa
ntaxa(genus.table)

# Prune data to keep ~70% of genera (>1% relative abundance in >1% of samples)
filtered.genus.table = genefilter_sample(genus.table, filterfun_sample(function(x) x > 1), A = 0.01 * nsamples(genus.table))
overall.pruned.genus.table = prune_taxa(filtered.genus.table, genus.table)

# Check number of taxa after pruning
ntaxa(overall.pruned.genus.table)

# Save
save.image(file="pemax.deseq.RData")

# Load
load(file="pemax.deseq.RData")

# Figure 3E: Overall - untreated vs. PMA
# Subset genus table
pruned.genus.table = subset_samples(overall.pruned.genus.table, Cohort %in% c("A"))
pruned.genus.table = subset_samples(pruned.genus.table, Treatment %in% c("untreated", "pma"))

# Check
rownames(sample_data(pruned.genus.table))

# Create relative abundance table
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# Convert phyloseq object to DESEq object
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ Treatment)

# Calculate geometric means prior to estimate size factor
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$Treatment <- droplevels(diagdds$Treatment)

# Choose reference variable
diagdds$Treatment <- relevel(diagdds$Treatment, ref ="untreated")

# Run the differential analysis
diagdds<- DESeq(diagdds)

# Output results to a table
res <- results(diagdds)

# Reorder results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant taxa
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Kingdom, res$Phylum, res$Class, res$Order, res$Family, res$Genus, res$Species, rownames(res))

# Replace spaces with .
res$row2 <- gsub('\\s+', '|', res$row2)

# Convert Results table into a data.frame
res <- as.data.frame(res)

# Set Names of Results Table
res <- as.data.frame(res)
res <- setNames(cbind(rownames(res), res, row.names = NULL), 
                c("Taxa", "baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "row2"))
res$names <- res$Taxa
res$Taxa <- res$row2

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
test.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Treatment_c %in% c("2"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Treatment_c %in% c("0"))

# Calculate mean across the rows of the feature table
test.pruned.genus.rel.table.df <- data.frame(otu_table(test.pruned.genus.rel.table))
test.pruned.genus.rel.table.df.meanRA <- rowMeans(test.pruned.genus.rel.table.df)
control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)

# Subset AND reorder just the ASVs we have
test.pruned.genus.rel.table.df.meanRA.save <- test.pruned.genus.rel.table.df.meanRA[otu.to.save]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res data frame
res$abundance.test <- test.pruned.genus.rel.table.df.meanRA.save
res$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need
res.1 <- res[,c("Taxa", "abundance.test", "abundance.control", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Output table
write.table(res.1,file="Figure3E.abundances.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Compute FDR in log control
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "firebrick1"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "mediumorchid3"

# Create a variable for the size of the dots in the volcano plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="Figure3E.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Taxa)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 200 * res$abundance.test, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 200 * res$abundance.control,2)),alpha=0.8) + 
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=4,force=25, segment.colour="grey",segment.alpha=0.2) + 
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + 
  xlab("Effect size: log2(fold-change)") + 
  ylab("-log10(adjusted p-value)") + 
  theme 
dev.off()

# Figure 3F: Overall - untreated vs. PEMAX
# Subset genus table
pruned.genus.table = subset_samples(overall.pruned.genus.table, Cohort %in% c("A"))
pruned.genus.table = subset_samples(pruned.genus.table, Treatment %in% c("untreated", "pemax"))

# Check
rownames(sample_data(pruned.genus.table))

# Create relative abundance table
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# Convert phyloseq object to DESEq object
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ Treatment)

# Calculate geometric means prior to estimate size factor
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$Treatment <- droplevels(diagdds$Treatment)

# Choose reference variable
diagdds$Treatment <- relevel(diagdds$Treatment, ref ="untreated")

# Run the differential analysis
diagdds<- DESeq(diagdds)

# Output results to a table
res <- results(diagdds)

# Reorder results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant taxa
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Kingdom, res$Phylum, res$Class, res$Order, res$Family, res$Genus, res$Species, rownames(res))

# Replace spaces with .
res$row2 <- gsub('\\s+', '|', res$row2)

# Convert Results table into a data.frame
res <- as.data.frame(res)

# Set Names of Results Table
res <- as.data.frame(res)
res <- setNames(cbind(rownames(res), res, row.names = NULL), 
                c("Taxa", "baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "row2"))
res$names <- res$Taxa
res$Taxa <- res$row2

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
test.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Treatment_c %in% c("3"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Treatment_c %in% c("0"))

# Calculate mean across the rows of the feature table
test.pruned.genus.rel.table.df <- data.frame(otu_table(test.pruned.genus.rel.table))
test.pruned.genus.rel.table.df.meanRA <- rowMeans(test.pruned.genus.rel.table.df)
control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)

# Subset AND reorder just the ASVs we have
test.pruned.genus.rel.table.df.meanRA.save <- test.pruned.genus.rel.table.df.meanRA[otu.to.save]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res data frame
res$abundance.test <- test.pruned.genus.rel.table.df.meanRA.save
res$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need
res.1 <- res[,c("Taxa", "abundance.test", "abundance.control", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Output table
write.table(res.1,file="Figure3F.abundances.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Compute FDR in log control
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "forestgreen"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "mediumorchid3"

# Create a variable for the size of the dots in the volcano plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="Figure3F.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Taxa)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 200 * res$abundance.test, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 200 * res$abundance.control,2)),alpha=0.8) + 
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=4,force=25, segment.colour="grey",segment.alpha=0.2) + 
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + 
  xlab("Effect size: log2(fold-change)") + 
  ylab("-log10(adjusted p-value)") + 
  theme 
dev.off()

# Overall - untreated vs. DNaseI
# Subset genus table
pruned.genus.table = subset_samples(overall.pruned.genus.table, Cohort %in% c("A"))
pruned.genus.table = subset_samples(pruned.genus.table, Treatment %in% c("untreated", "dnase"))

# Check
rownames(sample_data(pruned.genus.table))

# Create relative abundance table
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# Convert phyloseq object to DESEq object
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ Treatment)

# Calculate geometric means prior to estimate size factor
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$Treatment <- droplevels(diagdds$Treatment)

# Choose reference variable
diagdds$Treatment <- relevel(diagdds$Treatment, ref ="untreated")

# Run the differential analysis
diagdds<- DESeq(diagdds)

# Output results to a table
res <- results(diagdds)

# Reorder results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant taxa
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Kingdom, res$Phylum, res$Class, res$Order, res$Family, res$Genus, res$Species, rownames(res))

# Replace spaces with .
res$row2 <- gsub('\\s+', '|', res$row2)

# Convert Results table into a data.frame
res <- as.data.frame(res)

# Set Names of Results Table
res <- as.data.frame(res)
res <- setNames(cbind(rownames(res), res, row.names = NULL), 
                c("Taxa", "baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "row2"))
res$names <- res$Taxa
res$Taxa <- res$row2

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
test.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Treatment_c %in% c("1"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Treatment_c %in% c("0"))

# Calculate mean across the rows of the feature table
test.pruned.genus.rel.table.df <- data.frame(otu_table(test.pruned.genus.rel.table))
test.pruned.genus.rel.table.df.meanRA <- rowMeans(test.pruned.genus.rel.table.df)
control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)

# Subset AND reorder just the ASVs we have
test.pruned.genus.rel.table.df.meanRA.save <- test.pruned.genus.rel.table.df.meanRA[otu.to.save]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res data frame
res$abundance.test <- test.pruned.genus.rel.table.df.meanRA.save
res$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need
res.1 <- res[,c("Taxa", "abundance.test", "abundance.control", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Output table
write.table(res.1,file="dnase.vs.untreated.abundances.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Compute FDR in log control
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "orange"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "mediumorchid3"

# Create a variable for the size of the dots in the volcano plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="dnase.vs.untreated.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Taxa)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 200 * res$abundance.test, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 200 * res$abundance.control,2)),alpha=0.8) + 
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=4,force=25, segment.colour="grey",segment.alpha=0.2) + 
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + 
  xlab("Effect size: log2(fold-change)") + 
  ylab("-log10(adjusted p-value)") + 
  theme 
dev.off()

# Supplementary Figure 8A - TB-NEGATIVES: UNTREATED VS PMA
# Subset genus table
pruned.genus.table = subset_samples(overall.pruned.genus.table, Cohort %in% c("A"))
pruned.genus.table = subset_samples(pruned.genus.table, Treatment %in% c("untreated", "pma"))
pruned.genus.table = subset_samples(pruned.genus.table, Group %in% c("TB-negative"))

# Check
rownames(sample_data(pruned.genus.table))

# Create relative abundance table
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# Convert phyloseq object to DESEq object
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ Treatment)

# Calculate geometric means prior to estimate size factor
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$Treatment <- droplevels(diagdds$Treatment)

# Choose reference variable
diagdds$Treatment <- relevel(diagdds$Treatment, ref ="untreated")

# Run the differential analysis
diagdds<- DESeq(diagdds)

# Output results to a table
res <- results(diagdds)

# Reorder results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant taxa
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Kingdom, res$Phylum, res$Class, res$Order, res$Family, res$Genus, res$Species, rownames(res))

# Replace spaces with .
res$row2 <- gsub('\\s+', '|', res$row2)

# Convert Results table into a data.frame
res <- as.data.frame(res)

# Set Names of Results Table
res <- as.data.frame(res)
res <- setNames(cbind(rownames(res), res, row.names = NULL), 
                c("Taxa", "baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "row2"))
res$names <- res$Taxa
res$Taxa <- res$row2

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
test.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Treatment_c %in% c("2"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Treatment_c %in% c("0"))

# Calculate mean across the rows of the feature table
test.pruned.genus.rel.table.df <- data.frame(otu_table(test.pruned.genus.rel.table))
test.pruned.genus.rel.table.df.meanRA <- rowMeans(test.pruned.genus.rel.table.df)
control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)

# Subset AND reorder just the ASVs we have
test.pruned.genus.rel.table.df.meanRA.save <- test.pruned.genus.rel.table.df.meanRA[otu.to.save]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res data frame
res$abundance.test <- test.pruned.genus.rel.table.df.meanRA.save
res$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need
res.1 <- res[,c("Taxa", "abundance.test", "abundance.control", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Output table
write.table(res.1,file="SuppFigure8A.abundances.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Compute FDR in log control
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "firebrick1"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "mediumorchid3"

# Create a variable for the size of the dots in the volcano plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="SuppFigure8A.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Taxa)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 200 * res$abundance.test, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 200 * res$abundance.control,2)),alpha=0.8) + 
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=4,force=25, segment.colour="grey",segment.alpha=0.2) + 
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + 
  xlab("Effect size: log2(fold-change)") + 
  ylab("-log10(adjusted p-value)") + 
  theme 
dev.off()

# Supplementary Figure 8B - TB-POSITIVES: UNTREATED VS PMA
# Subset genus table
pruned.genus.table = subset_samples(overall.pruned.genus.table, Cohort %in% c("A"))
pruned.genus.table = subset_samples(pruned.genus.table, Treatment %in% c("untreated", "pma"))
pruned.genus.table = subset_samples(pruned.genus.table, Group %in% c("TB-positive"))

# Check
rownames(sample_data(pruned.genus.table))

# Create relative abundance table
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# Convert phyloseq object to DESEq object
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ Treatment)

# Calculate geometric means prior to estimate size factor
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$Treatment <- droplevels(diagdds$Treatment)

# Choose reference variable
diagdds$Treatment <- relevel(diagdds$Treatment, ref ="untreated")

# Run the differential analysis
diagdds<- DESeq(diagdds)

# Output results to a table
res <- results(diagdds)

# Reorder results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant taxa
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Kingdom, res$Phylum, res$Class, res$Order, res$Family, res$Genus, res$Species, rownames(res))

# Replace spaces with .
res$row2 <- gsub('\\s+', '|', res$row2)

# Convert Results table into a data.frame
res <- as.data.frame(res)

# Set Names of Results Table
res <- as.data.frame(res)
res <- setNames(cbind(rownames(res), res, row.names = NULL), 
                c("Taxa", "baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "row2"))
res$names <- res$Taxa
res$Taxa <- res$row2

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
test.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Treatment_c %in% c("2"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Treatment_c %in% c("0"))

# Calculate mean across the rows of the feature table
test.pruned.genus.rel.table.df <- data.frame(otu_table(test.pruned.genus.rel.table))
test.pruned.genus.rel.table.df.meanRA <- rowMeans(test.pruned.genus.rel.table.df)
control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)

# Subset AND reorder just the ASVs we have
test.pruned.genus.rel.table.df.meanRA.save <- test.pruned.genus.rel.table.df.meanRA[otu.to.save]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res data frame
res$abundance.test <- test.pruned.genus.rel.table.df.meanRA.save
res$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need
res.1 <- res[,c("Taxa", "abundance.test", "abundance.control", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Output table
write.table(res.1,file="SuppFigure8B.abundances.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Compute FDR in log control
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "firebrick1"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "mediumorchid3"

# Create a variable for the size of the dots in the volcano plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="SuppFigure8B.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Taxa)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 200 * res$abundance.test, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 200 * res$abundance.control,2)),alpha=0.8) + 
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=4,force=25, segment.colour="grey",segment.alpha=0.2) + 
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + 
  xlab("Effect size: log2(fold-change)") + 
  ylab("-log10(adjusted p-value)") + 
  theme 
dev.off()

# TB-NEGATIVES: UNTREATED VS PEMAX
# Subset genus table
pruned.genus.table = subset_samples(overall.pruned.genus.table, Cohort %in% c("A"))
pruned.genus.table = subset_samples(pruned.genus.table, Treatment %in% c("untreated", "pemax"))
pruned.genus.table = subset_samples(pruned.genus.table, Group %in% c("TB-negative"))

# Check
rownames(sample_data(pruned.genus.table))

# Create relative abundance table
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# Convert phyloseq object to DESEq object
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ Treatment)

# Calculate geometric means prior to estimate size factor
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$Treatment <- droplevels(diagdds$Treatment)

# Choose reference variable
diagdds$Treatment <- relevel(diagdds$Treatment, ref ="untreated")

# Run the differential analysis
diagdds<- DESeq(diagdds)

# Output results to a table
res <- results(diagdds)

# Reorder results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant taxa
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Kingdom, res$Phylum, res$Class, res$Order, res$Family, res$Genus, res$Species, rownames(res))

# Replace spaces with .
res$row2 <- gsub('\\s+', '|', res$row2)

# Convert Results table into a data.frame
res <- as.data.frame(res)

# Set Names of Results Table
res <- as.data.frame(res)
res <- setNames(cbind(rownames(res), res, row.names = NULL), 
                c("Taxa", "baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "row2"))
res$names <- res$Taxa
res$Taxa <- res$row2

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
test.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Treatment_c %in% c("3"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Treatment_c %in% c("0"))

# Calculate mean across the rows of the feature table
test.pruned.genus.rel.table.df <- data.frame(otu_table(test.pruned.genus.rel.table))
test.pruned.genus.rel.table.df.meanRA <- rowMeans(test.pruned.genus.rel.table.df)
control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)

# Subset AND reorder just the ASVs we have
test.pruned.genus.rel.table.df.meanRA.save <- test.pruned.genus.rel.table.df.meanRA[otu.to.save]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res data frame
res$abundance.test <- test.pruned.genus.rel.table.df.meanRA.save
res$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need
res.1 <- res[,c("Taxa", "abundance.test", "abundance.control", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Output table
write.table(res.1,file="TB-negatives.untreated.vs.pemax.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Compute FDR in log control
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "forestgreen"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "mediumorchid3"

# Create a variable for the size of the dots in the volcano plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="TB-negatives.untreated.vs.pemax.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Taxa)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 200 * res$abundance.test, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 200 * res$abundance.control,2)),alpha=0.7) + 
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=4,force=25, segment.colour="grey",segment.alpha=0.2) + 
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + 
  xlab("Effect size: log2(fold-change)") + 
  ylab("-log10(adjusted p-value)") + 
  theme 
dev.off()

# Supplementary Figure 8C - TB-POSITIVES: UNTREATED VS PEMAX
# Subset genus table
pruned.genus.table = subset_samples(overall.pruned.genus.table, Cohort %in% c("A"))
pruned.genus.table = subset_samples(pruned.genus.table, Treatment %in% c("untreated", "pemax"))
pruned.genus.table = subset_samples(pruned.genus.table, Group %in% c("TB-positive"))

# Check
rownames(sample_data(pruned.genus.table))

# Create relative abundance table
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# Convert phyloseq object to DESEq object
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ Treatment)

# Calculate geometric means prior to estimate size factor
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$Treatment <- droplevels(diagdds$Treatment)

# Choose reference variable
diagdds$Treatment <- relevel(diagdds$Treatment, ref ="untreated")

# Run the differential analysis
diagdds<- DESeq(diagdds)

# Output results to a table
res <- results(diagdds)

# Reorder results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant taxa
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Kingdom, res$Phylum, res$Class, res$Order, res$Family, res$Genus, res$Species, rownames(res))

# Replace spaces with .
res$row2 <- gsub('\\s+', '|', res$row2)

# Convert Results table into a data.frame
res <- as.data.frame(res)

# Set Names of Results Table
res <- as.data.frame(res)
res <- setNames(cbind(rownames(res), res, row.names = NULL), 
                c("Taxa", "baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "row2"))
res$names <- res$Taxa
res$Taxa <- res$row2

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
test.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Treatment_c %in% c("3"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Treatment_c %in% c("0"))

# Calculate mean across the rows of the feature table
test.pruned.genus.rel.table.df <- data.frame(otu_table(test.pruned.genus.rel.table))
test.pruned.genus.rel.table.df.meanRA <- rowMeans(test.pruned.genus.rel.table.df)
control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)

# Subset AND reorder just the ASVs we have
test.pruned.genus.rel.table.df.meanRA.save <- test.pruned.genus.rel.table.df.meanRA[otu.to.save]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res data frame
res$abundance.test <- test.pruned.genus.rel.table.df.meanRA.save
res$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need
res.1 <- res[,c("Taxa", "abundance.test", "abundance.control", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Output table
write.table(res.1,file="SuppFigure8C.abundances.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Compute FDR in log control
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "forestgreen"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "mediumorchid3"

# Create a variable for the size of the dots in the volcano plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="SuppFigure8C.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Taxa)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 200 * res$abundance.test, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 200 * res$abundance.control,2)),alpha=0.8) + 
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=4,force=25, segment.colour="grey",segment.alpha=0.2) + 
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + 
  xlab("Effect size: log2(fold-change)") + 
  ylab("-log10(adjusted p-value)") + 
  theme 
dev.off()

# TB-NEGATIVES: UNTREATED VS DNASE I
# Subset genus table
pruned.genus.table = subset_samples(overall.pruned.genus.table, Cohort %in% c("A"))
pruned.genus.table = subset_samples(pruned.genus.table, Treatment %in% c("untreated", "dnase"))
pruned.genus.table = subset_samples(pruned.genus.table, Group %in% c("TB-negative"))

# Check
rownames(sample_data(pruned.genus.table))

# Create relative abundance table
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# Convert phyloseq object to DESEq object
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ Treatment)

# Calculate geometric means prior to estimate size factor
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$Treatment <- droplevels(diagdds$Treatment)

# Choose reference variable
diagdds$Treatment <- relevel(diagdds$Treatment, ref ="untreated")

# Run the differential analysis
diagdds<- DESeq(diagdds)

# Output results to a table
res <- results(diagdds)

# Reorder results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant taxa
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Kingdom, res$Phylum, res$Class, res$Order, res$Family, res$Genus, res$Species, rownames(res))

# Replace spaces with .
res$row2 <- gsub('\\s+', '|', res$row2)

# Convert Results table into a data.frame
res <- as.data.frame(res)

# Set Names of Results Table
res <- as.data.frame(res)
res <- setNames(cbind(rownames(res), res, row.names = NULL), 
                c("Taxa", "baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "row2"))
res$names <- res$Taxa
res$Taxa <- res$row2

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
test.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Treatment_c %in% c("1"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Treatment_c %in% c("0"))

# Calculate mean across the rows of the feature table
test.pruned.genus.rel.table.df <- data.frame(otu_table(test.pruned.genus.rel.table))
test.pruned.genus.rel.table.df.meanRA <- rowMeans(test.pruned.genus.rel.table.df)
control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)

# Subset AND reorder just the ASVs we have
test.pruned.genus.rel.table.df.meanRA.save <- test.pruned.genus.rel.table.df.meanRA[otu.to.save]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res data frame
res$abundance.test <- test.pruned.genus.rel.table.df.meanRA.save
res$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need
res.1 <- res[,c("Taxa", "abundance.test", "abundance.control", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Output table
write.table(res.1,file="tb-negatives.untreated.vs.dnase.abundances.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Compute FDR in log control
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "orange"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "mediumorchid3"

# Create a variable for the size of the dots in the volcano plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="tb-negatives.untreated.vs.dnase.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Taxa)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 200 * res$abundance.test, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 200 * res$abundance.control,2)),alpha=0.8) + 
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=4,force=25, segment.colour="grey",segment.alpha=0.2) + 
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + 
  xlab("Effect size: log2(fold-change)") + 
  ylab("-log10(adjusted p-value)") + 
  theme 
dev.off()

# TB-POSITIVES: UNTREATED VS DNASE I
# Subset genus table
pruned.genus.table = subset_samples(overall.pruned.genus.table, Cohort %in% c("A"))
pruned.genus.table = subset_samples(pruned.genus.table, Treatment %in% c("untreated", "dnase"))
pruned.genus.table = subset_samples(pruned.genus.table, Group %in% c("TB-positive"))

# Check
rownames(sample_data(pruned.genus.table))

# Create relative abundance table
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# Convert phyloseq object to DESEq object
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ Treatment)

# Calculate geometric means prior to estimate size factor
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$Treatment <- droplevels(diagdds$Treatment)

# Choose reference variable
diagdds$Treatment <- relevel(diagdds$Treatment, ref ="untreated")

# Run the differential analysis
diagdds<- DESeq(diagdds)

# Output results to a table
res <- results(diagdds)

# Reorder results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant taxa
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Kingdom, res$Phylum, res$Class, res$Order, res$Family, res$Genus, res$Species, rownames(res))

# Replace spaces with .
res$row2 <- gsub('\\s+', '|', res$row2)

# Convert Results table into a data.frame
res <- as.data.frame(res)

# Set Names of Results Table
res <- as.data.frame(res)
res <- setNames(cbind(rownames(res), res, row.names = NULL), 
                c("Taxa", "baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "row2"))
res$names <- res$Taxa
res$Taxa <- res$row2

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
test.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Treatment_c %in% c("1"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Treatment_c %in% c("0"))

# Calculate mean across the rows of the feature table
test.pruned.genus.rel.table.df <- data.frame(otu_table(test.pruned.genus.rel.table))
test.pruned.genus.rel.table.df.meanRA <- rowMeans(test.pruned.genus.rel.table.df)
control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)

# Subset AND reorder just the ASVs we have
test.pruned.genus.rel.table.df.meanRA.save <- test.pruned.genus.rel.table.df.meanRA[otu.to.save]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res data frame
res$abundance.test <- test.pruned.genus.rel.table.df.meanRA.save
res$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need
res.1 <- res[,c("Taxa", "abundance.test", "abundance.control", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Output table
write.table(res.1,file="tb-positives.untreated.vs.dnase.abundances.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Compute FDR in log control
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "orange"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "mediumorchid3"

# Create a variable for the size of the dots in the volcano plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="tb-positives.untreated.vs.dnase.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Taxa)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 200 * res$abundance.test, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 200 * res$abundance.control,2)),alpha=0.8) + 
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=4,force=25, segment.colour="grey",segment.alpha=0.2) + 
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + 
  xlab("Effect size: log2(fold-change)") + 
  ylab("-log10(adjusted p-value)") + 
  theme 
dev.off()

# Figure 4C - UNTREATED: TB-NEGATIVES VS. TB-POSITIVES
# Subset genus table
pruned.genus.table = subset_samples(overall.pruned.genus.table, Cohort %in% c("A"))
pruned.genus.table = subset_samples(pruned.genus.table, Treatment %in% c("untreated"))

# Check
rownames(sample_data(pruned.genus.table))

# Create relative abundance table
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# Convert phyloseq object to DESEq object
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ Group)

# Calculate geometric means prior to estimate size factor
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$Group <- droplevels(diagdds$Group)

# Choose reference variable
diagdds$Group <- relevel(diagdds$Group, ref ="TB-negative")

# Run the differential analysis
diagdds<- DESeq(diagdds)

# Output results to a table
res <- results(diagdds)

# Reorder results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant taxa
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Kingdom, res$Phylum, res$Class, res$Order, res$Family, res$Genus, res$Species, rownames(res))

# Replace spaces with .
res$row2 <- gsub('\\s+', '|', res$row2)

# Convert Results table into a data.frame
res <- as.data.frame(res)

# Set Names of Results Table
res <- as.data.frame(res)
res <- setNames(cbind(rownames(res), res, row.names = NULL), 
                c("Taxa", "baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "row2"))
res$names <- res$Taxa
res$Taxa <- res$row2

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
test.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Group_c %in% c("1"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Group_c %in% c("0"))

# Calculate mean across the rows of the feature table
test.pruned.genus.rel.table.df <- data.frame(otu_table(test.pruned.genus.rel.table))
test.pruned.genus.rel.table.df.meanRA <- rowMeans(test.pruned.genus.rel.table.df)
control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)

# Subset AND reorder just the ASVs we have
test.pruned.genus.rel.table.df.meanRA.save <- test.pruned.genus.rel.table.df.meanRA[otu.to.save]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res data frame
res$abundance.test <- test.pruned.genus.rel.table.df.meanRA.save
res$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need
res.1 <- res[,c("Taxa", "abundance.test", "abundance.control", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Output table
write.table(res.1,file="Figure4C.abundances.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Compute FDR in log control
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "#FF6F6F"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "#4DDDDE"

# Create a variable for the size of the dots in the volcano plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="Figure4C.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Taxa)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 200 * res$abundance.test, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 200 * res$abundance.control,2)),alpha=0.8) + 
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=4,force=25, segment.colour="grey",segment.alpha=0.2) + 
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + 
  xlab("Effect size: log2(fold-change)") + 
  ylab("-log10(adjusted p-value)") + 
  theme 
dev.off()

# Figure 4D - PMA: TB-NEGATIVES VS. TB-POSITIVES
# Subset genus table
pruned.genus.table = subset_samples(overall.pruned.genus.table, Cohort %in% c("A"))
pruned.genus.table = subset_samples(pruned.genus.table, Treatment %in% c("pma"))

# Check
rownames(sample_data(pruned.genus.table))

# Create relative abundance table
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# Convert phyloseq object to DESEq object
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ Group)

# Calculate geometric means prior to estimate size factor
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$Group <- droplevels(diagdds$Group)

# Choose reference variable
diagdds$Group <- relevel(diagdds$Group, ref ="TB-negative")

# Run the differential analysis
diagdds<- DESeq(diagdds)

# Output results to a table
res <- results(diagdds)

# Reorder results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant taxa
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Kingdom, res$Phylum, res$Class, res$Order, res$Family, res$Genus, res$Species, rownames(res))

# Replace spaces with .
res$row2 <- gsub('\\s+', '|', res$row2)

# Convert Results table into a data.frame
res <- as.data.frame(res)

# Set Names of Results Table
res <- as.data.frame(res)
res <- setNames(cbind(rownames(res), res, row.names = NULL), 
                c("Taxa", "baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "row2"))
res$names <- res$Taxa
res$Taxa <- res$row2

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
test.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Group_c %in% c("1"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Group_c %in% c("0"))

# Calculate mean across the rows of the feature table
test.pruned.genus.rel.table.df <- data.frame(otu_table(test.pruned.genus.rel.table))
test.pruned.genus.rel.table.df.meanRA <- rowMeans(test.pruned.genus.rel.table.df)
control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)

# Subset AND reorder just the ASVs we have
test.pruned.genus.rel.table.df.meanRA.save <- test.pruned.genus.rel.table.df.meanRA[otu.to.save]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res data frame
res$abundance.test <- test.pruned.genus.rel.table.df.meanRA.save
res$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need
res.1 <- res[,c("Taxa", "abundance.test", "abundance.control", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Output table
write.table(res.1,file="Figure4D.abundances.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Compute FDR in log control
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "#FF6F6F"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "#4DDDDE"

# Create a variable for the size of the dots in the volcano plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="Figure4D.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Taxa)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 200 * res$abundance.test, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 200 * res$abundance.control,2)),alpha=0.8) + 
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=4,force=25, segment.colour="grey",segment.alpha=0.2) + 
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + 
  xlab("Effect size: log2(fold-change)") + 
  ylab("-log10(adjusted p-value)") + 
  theme 
dev.off()

# Figure 4E - PEMAX: TB-NEGATIVES VS. TB-POSITIVES
# Subset genus table
pruned.genus.table = subset_samples(overall.pruned.genus.table, Cohort %in% c("A"))
pruned.genus.table = subset_samples(pruned.genus.table, Treatment %in% c("pemax"))

# Check
rownames(sample_data(pruned.genus.table))

# Create relative abundance table
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# Convert phyloseq object to DESEq object
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ Group)

# Calculate geometric means prior to estimate size factor
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$Group <- droplevels(diagdds$Group)

# Choose reference variable
diagdds$Group <- relevel(diagdds$Group, ref ="TB-negative")

# Run the differential analysis
diagdds<- DESeq(diagdds)

# Output results to a table
res <- results(diagdds)

# Reorder results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant taxa
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Kingdom, res$Phylum, res$Class, res$Order, res$Family, res$Genus, res$Species, rownames(res))

# Replace spaces with .
res$row2 <- gsub('\\s+', '|', res$row2)

# Convert Results table into a data.frame
res <- as.data.frame(res)

# Set Names of Results Table
res <- as.data.frame(res)
res <- setNames(cbind(rownames(res), res, row.names = NULL), 
                c("Taxa", "baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "row2"))
res$names <- res$Taxa
res$Taxa <- res$row2

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
test.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Group_c %in% c("1"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Group_c %in% c("0"))

# Calculate mean across the rows of the feature table
test.pruned.genus.rel.table.df <- data.frame(otu_table(test.pruned.genus.rel.table))
test.pruned.genus.rel.table.df.meanRA <- rowMeans(test.pruned.genus.rel.table.df)
control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)

# Subset AND reorder just the ASVs we have
test.pruned.genus.rel.table.df.meanRA.save <- test.pruned.genus.rel.table.df.meanRA[otu.to.save]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res data frame
res$abundance.test <- test.pruned.genus.rel.table.df.meanRA.save
res$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need
res.1 <- res[,c("Taxa", "abundance.test", "abundance.control", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Output table
write.table(res.1,file="Figure4E.abundances.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Compute FDR in log control
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "#FF6F6F"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "#4DDDDE"

# Create a variable for the size of the dots in the volcano plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="Figure4E.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Taxa)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 200 * res$abundance.test, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 200 * res$abundance.control,2)),alpha=0.8) + 
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=4,force=25, segment.colour="grey",segment.alpha=0.2) + 
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + 
  xlab("Effect size: log2(fold-change)") + 
  ylab("-log10(adjusted p-value)") + 
  theme 
dev.off()

# DNASE I: TB-NEGATIVES VS. TB-POSITIVES
# Subset genus table
pruned.genus.table = subset_samples(overall.pruned.genus.table, Cohort %in% c("A"))
pruned.genus.table = subset_samples(pruned.genus.table, Treatment %in% c("dnase"))

# Check
rownames(sample_data(pruned.genus.table))

# Create relative abundance table
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# Convert phyloseq object to DESEq object
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ Group)

# Calculate geometric means prior to estimate size factor
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$Group <- droplevels(diagdds$Group)

# Choose reference variable
diagdds$Group <- relevel(diagdds$Group, ref ="TB-negative")

# Run the differential analysis
diagdds<- DESeq(diagdds)

# Output results to a table
res <- results(diagdds)

# Reorder results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant taxa
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Kingdom, res$Phylum, res$Class, res$Order, res$Family, res$Genus, res$Species, rownames(res))

# Replace spaces with .
res$row2 <- gsub('\\s+', '|', res$row2)

# Convert Results table into a data.frame
res <- as.data.frame(res)

# Set Names of Results Table
res <- as.data.frame(res)
res <- setNames(cbind(rownames(res), res, row.names = NULL), 
                c("Taxa", "baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "row2"))
res$names <- res$Taxa
res$Taxa <- res$row2

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
test.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Group_c %in% c("1"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Group_c %in% c("0"))

# Calculate mean across the rows of the feature table
test.pruned.genus.rel.table.df <- data.frame(otu_table(test.pruned.genus.rel.table))
test.pruned.genus.rel.table.df.meanRA <- rowMeans(test.pruned.genus.rel.table.df)
control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)

# Subset AND reorder just the ASVs we have
test.pruned.genus.rel.table.df.meanRA.save <- test.pruned.genus.rel.table.df.meanRA[otu.to.save]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res data frame
res$abundance.test <- test.pruned.genus.rel.table.df.meanRA.save
res$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need
res.1 <- res[,c("Taxa", "abundance.test", "abundance.control", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Output table
write.table(res.1,file="dnase.tb-negatives.vs.tb-positives.abundances.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Compute FDR in log control
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "#FF6F6F"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "#4DDDDE"

# Create a variable for the size of the dots in the volcano plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="dnase.tb-negatives.vs.tb-positives.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Taxa)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 200 * res$abundance.test, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 200 * res$abundance.control,2)),alpha=0.8) + 
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=4,force=25, segment.colour="grey",segment.alpha=0.2) + 
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + 
  xlab("Effect size: log2(fold-change)") + 
  ylab("-log10(adjusted p-value)") + 
  theme 
dev.off()

# WEEK 0 - UNTREATED VS PEMAX
# Subset genus table
pruned.genus.table = subset_samples(overall.pruned.genus.table, Cohort %in% c("B"))
pruned.genus.table = subset_samples(pruned.genus.table, Timepoint %in% c("week0"))

# Check
rownames(sample_data(pruned.genus.table))

# Create relative abundance table
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# Convert phyloseq object to DESEq object
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ Timepoint_treatment)

# Calculate geometric means prior to estimate size factor
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$Timepoint_treatment <- droplevels(diagdds$Timepoint_treatment)

# Choose reference variable
diagdds$Timepoint_treatment <- relevel(diagdds$Timepoint_treatment, ref ="week0")

# Run the differential analysis
diagdds<- DESeq(diagdds)

# Output results to a table
res <- results(diagdds)

# Reorder results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant taxa
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Kingdom, res$Phylum, res$Class, res$Order, res$Family, res$Genus, res$Species, rownames(res))

# Replace spaces with .
res$row2 <- gsub('\\s+', '|', res$row2)

# Convert Results table into a data.frame
res <- as.data.frame(res)

# Set Names of Results Table
res <- as.data.frame(res)
res <- setNames(cbind(rownames(res), res, row.names = NULL), 
                c("Taxa", "baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "row2"))
res$names <- res$Taxa
res$Taxa <- res$row2

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
test.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Timepoint_treatment_c %in% c("1"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Timepoint_treatment_c %in% c("0"))

# Calculate mean across the rows of the feature table
test.pruned.genus.rel.table.df <- data.frame(otu_table(test.pruned.genus.rel.table))
test.pruned.genus.rel.table.df.meanRA <- rowMeans(test.pruned.genus.rel.table.df)
control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)

# Subset AND reorder just the ASVs we have
test.pruned.genus.rel.table.df.meanRA.save <- test.pruned.genus.rel.table.df.meanRA[otu.to.save]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res data frame
res$abundance.test <- test.pruned.genus.rel.table.df.meanRA.save
res$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need
res.1 <- res[,c("Taxa", "abundance.test", "abundance.control", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Output table
write.table(res.1,file="week0.untreated.vs.pemax.abundances.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Compute FDR in log control
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "#F253A9"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "#59BE86"

# Create a variable for the size of the dots in the volcano plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="week0.untreated.vs.pemax.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Taxa)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 200 * res$abundance.test, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 200 * res$abundance.control,2)),alpha=0.8) + 
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=4,force=25, segment.colour="grey",segment.alpha=0.2) + 
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + 
  xlab("Effect size: log2(fold-change)") + 
  ylab("-log10(adjusted p-value)") + 
  theme 
dev.off()

# FIGURE 5C: WEEK 24 - UNTREATED VS PEMAX
# Subset genus table
pruned.genus.table = subset_samples(overall.pruned.genus.table, Cohort %in% c("B"))
pruned.genus.table = subset_samples(pruned.genus.table, Timepoint %in% c("week24"))

# Check
rownames(sample_data(pruned.genus.table))

# Create relative abundance table
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# Convert phyloseq object to DESEq object
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ Timepoint_treatment)

# Calculate geometric means prior to estimate size factor
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$Timepoint_treatment <- droplevels(diagdds$Timepoint_treatment)

# Choose reference variable
diagdds$Timepoint_treatment <- relevel(diagdds$Timepoint_treatment, ref ="week24")

# Run the differential analysis
diagdds<- DESeq(diagdds)

# Output results to a table
res <- results(diagdds)

# Reorder results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant taxa
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Kingdom, res$Phylum, res$Class, res$Order, res$Family, res$Genus, res$Species, rownames(res))

# Replace spaces with .
res$row2 <- gsub('\\s+', '|', res$row2)

# Convert Results table into a data.frame
res <- as.data.frame(res)

# Set Names of Results Table
res <- as.data.frame(res)
res <- setNames(cbind(rownames(res), res, row.names = NULL), 
                c("Taxa", "baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "row2"))
res$names <- res$Taxa
res$Taxa <- res$row2

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
test.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Timepoint_treatment_c %in% c("3"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Timepoint_treatment_c %in% c("2"))

# Calculate mean across the rows of the feature table
test.pruned.genus.rel.table.df <- data.frame(otu_table(test.pruned.genus.rel.table))
test.pruned.genus.rel.table.df.meanRA <- rowMeans(test.pruned.genus.rel.table.df)
control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)

# Subset AND reorder just the ASVs we have
test.pruned.genus.rel.table.df.meanRA.save <- test.pruned.genus.rel.table.df.meanRA[otu.to.save]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res data frame
res$abundance.test <- test.pruned.genus.rel.table.df.meanRA.save
res$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need
res.1 <- res[,c("Taxa", "abundance.test", "abundance.control", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Output table
write.table(res.1,file="Figure5C.abundances.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Compute FDR in log control
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "#F5D27C"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "#FF964D"

# Create a variable for the size of the dots in the volcano plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="Figure5C.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Taxa)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 200 * res$abundance.test, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 200 * res$abundance.control,2)),alpha=0.8) + 
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=4,force=25, segment.colour="grey",segment.alpha=0.2) + 
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + 
  xlab("Effect size: log2(fold-change)") + 
  ylab("-log10(adjusted p-value)") + 
  theme 
dev.off()

# FIGURE 5D: UNTREATED - WEEK 0 VS. WEEK 24
# Subset genus table
pruned.genus.table = subset_samples(overall.pruned.genus.table, Cohort %in% c("B"))
pruned.genus.table = subset_samples(pruned.genus.table, Timepoint_treatment %in% c("week0", "week24"))

# Check
rownames(sample_data(pruned.genus.table))

# Create relative abundance table
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# Convert phyloseq object to DESEq object
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ Timepoint)

# Calculate geometric means prior to estimate size factor
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$Timepoint <- droplevels(diagdds$Timepoint)

# Choose reference variable
diagdds$Timepoint <- relevel(diagdds$Timepoint, ref ="week0")

# Run the differential analysis
diagdds<- DESeq(diagdds)

# Output results to a table
res <- results(diagdds)

# Reorder results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant taxa
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Kingdom, res$Phylum, res$Class, res$Order, res$Family, res$Genus, res$Species, rownames(res))

# Replace spaces with .
res$row2 <- gsub('\\s+', '|', res$row2)

# Convert Results table into a data.frame
res <- as.data.frame(res)

# Set Names of Results Table
res <- as.data.frame(res)
res <- setNames(cbind(rownames(res), res, row.names = NULL), 
                c("Taxa", "baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "row2"))
res$names <- res$Taxa
res$Taxa <- res$row2

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
test.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Timepoint_c %in% c("1"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Timepoint_c %in% c("0"))

# Calculate mean across the rows of the feature table
test.pruned.genus.rel.table.df <- data.frame(otu_table(test.pruned.genus.rel.table))
test.pruned.genus.rel.table.df.meanRA <- rowMeans(test.pruned.genus.rel.table.df)
control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)

# Subset AND reorder just the ASVs we have
test.pruned.genus.rel.table.df.meanRA.save <- test.pruned.genus.rel.table.df.meanRA[otu.to.save]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res data frame
res$abundance.test <- test.pruned.genus.rel.table.df.meanRA.save
res$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need
res.1 <- res[,c("Taxa", "abundance.test", "abundance.control", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Output table
write.table(res.1,file="Figure5D.abundances.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Compute FDR in log control
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "#FF964D"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "#59BE86"

# Create a variable for the size of the dots in the volcano plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="Figure5D.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Taxa)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 200 * res$abundance.test, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 200 * res$abundance.control,2)),alpha=0.8) + 
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=4,force=25, segment.colour="grey",segment.alpha=0.2) + 
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + 
  xlab("Effect size: log2(fold-change)") + 
  ylab("-log10(adjusted p-value)") + 
  theme 
dev.off()

# FIGURE 5E: PEMAX - WEEK 0 VS. WEEK 24
# Subset genus table
pruned.genus.table = subset_samples(overall.pruned.genus.table, Cohort %in% c("B"))
pruned.genus.table = subset_samples(pruned.genus.table, Timepoint_treatment %in% c("week0.PEMAX", "week24.PEMAX"))

# Check
rownames(sample_data(pruned.genus.table))

# Create relative abundance table
pruned.genus.rel.table = transformSampleCounts(pruned.genus.table, normalizeSample)

# Convert phyloseq object to DESEq object
diagdds <- phyloseq_to_deseq2(pruned.genus.table, ~ Timepoint)

# Calculate geometric means prior to estimate size factor
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)

# Drop rows with no data in your comparison variable
diagdds$Timepoint <- droplevels(diagdds$Timepoint)

# Choose reference variable
diagdds$Timepoint <- relevel(diagdds$Timepoint, ref ="week0")

# Run the differential analysis
diagdds<- DESeq(diagdds)

# Output results to a table
res <- results(diagdds)

# Reorder results based on FDR
res = res[order(res$padj, na.last = NA), ]

# Create list of top 50 significant taxa
select_genes = rownames(res[res$padj < alpha & !is.na(res$padj), ])[1:50]

# Get taxa names from phyloseq object
res = cbind(as(res, "data.frame"), as(tax_table(pruned.genus.table)[rownames(res), ], "matrix"))

# Replace OTU with taxa
res$row2 <- paste(res$Kingdom, res$Phylum, res$Class, res$Order, res$Family, res$Genus, res$Species, rownames(res))

# Replace spaces with .
res$row2 <- gsub('\\s+', '|', res$row2)

# Convert Results table into a data.frame
res <- as.data.frame(res)

# Set Names of Results Table
res <- as.data.frame(res)
res <- setNames(cbind(rownames(res), res, row.names = NULL), 
                c("Taxa", "baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "row2"))
res$names <- res$Taxa
res$Taxa <- res$row2

# Decide what otu to save 
otu.to.save <-as.character(res$names)

# Subset relative table
test.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Timepoint_c %in% c("1"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table, Timepoint_c %in% c("0"))

# Calculate mean across the rows of the feature table
test.pruned.genus.rel.table.df <- data.frame(otu_table(test.pruned.genus.rel.table))
test.pruned.genus.rel.table.df.meanRA <- rowMeans(test.pruned.genus.rel.table.df)
control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)

# Subset AND reorder just the ASVs we have
test.pruned.genus.rel.table.df.meanRA.save <- test.pruned.genus.rel.table.df.meanRA[otu.to.save]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save]

# Add the abundance data for the res data frame
res$abundance.test <- test.pruned.genus.rel.table.df.meanRA.save
res$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

# Keep only the variables you need
res.1 <- res[,c("Taxa", "abundance.test", "abundance.control", "baseMean", "logFC", "pvalue", "adj.P.Val")]

# Output table
write.table(res.1,file="Figure5E.abundances.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

# Compute FDR in log control
res$sig <- -log10(res$adj.P.Val)

# See how many are now infinite
sum(is.infinite(res$sig))

# Set the colors for your volcano plot
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "#F5D27C"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "#F253A9"

# Create a variable for the size of the dots in the volcano plot
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

# Plot
pdf(file="Figure5E.pdf", width=5, height=5)
ggplot(res, aes(x = logFC, y = sig,label=Taxa)) +
  geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 200 * res$abundance.test, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 200 * res$abundance.control,2)),alpha=0.7) + 
  geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=4,force=25, segment.colour="grey",segment.alpha=0.2) + 
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + 
  xlab("Effect size: log2(fold-change)") + 
  ylab("-log10(adjusted p-value)") + 
  theme 
dev.off()

## DECONTAM ANALYSIS

# Clear all objects
rm(list=ls(all.names=T))
gc()

## Open up the R script with all of the function created for the project
source("KW_functions_version1.R")

# Packages needed for the script
pkg_list <- c("devtools", "tidyverse", "readr", "readtext", "vegan", "ade4", "biomformat", "cachem", 
              "utf8", "backports", "colorspace", "rhdf5", "DelayedArray","Biobase", "Biostrings", "magick",
              "phyloseq", "ape", "phangorn", "ggpubr", "decontam", "ggplot2", "reshape2",
              "tidyr", "matrixStats", "DESeq2", "edgeR", "limma", "Glimma", "RColorBrewer",
              "pheatmap", "ggrepel", "scales", "data.table", "fBasics", "forcats", "maptools", 
              "lubridate", "boot", "table1", "stringr", "papaja", "flextable",
              "Gmisc", "glue", "htmlTable", "grid", "magrittr", "rmarkdown", "plotly",
              "microbiome", "DT", "webshot", "lubridate", "png", "RCurl", "cowplot", "janitor",
              "optmatch", "MatchIt", "decontam", "qdap", "stringr","openxlsx", "chisq.posthoc.test", 
              "FSA", "cobalt", "ggplotify", "grid", "gridExtra", "combinat", "mixOmics", "gplots", "plyr", 
              "readxl", "DESeq2", "mia", "microbiomeMarker", "jpeg", "openxlsx",
              "mia", "miaViz", "corrplot", "ggcorrplot", "cowplot", "gridGraphics",
              "ade4", "ggthemes", "Hmisc", "rdist", "rstatix")

# Install all packages needed  
install_packages(pkg_list)

# Loading packages
for (i in pkg_list){
  eval(bquote(library(.(i))))
}

# Install qiime2R from github directly
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
#devtools::install_github("vmikk/metagMisc")
library(metagMisc)
#devtools::install_github("Sebastien-Le/YesSiR")
library(YesSiR)
#devtools::install_github("thomasp85/patchwork")
library(patchwork)
#devtools::install_github("david-barnett/microViz")
library(microViz)
#devtools::install_github("zdk123/pulsar")
library(pulsar)
#devtools::install_github("zdk123/SpiecEasi")
library(SpiecEasi)
#install.packages("maptools", repos = "https://packagemanager.posit.co/cran/2023-10-13")
library(maptools)
library(openxlsx)

# Run analysis
decontaminant_KW(input_phyloseq=physeq, 
                 sample_type_var_name="Description",                              
                 sample_types=c("bkg", "sputum"), 
                 sample_type_color=c("darkgoldenrod1", "steelblue2"), 
                 sample_type_color_2nd=c("darkgoldenrod4","blue4"), 
                 negative_sample_type="bkg",
                 compare_type=c("sputum"), 
                 stat_option="mean", 
                 display_contam_method="none", 
                 graph_option="boxplot", 
                 test_threshold=0.5,
                 log_scale="yes", 
                 output_suffix="sputum", 
                 taxa_genus_output="no")

# Produces Supplementary Figures 6 and 7