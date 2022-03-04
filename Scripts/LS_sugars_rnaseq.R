
#### set up ####

getwd()
setwd("D:/Box Sync/Manuscripts/01 - Writing - 2021 LS sugars RNAseq/lipo11557_growth/Scripts")
# 
# x <- c("edgeR", "affycoretools", "limma", "gplots", "rgl", "Glimma", "sva", 
#        "WGCNA")
# BiocManager::install("PCAtools")

options(stringsAsFactors = F)
library(edgeR)
#library(biomaRt)
library(affycoretools)
library(limma)
library(gplots)
library(rgl)
library(Glimma)
library(rtracklayer)
library(sva)
library(ggplot2)
library(WGCNA)
library(PCAtools)

rm(list = ls())

#### Read Targets ####

targets <- readTargets("input/Targets_Final_LIPO.txt")

targets$col <- as.numeric(factor(targets$Trt))
targets$Label <- paste(targets$Species, targets$Trt, targets$Rep, sep=".")
rownames(targets) <- targets$Label


ReadFate <- dplyr::select(targets,QCfiltered:in.a.gene)
rownames(ReadFate) <- rownames(targets)
ReadFate <- ReadFate/targets$Total *100

rowSums(ReadFate)
#all 100

apply(ReadFate, 2, summary)
#           QCfiltered  unmapped multimapped not.in.gene  ambiguous in.a.gene
# Min.    7.062141e-05 0.2139187    0.957263    12.23379 0.05302467  84.79482
# 1st Qu. 1.025106e-04 0.2241951    1.128251    13.16891 0.05573314  85.04495
# Median  1.160725e-04 0.2686598    1.165954    13.28556 0.05743333  85.22861
# Mean    1.347283e-04 0.2738317    1.169015    13.12157 0.05767364  85.37777
# 3rd Qu. 1.859917e-04 0.3038976    1.248807    13.35724 0.05893356  85.41878
# Max.    2.122166e-04 0.3592119    1.339331    13.60661 0.06347608  86.39985

#### read fates plot ####

ReadFate$Sample = targets$Label
df <- tidyr::gather(ReadFate, "Fate", "pct", 1:6)
df$Fate <- factor(df$Fate, levels = unique(df$Fate), ordered = TRUE)

x11(10,6)
ggplot(data=df, aes(x=Sample, y=pct, fill=Fate)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("Read Fate")+ ylab("Percentage")+
  ylim(0, 100)


# x11(w = 10, h = 6)
setEPS()
postscript("results/Readfates.eps", width = 10, height = 6)
par(mar=c(6,4,4,1))
barplot(t(ReadFate[,1:6]), col=1:6,las=2, legend.text = F,beside=T,
        ylab = "Percent of total reads", main = "Read fates per sample", cex.names = 0.7,ylim=c(0,110),
        names.arg = targets$Label, add=F)
abline(h = seq(0,100, by = 20), col = "gray")
barplot(t(ReadFate[,1:6]), col=1:6,las=2, legend.text = F,beside=T,
        ylab = "Percent of total reads", main = "Read fates per sample", cex.names = 0.7,ylim=c(0,110),
        names.arg = targets$Label, add=T)
legend ("topright", legend = colnames( ReadFate[,1:6]), col = 1:6, fill = 1:6, ncol = 6 )


#### Read in Counts ####

d <- readDGE(targets$files, path="input/featcounts/" ,labels=targets$Label, 
             comment.char="#",header=T,columns = c(1,7), group=targets$Trt)
dim(d)
# 8192    9

head(d$counts)
# Tags                                       LIPO.Cel.1 LIPO.Cel.2 LIPO.Cel.3 LIPO.Glu.1 LIPO.Glu.2 LIPO.Glu.3 LIPO.Xyl.1 LIPO.Xyl.2 LIPO.Xyl.3
# fgenesh1_kg.1_-_1_-_Locus9988v2rpkm0.52           8         17         13          3          3          2          3          6          2
# fgenesh1_kg.1_-_4_-_Locus7012v1rpkm5.68         356        201        104        662        460        646        848        580        722
# CE30_10333                                        0          0          0          0          0          0          1          2          1
# estExt_Genewise1Plus.C_1_t10002                9951       2433        290      35809      33448      29821      21058      14942      22397
# CE166_153                                        27         21         16         57         80         68         66         51         50
# fgenesh1_kg.1_-_11_-_Locus6662v1rpkm7.09        418        602        476        932        683        706        754        593        701

tail(d$counts)
# Tags                                            LIPO.Cel.1 LIPO.Cel.2 LIPO.Cel.3 LIPO.Glu.1 LIPO.Glu.2 LIPO.Glu.3 LIPO.Xyl.1 LIPO.Xyl.2 LIPO.Xyl.3
# MIX6738_3604_100                                       4         14          9          0          5          1          2          2          1
# fgenesh1_kg.96_-_1_-_Locus692v12rpkm1.03_PRE         145        245        223        102        157        143         98         80         97
# fgenesh1_kg.97_-_2_-_Locus5497v13rpkm0.46_PRE          0          0          0          0          0          0          0          0          0
# gm1.7515_g                                             0          5          3          0          1          4          1          4          0
# fgenesh1_kg.98_-_1_-_Locus13884v2rpkm0.38              2          0          0          0          0          2          0          0          0
# fgenesh1_kg.99_-_4_-_Locus8153v2rpkm0.00_PRE           0          0          0          0          0          0          0          0          0
# 


#### Add gtf annot ####

gtf0 <- import("input/genome/Lipst1_1_GeneCatalog_genes_20110609.gtf")
gtf0
#GRanges object with 61632 ranges and 10 metadata columns:

table(gtf0$type)
# exon         CDS start_codon  stop_codon 
# 23318       22698        7818        7798

gtf1 <- gtf0[gtf0$type == "exon"]
gtf1
#GRanges object with 23318 ranges and 10 metadata columns:

length(unique(gtf1$gene_id))
#8192
length(unique(gtf1$transcript_id))
#8192

gtf2 <- gtf1[!duplicated(gtf1$gene_id)]
gtf2
#GRanges object with 8192 ranges and 10 metadata columns:

sum(is.na(gtf2$transcriptId))
#0

table(gtf2$source)
# JGI 
# 8192 

#also get protein IDs, which are in the CDS rows

gtf3 <- gtf0[gtf0$type == "CDS"]
gtf3 <- gtf3[!duplicated(gtf3$gene_id)]
gtf3

all.equal(gtf2$gene_id, gtf3$gene_id)
#TRUE

gtf2$proteinId <- gtf3$proteinId

#keep only gene_id and proteinId

all.equal(gtf2$gene_id, rownames(d$counts))
#TRUE

d$genes <- mcols(gtf2)[,c("gene_id", "proteinId")]

head(d$gene)


#### check KEGG and EC annotations ####

ec.kegg <- read.delim("input/genome/Lipst1_1_GeneCatalog_proteins_20110609_KEGG.tab")
dim(ec.kegg)
#3741    9
head(ec.kegg)

#check to see how many protein Ids from the annotation file are in the gtf 

sum(ec.kegg$X.proteinId %in% d$genes$proteinId)
#3741 - all protein IDs in annotation found in gtf

sum(!ec.kegg$X.proteinId %in% d$genes$proteinId)
#0

#how many genes have ec or kegg?

sum(d$genes$proteinId %in% ec.kegg$X.proteinId)
#1819  - not very many

sum(duplicated(ec.kegg$X.proteinId))
#1922 - multiple entries in annotation for some genes

View(ec.kegg)

#Remove any leading or trailing white spaces in the ecNum, definition or pathway

ec.kegg$ecNum <- gsub("^\\s+|\\s+$", "", ec.kegg$ecNum) 
ec.kegg$definition <- gsub("^\\s+|\\s+$", "", ec.kegg$definition) 
ec.kegg$pathway <- gsub("^\\s+|\\s+$", "", ec.kegg$pathway) 
ec.kegg$pathway_class <- gsub("^\\s+|\\s+$", "", ec.kegg$pathway_class) 
ec.kegg$pathway_type <- gsub("^\\s+|\\s+$", "", ec.kegg$pathway_type) 



table(ec.kegg$pathway_type)
#       \\N  METABOLIC REGULATORY 
#       566       3044        131

table(ec.kegg$pathway_class)


#separate out and reduce down

ec.part <- ec.kegg[,c(1:3)]
dim(ec.part)
#3741    3
ec.part <- ec.part[!duplicated(paste(ec.part$X.proteinId, ec.kegg$ecNum)),]
dim(ec.part)
#1819    3

sum(duplicated (ec.part$X.proteinId))
#0

ec.part$X.proteinId <- as.character(ec.part$X.proteinId)
  
d$genes <- as.data.frame(d$genes)
d$genes <- dplyr::left_join(d$genes, ec.part, by = c("proteinId" = "X.proteinId"))

head(d$genes)
sum(!is.na(d$genes$ecNum))
#1819

head(ec.kegg)

kegg.part <- ec.kegg[,c(1, 7:9)]
head(kegg.part)

sum(is.na(kegg.part$pathway))
#0
sum(kegg.part$pathway == "\\N")
#537


#remove these
kegg.part <- kegg.part[kegg.part$pathway != "\\N",]
dim(kegg.part)
#3204    4

#any protein-pathway duplicates?

sum(duplicated(paste(kegg.part$X.proteinId, kegg.part$pathway)))
#0

table(table(kegg.part$X.proteinId))
#  1   2   3   4   5   6   7   8   9  11  12  13  14  15 
#626 323 148  41  45   9   1  13  11  19  27  11   6   5


#collapse down multiple pathways per protein; check to see if any ; already in pathway names

length(grep(";", kegg.part$pathway))
# 0  - ok to use ; as separator

kegg.part2 <- tapply(kegg.part$pathway, list(kegg.part$X.proteinId), function(x) paste(unique(x), collapse = "; "))
head(kegg.part2)


d$genes$KEGG <- kegg.part2[d$genes$proteinId]
head(d$genes)
sum(!is.na(d$genes$KEGG))
#1285

#some proteins have EC but no KEGG pathways



#remove to save space

rm(gtf0, gtf1, gtf2, gtf3)


save.image("20210406.RData")


#checking for highest counts

range(rowSums(d$counts))
# 0 1247512
max(rowSums(d$counts))/sum(d$counts)
# 0.009955756

#maxium count percentage of each library
x11()
apply(d$counts,2,function(x){max(x)/sum(x)}) %>%
  barplot(las=2,col=targets$col,main="maxium count percentage of each library", cex.axis = 0.8)
abline(h=mean(apply(d$counts,2,function(x){max(x)/sum(x)})))
#Cel lowest, Xyl highest


#### library size ####

x11(20,20)
#jpeg("library_size.jpeg", width = 7, height = 4, units = "in", res = 300, quality = 100)
barplot(d$samples$lib.size/1e6,las=2, main = "Number of reads in genes (individual lib size)",
        ylab = "millions of reads", col=targets$col,names.arg = targets$Label, cex.names = 0.9)
abline(h=mean(d$samples$lib.size/1e6),col="blue")
#dev.off()
#variable but OK


# genes with no reads in all samples  
d$counts[apply(d$counts,1,sum)==0,] %>%nrow()
#  87
# total gene number:
nrow(d$counts)
#8192

# genes with reads in all samples  
d$counts[apply(d$counts,1,sum)!=0,] %>%nrow()
# 8105


# detecting trends between library size and gene with 0 counts
x11()
plot(d$samples$lib.size/1e6,apply(d$counts==0,2,sum), 
     xlab="Library Size (million)", ylab="number genes of zero count", 
     col=targets$col, pch=16)
#NO correlation, but Cel has fewer genes with 0 counts


#### Calculate TMM ####

d <- calcNormFactors(d) 
x11(20,15)
barplot(d$samples$norm.factors,col=targets$col,las=2,names.arg = targets$Label, main="normalization factor",cex.names=0.8)
abline(h=1, lty=2)
#Xyl has lower values, meaning their lib sizes are too large

summary(d$samples$lib.size/1e6)


#calculate counts per million of reads per kilobase per million
logCPM <- cpm(d, log=T, prior.count = 3)

#CPM density plot
x11(20,20)
plotDensities(logCPM,group=d$samples$group,col=as.numeric(unique(targets$col)),legend="topright")


#who has more unique expressing genes?
x11(10,6)
unique_index <- rowSums(logCPM>log2(0.5))==1
barplot(colSums((logCPM>log2(0.5))[unique_index,]), col = targets$col, names.arg = targets$Label, 
        las=2, cex.names=0.8, main=paste0("number of unique genes out of ", nrow(logCPM)))



#frequency of number of samples larger than 1 cpm
x11()
hist(rowSums(logCPM > log2(1)))
# the actual number
table(rowSums(logCPM > log2(1)))
#   0    1    2    3    4    5    6    7    8    9 
# 535   85   64   70   46   54   66   70  105 7097 

# how many genes has 1 cpms in at least 3 samples?
sum(rowSums(logCPM>log2(1))>=3)
#7508

#the precentage
mean(rowSums(logCPM>log2(1))>=3)
#0.9165039

#### Filtering ####

# use 0.5 cpm in at least 4 samples as a criteria
filter_index <- rowSums(logCPM>log2(1))>=3

#filtering d object
d.filt <- d[filter_index,,keep.lib.size=F]
dim(d)
#8192    9
dim(d.filt)
#7508    9


# make sure libs was not messed up
range(d.filt$samples$lib.size/d$samples$lib.size)
# 0.9997932 0.9999144

#redo TMM normalization
d.filt <- calcNormFactors(d.filt)  

x11()
barplot(d.filt$samples$norm.factors,col=targets$col,las=2,names.arg = targets$Label, main="filterd normalization factor",cex.names=0.9)
abline(h=1, lty=2)
#same, Xyl lower

all.equal(d.filt$samples$norm.factors,d$samples$norm.factors)
#"Mean relative difference: 0.002573549"

x11()
plot(d.filt$samples$norm.factors,d$samples$norm.factors, col=targets$col)
abline(0,1)
cor(d.filt$samples$norm.factors,d$samples$norm.factors)
#  0.9993071


#recalculate CPM ####
logCPM.filt <- cpm(d.filt, log = T, prior.count = 3)

#
#density plot
x11()
plotDensities(logCPM.filt,group=d$samples$group,col=as.numeric(unique(targets$col)) )
#Cel has higher peak > 5


# clustering analysis
x11(width = 20 , height = 10 )
hclust( dist( t(logCPM.filt ) ), method = "average" ) %>%
  plot( hang = -1, main = "CPM after filtering", sub = "", xlab = "", cex = 0.9 )
# no issue

# PCA ####
x11( w =5, h = 6)
test_new <- pca(logCPM.filt)
postscript("results/PCA.eps", width = 5, height = 6)
biplot(test_new, xlim = c(-75,50), ylim = c(-50,50), 
       lab = d.filt$samples$group, 
       pointSize = 4, gridlines.major = FALSE, gridlines.minor = FALSE)
dev.off()


#### Design matrix ####

#one-way anova
group <- factor(d$samples$group, levels=c("Glu", "Cel", "Xyl"))
design <- model.matrix(~group)
colnames(design)[-1] <- paste0(levels(group)[-1], "_vs_", levels(group)[1])
design
#   (Intercept) Cel_vs_Glu Xyl_vs_Glu
# 1           1          1          0
# 2           1          1          0
# 3           1          1          0
# 4           1          0          0
# 5           1          0          0
# 6           1          0          0
# 7           1          0          1
# 8           1          0          1
# 9           1          0          1


fit <- lmFit(logCPM.filt, design = design)

#Add in genes info

fit$genes <- d.filt$genes

cont.matrix <- makeContrasts(Cel_vs_Glu,
                             Xyl_vs_Glu,
                             Xyl_vs_Cel = Xyl_vs_Glu - Cel_vs_Glu,
                             levels = design)
fit2 <- contrasts.fit(fit, contrasts = cont.matrix)
fit2 <- eBayes(fit2, trend = T)
summary(decideTests(fit2))
#        Cel_vs_Glu Xyl_vs_Glu Xyl_vs_Cel
# Down         1504        779       1589
# NotSig       4341       5822       4389
# Up           1663        907       1530


#### 1way test ####

res.anova <- topTable(fit2, coef = 1:2, num = Inf, sort.by = "none" )
sum(res.anova$adj.P.Val < 0.05)
#4089

Sig_index <- topTable(fit2, coef = 1:2, n=Inf, sort.by = "none")$adj.P.Val < 0.05
table(Sig_index)
  
# Sig_index
# FALSE  TRUE 
#  3419  4089


source("input/summarizeFit.R")

#### Get pairwise ####

out.conts <- summarizeFit(fit2, addAnova = TRUE)
names(out.conts)
# [1] "gene_id"         "proteinId"       "ecNum"           "definition"     
# [5] "KEGG"            "Amean"           "FC.Cel_vs_Glu"   "rawP.Cel_vs_Glu"
# [9] "FDR.Cel_vs_Glu"  "FC.Xyl_vs_Glu"   "rawP.Xyl_vs_Glu" "FDR.Xyl_vs_Glu" 
# [13] "FC.Xyl_vs_Cel"   "rawP.Xyl_vs_Cel" "FDR.Xyl_vs_Cel"  "F"              
# [17] "F.rawP"          "F.FDR"       

sum(out.conts$F.FDR < 0.05)
#4089  - the correct numberr

colnames(out.conts)[16:18] <- c("Fstat.ANOVA","rawP.ANOVA","FDR.ANOVA")


##heatmap ####

logCPM.heat <- t(scale(t(logCPM.filt[Sig_index,])))
dim(logCPM.heat)
col.pan <- colorpanel(100, "blue", "white", "red")

x.cluster.h1 <- hclust(dist(logCPM.heat))

# x11(10,15)

setEPS()
postscript("results/Heatmap.eps", width = 6, height = 10)
heatmap.2(logCPM.heat, col=col.pan, Rowv=as.dendrogram(x.cluster.h1),Colv=FALSE, scale="none",
          trace="none", dendrogram="row", cexRow=1, cexCol=0.9,labRow ="",
          keysize = 0.8,lwid = c(1.4,4), lhei = c(0.8,4), colsep = c(3,6), sepcolor = "black",
          margins = c(7,2),density.info = "none", 
          main = paste0("1-way ANOVA\nFDR < 0.05","\n",nrow(logCPM.heat), " genes" ))
dev.off()

x11(12,6)
plot(x.cluster.h1, labels = FALSE)
abline(h = 4.5 )

rowcols.h1 <- cutree(x.cluster.h1, h = 4.5)
table(rowcols.h1)
#  1   2   3   4   5   6   7 
#999 552 426 896 318 368 530


#put in results and add on individual values

names(out.conts)
# [1] "gene_id"         "proteinId"       "ecNum"           "definition"     
# [5] "KEGG"            "Amean"           "FC.Cel_vs_Glu"   "rawP.Cel_vs_Glu"
# [9] "FDR.Cel_vs_Glu"  "FC.Xyl_vs_Glu"   "rawP.Xyl_vs_Glu" "FDR.Xyl_vs_Glu" 
# [13] "FC.Xyl_vs_Cel"   "rawP.Xyl_vs_Cel" "FDR.Xyl_vs_Cel"  "Fstat.ANOVA"    
# [17] "rawP.ANOVA"      "FDR.ANOVA"

all.equal(out.conts$gene_id, rownames(logCPM.filt))
#TRUE
out.conts <- cbind(out.conts, logCPM.filt)
names(out.conts)
# [1] "gene_id"         "proteinId"       "ecNum"           "definition"     
# [5] "KEGG"            "Amean"           "FC.Cel_vs_Glu"   "rawP.Cel_vs_Glu"
# [9] "FDR.Cel_vs_Glu"  "FC.Xyl_vs_Glu"   "rawP.Xyl_vs_Glu" "FDR.Xyl_vs_Glu" 
# [13] "FC.Xyl_vs_Cel"   "rawP.Xyl_vs_Cel" "FDR.Xyl_vs_Cel"  "Fstat.ANOVA"    
# [17] "rawP.ANOVA"      "FDR.ANOVA"       "LIPO.Cel.1"      "LIPO.Cel.2"     
# [21] "LIPO.Cel.3"      "LIPO.Glu.1"      "LIPO.Glu.2"      "LIPO.Glu.3"     
# [25] "LIPO.Xyl.1"      "LIPO.Xyl.2"      "LIPO.Xyl.3"

write.table(out.conts[,c(1:6, 16:18, 7:15, 19:27)], 
            file = "results/gene_expression.txt", row.names = FALSE, sep = "\t")




save.image("20210406.RData")



sessionInfo()
# R version 4.0.5 (2021-03-31)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19042)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
# [3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.1252    
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods  
# [9] base     
# 
# other attached packages:
#   [1] PCAtools_2.2.0        ggrepel_0.9.1         WGCNA_1.70-3          fastcluster_1.1.25   
# [5] dynamicTreeCut_1.63-1 sva_3.38.0            BiocParallel_1.24.1   genefilter_1.72.1    
# [9] mgcv_1.8-34           nlme_3.1-152          Glimma_2.0.0          rgl_0.105.22         
# [13] gplots_3.1.1          affycoretools_1.62.0  Biobase_2.50.0        edgeR_3.32.1         
# [17] limma_3.46.0          ggplot2_3.3.3         rtracklayer_1.49.5    GenomicRanges_1.42.0 
# [21] GenomeInfoDb_1.26.4   IRanges_2.24.1        S4Vectors_0.28.1      BiocGenerics_0.36.0  
# 
# loaded via a namespace (and not attached):
#   [1] utf8_1.2.1                  R.utils_2.10.1              tidyselect_1.1.0           
# [4] RSQLite_2.2.5               AnnotationDbi_1.52.0        htmlwidgets_1.5.3          
# [7] grid_4.0.5                  munsell_0.5.0               codetools_0.2-18           
# [10] preprocessCore_1.52.1       miniUI_0.1.1.1              withr_2.4.1                
# [13] colorspace_2.0-0            Category_2.56.0             OrganismDbi_1.32.0         
# [16] knitr_1.31                  rstudioapi_0.13             labeling_0.4.2             
# [19] MatrixGenerics_1.2.1        GenomeInfoDbData_1.2.4      hwriter_1.3.2              
# [22] farver_2.1.0                bit64_4.0.5                 vctrs_0.3.6                
# [25] generics_0.1.0              xfun_0.22                   biovizBase_1.38.0          
# [28] BiocFileCache_1.14.0        doParallel_1.0.16           R6_2.5.0                   
# [31] rsvd_1.0.3                  locfit_1.5-9.4              AnnotationFilter_1.14.0    
# [34] manipulateWidget_0.10.1     bitops_1.0-6                cachem_1.0.4               
# [37] reshape_0.8.8               DelayedArray_0.16.3         assertthat_0.2.1           
# [40] promises_1.2.0.1            scales_1.1.1                nnet_7.3-15                
# [43] gtable_0.3.0                beachmat_2.6.4              affy_1.68.0                
# [46] ggbio_1.38.0                ensembldb_2.14.0            rlang_0.4.10               
# [49] splines_4.0.5               lazyeval_0.2.2              impute_1.64.0              
# [52] dichromat_2.0-0             checkmate_2.0.0             BiocManager_1.30.12        
# [55] reshape2_1.4.4              GenomicFeatures_1.42.2      crosstalk_1.1.1            
# [58] backports_1.2.1             httpuv_1.5.5                Hmisc_4.5-0                
# [61] RBGL_1.66.0                 tools_4.0.5                 affyio_1.60.0              
# [64] ellipsis_0.3.1              ff_4.0.4                    RColorBrewer_1.1-2         
# [67] Rcpp_1.0.6                  plyr_1.8.6                  sparseMatrixStats_1.2.1    
# [70] base64enc_0.1-3             progress_1.2.2              zlibbioc_1.36.0            
# [73] purrr_0.3.4                 RCurl_1.98-1.3              prettyunits_1.1.1          
# [76] rpart_4.1-15                openssl_1.4.3               cowplot_1.1.1              
# [79] SummarizedExperiment_1.20.0 cluster_2.1.1               tinytex_0.31               
# [82] magrittr_2.0.1              data.table_1.14.0           ProtGenerics_1.22.0        
# [85] matrixStats_0.58.0          hms_1.0.0                   mime_0.10                  
# [88] xtable_1.8-4                XML_3.99-0.6                jpeg_0.1-8.1               
# [91] gcrma_2.62.0                gridExtra_2.3               compiler_4.0.5             
# [94] biomaRt_2.46.3              tibble_3.1.0                KernSmooth_2.23-18         
# [97] crayon_1.4.1                ReportingTools_2.30.2       R.oo_1.24.0                
# [100] htmltools_0.5.1.1           GOstats_2.56.0              later_1.1.0.1              
# [103] Formula_1.2-4               tidyr_1.1.3                 geneplotter_1.68.0         
# [106] DBI_1.1.1                   dbplyr_2.1.0                rappdirs_0.3.3             
# [109] Matrix_1.3-2                R.methodsS3_1.8.1           pkgconfig_2.0.3            
# [112] GenomicAlignments_1.26.0    foreign_0.8-81              xml2_1.3.2                 
# [115] foreach_1.5.1               annotate_1.68.0             dqrng_0.2.1                
# [118] webshot_0.5.2               XVector_0.30.0              AnnotationForge_1.32.0     
# [121] stringr_1.4.0               VariantAnnotation_1.36.0    digest_0.6.27              
# [124] graph_1.68.0                Biostrings_2.58.0           htmlTable_2.1.0            
# [127] DelayedMatrixStats_1.12.3   GSEABase_1.52.1             curl_4.3                   
# [130] shiny_1.6.0                 Rsamtools_2.6.0             gtools_3.8.2               
# [133] lifecycle_1.0.0             jsonlite_1.7.2              PFAM.db_3.12.0             
# [136] viridisLite_0.3.0           askpass_1.1                 BSgenome_1.58.0            
# [139] fansi_0.4.2                 pillar_1.5.1                lattice_0.20-41            
# [142] GGally_2.1.1                fastmap_1.1.0               httr_1.4.2                 
# [145] survival_3.2-10             GO.db_3.12.1                glue_1.4.2                 
# [148] png_0.1-7                   iterators_1.0.13            bit_4.0.4                  
# [151] Rgraphviz_2.34.0            stringi_1.5.3               blob_1.2.1                 
# [154] oligoClasses_1.52.0         BiocSingular_1.6.0          DESeq2_1.30.1              
# [157] latticeExtra_0.6-29         caTools_1.18.2              memoise_2.0.0              
# [160] dplyr_1.0.5                 irlba_2.3.3   
