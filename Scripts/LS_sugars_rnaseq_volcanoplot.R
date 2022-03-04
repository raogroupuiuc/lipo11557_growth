
#### set up ####

getwd()
setwd("D:/Box Sync/Manuscripts/01 - Writing - 2021 LS sugars RNAseq/lipo11557_growth/Script")
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


load("20210406.RData")

rm(list = ls()[!ls() %in% c("out.conts")])

head(out.conts)
colnames(out.conts)
# [1] "gene_id"         "proteinId"       "ecNum"           "definition"      "KEGG"            "Amean"           "FC.Cel_vs_Glu"   "rawP.Cel_vs_Glu"
# [9] "FDR.Cel_vs_Glu"  "FC.Xyl_vs_Glu"   "rawP.Xyl_vs_Glu" "FDR.Xyl_vs_Glu"  "FC.Xyl_vs_Cel"   "rawP.Xyl_vs_Cel" "FDR.Xyl_vs_Cel"  "Fstat.ANOVA"    
# [17] "rawP.ANOVA"      "FDR.ANOVA"       "LIPO.Cel.1"      "LIPO.Cel.2"      "LIPO.Cel.3"      "LIPO.Glu.1"      "LIPO.Glu.2"      "LIPO.Glu.3"     
# [25] "LIPO.Xyl.1"      "LIPO.Xyl.2"      "LIPO.Xyl.3" 

head(rownames(out.conts))

celvsglu <- out.conts[,c("FC.Cel_vs_Glu","rawP.Cel_vs_Glu","proteinId","ecNum","definition")]

cel.filter <- celvsglu$FC.Cel_vs_Glu >= 0

celvsglu$logFC[cel.filter] <- log10(celvsglu$FC.Cel_vs_Glu[cel.filter])
celvsglu$logFC[!cel.filter] <- log10(-1/celvsglu$FC.Cel_vs_Glu[!cel.filter])
celvsglu$neg.log10.P <- -log10(celvsglu$rawP.Cel_vs_Glu)

head(celvsglu)
tail(celvsglu)

FC.low <- -2
FC.high <- 2
pval.cut <- 0.05

celvsglu$diffexpressed <- "NO"
celvsglu$diffexpressed[celvsglu$FC.Cel_vs_Glu > FC.high & celvsglu$rawP.Cel_vs_Glu < pval.cut] <- "UP"
celvsglu$diffexpressed[celvsglu$FC.Cel_vs_Glu < FC.low  & celvsglu$rawP.Cel_vs_Glu < pval.cut] <- "DOWN"

celvsglu$labels <- "."

filt <- (celvsglu$logFC > 3 | celvsglu$logFC < -2 | celvsglu$neg.log10.P > 12)
celvsglu$labels[filt] <- celvsglu$proteinId[filt]


xylvsglu <- out.conts[,c("FC.Xyl_vs_Glu","rawP.Xyl_vs_Glu","proteinId","ecNum","definition")]
xyl.filter <- xylvsglu$FC.Xyl_vs_Glu >= 0

xylvsglu$logFC[xyl.filter] <- log10(xylvsglu$FC.Xyl_vs_Glu[xyl.filter])
xylvsglu$logFC[!xyl.filter] <- log10(-1/xylvsglu$FC.Xyl_vs_Glu[!xyl.filter])
xylvsglu$neg.log10.P <- -log10(xylvsglu$rawP.Xyl_vs_Glu)

head(xylvsglu)
tail(xylvsglu)

xylvsglu$diffexpressed <- "NO"
xylvsglu$diffexpressed[xylvsglu$FC.Xyl_vs_Glu > FC.high & xylvsglu$rawP.Xyl_vs_Glu < pval.cut] <- "UP"
xylvsglu$diffexpressed[xylvsglu$FC.Xyl_vs_Glu < FC.low  & xylvsglu$rawP.Xyl_vs_Glu < pval.cut] <- "DOWN"

xylvsglu$labels <- "."

filt <- (xylvsglu$logFC > 2 | xylvsglu$logFC < -2 | xylvsglu$neg.log10.P > 10)
xylvsglu$labels[filt] <- xylvsglu$proteinId[filt]

# plot 
cel.p <- ggplot(data=celvsglu, aes(x=logFC, y=neg.log10.P, col=diffexpressed, label = labels)) + 
  theme_minimal() + geom_point() #+ geom_text()

cel.p2 <- cel.p + geom_vline(xintercept=c(log10(-1/FC.low), log10(FC.high)), col="grey") + 
  geom_hline(yintercept=-log10(pval.cut), col="grey") + labs(col = "Differential expression")

cel.p3 <- cel.p2 + scale_color_manual(values=c("blue", "black", "red")) + xlab("log(Fold Change)") + 
  ylab("-log(pValue)") + ggtitle("Cellobiose vs Glucose") + xlim(-4,4) + ylim(0,20) +
  theme(plot.title = element_text(size=12, face="bold", margin = margin(20, 20, 20, 20)), 
        legend.position = c(1, 1), legend.justification = c("right", "top")) + 
  annotate("text", x = 3, y = 7.5, label = table(celvsglu$diffexpressed)[3], col = "red") +
  annotate("text", x = -2.5, y = 7.5, label = table(celvsglu$diffexpressed)[1], col = "blue")


x11(width = 6, height = 6)
plot(cel.p3)

xyl.p <- ggplot(data=xylvsglu, aes(x=logFC, y=neg.log10.P, col=diffexpressed, label = labels)) + 
  theme_minimal()+ geom_point() #+ geom_text()

xyl.p2 <- xyl.p + geom_vline(xintercept=c(log10(-1/FC.low), log10(FC.high)), col="grey") + 
  geom_hline(yintercept=-log10(pval.cut), col="grey") + labs(col = "Differential expression")

xyl.p3 <- xyl.p2 + scale_color_manual(values=c("blue", "black", "red")) + xlab("log(Fold Change)") + 
  ylab("-log(pValue)") + ggtitle("Xylose vs Glucose") + xlim(-4,4) + ylim(0,20) + 
  theme(plot.title = element_text(size=12, face="bold", margin = margin(20, 20, 20, 20)), 
        legend.position = c(1, 1), legend.justification = c("right", "top")) + 
  annotate("text", x = 3, y = 7.5, label = table(xylvsglu$diffexpressed)[3], col = "red") +
  annotate("text", x = -2, y = 7.5, label = table(xylvsglu$diffexpressed)[1], col = "blue")

x11(width = 6, height = 6)
plot(xyl.p3)


summary(celvsglu$logFC)
summary(xylvsglu$logFC)


table(celvsglu$diffexpressed)
# DOWN   NO   UP 
# 374 6283  851 
table(xylvsglu$diffexpressed) # [1] - down, [3] - up
# DOWN   NO   UP 
# 120 6996  392 
xylvsglu[xylvsglu$logFC == -Inf,]

out.conts[out.conts$proteinId == 64475,]
# no expression in glucose and xylose.

save.image("20210708_Volcanoplot.Rdata")

sessionInfo()
# R version 4.0.5 (2021-03-31)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19042)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    
# 
# attached base packages:
#   [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] PCAtools_2.2.0        ggrepel_0.9.1         WGCNA_1.70-3          fastcluster_1.1.25    dynamicTreeCut_1.63-1 ggplot2_3.3.3        
# [7] sva_3.38.0            BiocParallel_1.24.1   genefilter_1.72.1     mgcv_1.8-34           nlme_3.1-152          rtracklayer_1.49.5   
# [13] GenomicRanges_1.42.0  GenomeInfoDb_1.26.4   IRanges_2.24.1        S4Vectors_0.28.1      Glimma_2.0.0          rgl_0.105.22         
# [19] gplots_3.1.1          affycoretools_1.62.0  Biobase_2.50.0        BiocGenerics_0.36.0   edgeR_3.32.1          limma_3.46.0         
# 
# loaded via a namespace (and not attached):
#   [1] utf8_1.2.1                  R.utils_2.10.1              tidyselect_1.1.0            RSQLite_2.2.5               AnnotationDbi_1.52.0       
# [6] htmlwidgets_1.5.3           grid_4.0.5                  munsell_0.5.0               codetools_0.2-18            preprocessCore_1.52.1      
# [11] miniUI_0.1.1.1              withr_2.4.1                 colorspace_2.0-0            Category_2.56.0             OrganismDbi_1.32.0         
# [16] knitr_1.31                  rstudioapi_0.13             labeling_0.4.2              MatrixGenerics_1.2.1        GenomeInfoDbData_1.2.4     
# [21] hwriter_1.3.2               farver_2.1.0                bit64_4.0.5                 vctrs_0.3.6                 generics_0.1.0             
# [26] xfun_0.22                   biovizBase_1.38.0           BiocFileCache_1.14.0        doParallel_1.0.16           R6_2.5.0                   
# [31] rsvd_1.0.3                  locfit_1.5-9.4              AnnotationFilter_1.14.0     manipulateWidget_0.10.1     bitops_1.0-6               
# [36] cachem_1.0.4                reshape_0.8.8               DelayedArray_0.16.3         assertthat_0.2.1            promises_1.2.0.1           
# [41] scales_1.1.1                nnet_7.3-15                 gtable_0.3.0                beachmat_2.6.4              affy_1.68.0                
# [46] ggbio_1.38.0                ensembldb_2.14.0            rlang_0.4.10                splines_4.0.5               lazyeval_0.2.2             
# [51] impute_1.64.0               dichromat_2.0-0             checkmate_2.0.0             BiocManager_1.30.12         reshape2_1.4.4             
# [56] GenomicFeatures_1.42.2      crosstalk_1.1.1             backports_1.2.1             httpuv_1.5.5                Hmisc_4.5-0                
# [61] RBGL_1.66.0                 tools_4.0.5                 affyio_1.60.0               ellipsis_0.3.1              ff_4.0.4                   
# [66] RColorBrewer_1.1-2          Rcpp_1.0.6                  plyr_1.8.6                  sparseMatrixStats_1.2.1     base64enc_0.1-3            
# [71] progress_1.2.2              zlibbioc_1.36.0             purrr_0.3.4                 RCurl_1.98-1.3              prettyunits_1.1.1          
# [76] rpart_4.1-15                openssl_1.4.3               cowplot_1.1.1               SummarizedExperiment_1.20.0 cluster_2.1.1              
# [81] tinytex_0.31                magrittr_2.0.1              data.table_1.14.0           ProtGenerics_1.22.0         matrixStats_0.58.0         
# [86] hms_1.0.0                   mime_0.10                   xtable_1.8-4                XML_3.99-0.6                jpeg_0.1-8.1               
# [91] gcrma_2.62.0                gridExtra_2.3               compiler_4.0.5              biomaRt_2.46.3              tibble_3.1.0               
# [96] KernSmooth_2.23-18          crayon_1.4.1                ReportingTools_2.30.2       R.oo_1.24.0                 htmltools_0.5.1.1          
# [101] GOstats_2.56.0              later_1.1.0.1               Formula_1.2-4               geneplotter_1.68.0          DBI_1.1.1                  
# [106] dbplyr_2.1.0                rappdirs_0.3.3              Matrix_1.3-2                R.methodsS3_1.8.1           pkgconfig_2.0.3            
# [111] GenomicAlignments_1.26.0    foreign_0.8-81              xml2_1.3.2                  foreach_1.5.1               annotate_1.68.0            
# [116] dqrng_0.2.1                 webshot_0.5.2               XVector_0.30.0              AnnotationForge_1.32.0      stringr_1.4.0              
# [121] VariantAnnotation_1.36.0    digest_0.6.27               graph_1.68.0                Biostrings_2.58.0           htmlTable_2.1.0            
# [126] DelayedMatrixStats_1.12.3   GSEABase_1.52.1             curl_4.3                    shiny_1.6.0                 Rsamtools_2.6.0            
# [131] gtools_3.8.2                lifecycle_1.0.0             jsonlite_1.7.2              PFAM.db_3.12.0              askpass_1.1                
# [136] BSgenome_1.58.0             fansi_0.4.2                 pillar_1.5.1                lattice_0.20-41             GGally_2.1.1               
# [141] fastmap_1.1.0               httr_1.4.2                  survival_3.2-10             GO.db_3.12.1                glue_1.4.2                 
# [146] png_0.1-7                   iterators_1.0.13            bit_4.0.4                   Rgraphviz_2.34.0            stringi_1.5.3              
# [151] blob_1.2.1                  oligoClasses_1.52.0         BiocSingular_1.6.0          DESeq2_1.30.1               latticeExtra_0.6-29        
# [156] caTools_1.18.2              memoise_2.0.0               dplyr_1.0.5                 irlba_2.3.3  
