#################################################################

                  ### Differential abundance ###

#################################################################


da_wilcox <- function(ps.x=NULL, prev=NULL, agglom=NULL, variable=NULL, groups=NULL, mc_samples=1000, normalise=c("TSS","CLR","rare"), rnd.seed=101, rare_depth=NULL){
  library(zCompositions)
  library(compositions)
  #otu_table orientation
  if(taxa_are_rows(ps.x)){
    otu_table(ps.x) <- t(otu_table(ps.x))
  }
  
  if(length(variable) > 1){
    stop("Error - only one variable may be analysed at a time")
  }
  if(length(groups) != 2){
    stop("Error - 'groups' must contain exactly 2 conditions")
  }
  if(normalise == "rare"){
    if(!is.null(rare_depth)){
      ps.rare <- rarefy_even_depth(ps.x, sample.size = rare_depth, rngseed = rnd.seed)
      otus <- assign_taxa_names(ps.x=ps.rare, agglom=agglom)
      #prevalence filtering
      if(!is.null(prev)){
        otus <- otus[,colSums(otus == 0)/nrow(otus) <= prev]
      }
    } else {
      if(min(sample_sums(ps.x)) > 2000){
        ps.rare <- rarefy_even_depth(ps.x, rngseed = rnd.seed)
        otus <- assign_taxa_names(ps.x=ps.rare, agglom=agglom)
        #prevalence filtering
        if(!is.null(prev)){
          otus <- otus[,colSums(otus == 0)/nrow(otus) <= prev]
        }
      } else {
        stop("Minimum sample depth after filtering is < 2000 reads. Please remove prior to rarefaction.")
      }
    }
  } else {
      #assign taxa names
      otus <- assign_taxa_names(ps.x=ps.x, agglom=agglom)
      #prevalence filtering
      if(!is.null(prev)){
        otus <- otus[,colSums(otus == 0)/nrow(otus) <= prev]
      }
      if(normalise == "TSS"){
        otus <- sweep(otus, 1, rowSums(otus), "/")
      } else if(normalise == "CLR"){
        otus <- cmultRepl(otus, method = "CZM")
        otus <- clr(otus)
      } else {
        stop("Error - 'normalise' must be one of 'TSS' or 'CLR'")
      }
    }
  if(sum(rowSums(otus) == 0) == ncol(otus)){
    stop("Error - sample with no reads following filtering. Please check filtering and read count thresholds.")
  }
  meta_data <- data.frame(sample_data(ps.x))
  otus <- otus[meta_data[,variable] %in% groups,]
  meta_data_s <- meta_data[meta_data[,variable] %in% groups,]
  if(!all(rownames(meta_data_s) == rownames(otus))){
    stop("Error - sample names do not match")
  }
  meta_data_s[,variable] <- factor(as.character(meta_data_s[,variable]), levels = rev(groups)) # reverse groups to keep sign of output consistent with ALDEx2/DESeq2
  p.val.vec <- apply(otus, 2, function(x) wilcox.test(x ~ meta_data_s[,variable], conf.int=TRUE)$p.value)
  estimate.val.vec <- apply(otus, 2, function(x) wilcox.test(x ~ meta_data_s[,variable], conf.int=TRUE)$estimate)
  wilcox.df <- data.frame(estimate=estimate.val.vec, pval=p.val.vec, padj=p.adjust(p.val.vec, "fdr"))
  return(wilcox.df)
  
  
}



da_aldex2 <- function(ps.x=NULL, prev=NULL, agglom=NULL, variable=NULL, groups=NULL, mc_samples=1000, rnd.seed=101){
  library(ALDEx2)
  #otu_table orientation
  if(taxa_are_rows(ps.x)){
    otu_table(ps.x) <- t(otu_table(ps.x))
  }
  
  #assign taxa names
  otus <- assign_taxa_names(ps.x=ps.x, agglom=agglom)
  if(sum(rowSums(otus) == 0) > 0){
    stop("Error - sample(s) with no reads following filtering. Please check filtering and read count thresholds.")
  }
  #prevalence filtering
  if(!is.null(prev)){
    otus <- otus[,colSums(otus == 0)/nrow(otus) <= prev]
  }
  otus <- t(otus)
  meta_data <- data.frame(sample_data(ps.x))
  otus <- otus[,meta_data[,variable] %in% groups]
  meta_data_s <- meta_data[meta_data[,variable] %in% groups,]
  if(!all(rownames(meta_data_s) == colnames(otus))){
    stop("Error - sample names do not match")
  }
  meta_data_s[,variable] <- factor(as.character(meta_data_s[,variable]), levels = groups)
  set.seed(rnd.seed)
  x <- aldex.clr(otus, meta_data_s[,variable], mc.samples=mc_samples, denom="all", verbose=F)
  x.tt <- aldex.ttest(x, paired.test=FALSE, verbose=FALSE)
  x.effect <- aldex.effect(x, CI=T, verbose=FALSE)
  x.combined <- data.frame(x.tt,x.effect)
  return(x.combined)
}





da_deseq2 <- function(ps.x=NULL, prev=NULL, agglom=NULL, variable=NULL, groups=NULL, mc_samples=1000, rnd.seed=101){
  library(DESeq2)
  #agglomeration
  if(!is.null(agglom)){
    ps.x <- tax_glom(ps.x, agglom, NArm = FALSE)
  }
  #prevalence filtering
  if(!is.null(prev)){
    ps.x <- prune_taxa(colSums(otu_table(ps.x) == 0)/nrow(otu_table(ps.x)) <= prev, ps.x)
  }
  #add taxa labels from 'clean_taxa_names()' to keep consistent across DA techniques
  tax_table(ps.x) <- tax_table(cbind(tax_table(ps.x),Taxa=clean_taxa_names(ps.x)))
  #convert 'NAs' to NA.x
  if(sum(is.na(sample_data(ps.x)[,variable])) > 0){
    temp.vec <- as.vector(unlist(sample_data(ps.x)[,variable]))
    temp.vec[is.na(temp.vec)] <- "NA.x"
    sample_data(ps.x)[,variable] <- temp.vec
  }
  
  #Run DESeq2
  design.f <- paste("~",eval(variable), sep = " ")
  design.f2 <- formula(design.f)
  deseq2 = phyloseq_to_deseq2(ps.x, design = design.f2)
  
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  deseq2_geoMeans = apply(counts(deseq2), 1, gm_mean)
  deseq2 = estimateSizeFactors(deseq2, geoMeans = deseq2_geoMeans)
  deseq2 = DESeq(deseq2, fitType="local")
  #Extract relevant results (N0 versus Nx)
  ax <- variable
  bx <- groups[2]
  cx <- groups[1]
  prac <- c(ax, bx, cx)
  res = results(deseq2, contrast=c(ax, bx, cx))
  res = cbind(as(res, "data.frame"), as(tax_table(ps.x)[rownames(res), ], "matrix"))
  
  #Order results at chosen level of (FDR-adjusted) significance
  res = res[order(res$padj, na.last=NA), ]
  #remove sequences as rownames
  rownames(res) <- res$Taxa
  return(res)
}



da_ensemble <- function(ps.x=NULL, prev=NULL, agglom=NULL, variable=NULL, groups=NULL, mc_samples=1000, normalise=c("TSS","CLR","rare"), rnd.seed=101, rare_depth=NULL){
  library(pheatmap)
  library(RColorBrewer)
  library(ggplot2)
  res.list <- vector("list",5)
  names(res.list) <- c("Wilcoxon", "ALDEx2", "DESeq2", "heatmap", "effectplot")
  da.test1 <- da_wilcox(ps.x, prev = prev, agglom = agglom, variable = variable, groups = groups, normalise = normalise)
  da.test2 <- da_aldex2(ps.x, prev = prev, agglom = agglom, variable = variable, groups = groups)
  da.test3 <- da_deseq2(ps.x, prev = prev, agglom = agglom, variable = variable, groups = groups)
  res.list$Wilcoxon <- da.test1
  res.list$ALDEx2 <- da.test2
  res.list$DESeq2 <- da.test3  
  da.test1 <- da.test1[da.test1$padj < 0.2, ]
  da.test2 <- da.test2[da.test2$wi.eBH < 0.2, ]
  da.test3 <- da.test3[da.test3$padj < 0.2, ]
  
  #plot consensus as heatmap
  da.test.df <- data.frame(matrix(nrow = length(unique(c(rownames(da.test1), rownames(da.test2), rownames(da.test3)))), ncol = 3))
  rownames(da.test.df) <- unique(c(rownames(da.test1), rownames(da.test2), rownames(da.test3)))
  colnames(da.test.df) <- c("Wilcoxon", "ALDEx2", "DESeq2")
  
  da.test.df[rownames(da.test1), 1] <- da.test1$estimate
  da.test.df[rownames(da.test2), 2] <- da.test2$effect
  da.test.df[rownames(da.test3), 3] <- da.test3$log2FoldChange
  da.test.df[is.na(da.test.df)] <- 0
  da.test.df[da.test.df < 0] <- -1
  da.test.df[da.test.df > 0] <- 1
  if(all(da.test.df == 0)){
    res.list$heatmap <- "No DA taxa by any method"
  } else {
    if(min(da.test.df) == 0){
      res.list$heatmap <- pheatmap(as.matrix(da.test.df), cluster_cols = F, cluster_rows = F, color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                                                                        "Reds")))(100), silent = T)
    } else if(max(da.test.df) == 0){
      res.list$heatmap <- pheatmap(as.matrix(da.test.df), cluster_cols = F, cluster_rows = F, color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                                                                        "Blues")))(100), silent = T)
    } else{
      res.list$heatmap <- pheatmap(as.matrix(da.test.df), cluster_cols = F, cluster_rows = F, color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                                                                        "RdBu")))(100), silent = T)
    }
  }
  
  

  
  da.taxa.intersect <- rownames(da.test.df[abs(rowSums(da.test.df)) == 3,])
  if(length(da.taxa.intersect) == 0){
    res.list$effectplot <- "No DA taxa detected by all methods"
  } else {
    sigtab <- da.test3[da.test3$Taxa %in% da.taxa.intersect,]
    sigtab <- sigtab[order(sigtab$log2FoldChange),]
    sigtab$Taxa = factor(sigtab$Taxa, levels = unique(sigtab$Taxa))
    c1 <- groups[1]
    c2 <- groups[2]
    sigtab$diff <- unlist(lapply(sigtab$log2FoldChange, function(y) if(y > 0){c2} else {c1}))
    sigtab$FDR <- unlist(lapply(sigtab$padj, function(y) if(y < 0.05){"[fdr<0.05]"} else {"[fdr<0.2]"}))
    sigtab$fill <- paste(sigtab$diff, sigtab$FDR, sep = " ")
    
    c1.leg.05 <- paste(c1, "[fdr<0.05]", sep = " ")
    c2.leg.05 <- paste(c2, "[fdr<0.05]", sep = " ")
    c1.leg.2 <- paste(c1, "[fdr<0.2]", sep = " ")
    c2.leg.2 <- paste(c2, "[fdr<0.2]", sep = " ")
    sigtab$fill <- factor(sigtab$fill, levels = c(c1.leg.05, c2.leg.05, c1.leg.2, c2.leg.2))
    res.list$effectplot <- ggplot(sigtab, aes(x=log2FoldChange, y=Taxa, fill=fill)) +
      geom_bar(stat = "identity",colour="black", width = 0.7) + scale_y_discrete(limits=sigtab$Taxa) + theme_bw() +
      theme(legend.position = "none", axis.title.y =  element_blank(), text = element_text(size = 18), 
            plot.title = element_text(hjust = 0.5)) + coord_fixed(ratio = 0.3) + ggtitle(paste(c1, c2, sep = " vs ")) + xlab("Log2-fold change")
  }
  
  return(res.list)
  
}




cancer_ibd <- da_ensemble(ps2, prev = 0.9, agglom = "Genus", variable = "Neoplasia3", groups = c("IBD.N0", "IBD.N3"), normalise = "TSS")

pdf("ensemble_DA/heatmap_IBD_N3_N0.pdf", width = 3.1, height = 2.3)
cancer_ibd$heatmap
dev.off()

pdf("ensemble_DA/effect_IBD_N3_N0.pdf", width = 7, height = 5)
cancer_ibd$effectplot + xlim(c(-8,1.5)) + scale_fill_manual("Legend",values = c("IBD.N0 [fdr<0.05]"= alpha("dodgerblue",1),
                                                                                "IBD.Nx [fdr<0.05]"= alpha("firebrick",1),
                                                                                "IBD.N0 [fdr<0.2]"= alpha("dodgerblue",0.5),
                                                                                "IBD.Nx [fdr<0.2]"= alpha("firebrick",0.5)))
dev.off()



subtype_cntrl_uc <- da_ensemble(ps2, prev = 0.9, agglom = "Genus", variable = "Subtype", groups = c("Control", "UC"), normalise = "TSS")

pdf("ensemble_DA/heatmap_subtype_cntrl_uc.pdf", width = 3.3, height = 2.3)
subtype_cntrl_uc$heatmap
dev.off()

pdf("ensemble_DA/effect_subtype_cntrl_uc.pdf", width = 5, height = 3.5)
subtype_cntrl_uc$effectplot + xlim(c(-8,1.5)) + scale_fill_manual("Legend",values = c("Control [fdr<0.05]"= alpha("dodgerblue",1),
                                                                                      "UC [fdr<0.05]"= alpha("firebrick",1),
                                                                                      "Control [fdr<0.2]"= alpha("dodgerblue",0.5),
                                                                                      "UC [fdr<0.2]"= alpha("firebrick",0.5)))
dev.off()



subtype_cntrl_cd <- da_ensemble(ps2, prev = 0.9, agglom = "Genus", variable = "Subtype", groups = c("Control", "CD"), normalise = "TSS")

pdf("ensemble_DA/heatmap_subtype_cntrl_cd.pdf", width = 4, height = 5)
subtype_cntrl_cd$heatmap
dev.off()

pdf("ensemble_DA/effect_subtype_cntrl_cd.pdf", width = 6, height = 4)
subtype_cntrl_cd$effectplot + xlim(c(-1.5,3)) + scale_fill_manual("Legend",values = c("Control [fdr<0.05]"= alpha("dodgerblue",1),
                                                                                    "CD [fdr<0.05]"= alpha("firebrick",1),
                                                                                    "Control [fdr<0.2]"= alpha("dodgerblue",0.5),
                                                                                    "CD [fdr<0.2]"= alpha("firebrick",0.5)))
dev.off()




subtype_cd_uc <- da_ensemble(ps2, prev = 0.9, agglom = "Genus", variable = "Subtype", groups = c("CD", "UC"), normalise = "TSS")

pdf("ensemble_DA/heatmap_subtype_cd_uc.pdf", width = 4, height = 4)
subtype_cd_uc$heatmap
dev.off()

pdf("ensemble_DA/effect_subtype_cd_uc.pdf", width = 5, height = 3.5)
subtype_cd_uc$effectplot + xlim(c(-1.5,3)) + scale_fill_manual("Legend",values = c("CD [fdr<0.05]"= alpha("dodgerblue",1),
                                                                                      "UC [fdr<0.05]"= alpha("firebrick",1),
                                                                                      "CD [fdr<0.2]"= alpha("dodgerblue",0.5),
                                                                                      "UC [fdr<0.2]"= alpha("firebrick",0.5)))
dev.off()



neoplasia_ibd <- da_ensemble(ps2, prev = 0.9, agglom = "Genus", variable = "Neoplasia.st", groups = c("N0.IBD", "Nx.IBD"), normalise = "TSS")

pdf("ensemble_DA/heatmap_neoplasia_ibd.pdf", width = 4, height = 4)
neoplasia_ibd$heatmap
dev.off()

pdf("ensemble_DA/effect_neoplasia_ibd.pdf", width = 5, height = 3.5)
neoplasia_ibd$effectplot + xlim(c(-1.5,3)) + scale_fill_manual("Legend",values = c("N0.IBD [fdr<0.05]"= alpha("dodgerblue",1),
                                                                                   "Nx.IBD [fdr<0.05]"= alpha("firebrick",1),
                                                                                   "N0.IBD [fdr<0.2]"= alpha("dodgerblue",0.5),
                                                                                   "Nx.IBD [fdr<0.2]"= alpha("firebrick",0.5)))
dev.off()


neoplasia_uc <- da_ensemble(ps2, prev = 0.9, agglom = "Genus", variable = "Subtype_neoplasia", groups = c("UC.N0", "UC.Nx"), normalise = "TSS")

pdf("ensemble_DA/heatmap_neoplasia_uc.pdf", width = 3.0, height = 2.4)
neoplasia_uc$heatmap
dev.off()

pdf("ensemble_DA/effect_neoplasia_uc.pdf", width = 5, height = 3.5)
neoplasia_uc$effectplot + xlim(c(-1.5,3)) + scale_fill_manual("Legend",values = c("UC.N0 [fdr<0.05]"= alpha("dodgerblue",1),
                                                                                   "UC.Nx [fdr<0.05]"= alpha("firebrick",1),
                                                                                   "UC.N0 [fdr<0.2]"= alpha("dodgerblue",0.5),
                                                                                   "UC.Nx [fdr<0.2]"= alpha("firebrick",0.5)))
dev.off()


endo_act <- da_ensemble(ps2.uc, prev = 0.9, agglom = "Genus", variable = "Combined_endoscopic_activity.bin", groups = c("Low", "High"), normalise = "TSS")

pdf("ensemble_DA/heatmap_endo_act.pdf", width = 3.1, height = 2.3)
endo_act$heatmap
dev.off()

pdf("ensemble_DA/effect_endo_act.pdf", width = 5, height = 2.5)
endo_act$effectplot + xlim(c(-1.5,1)) + scale_fill_manual("Legend",values = c("Low [fdr<0.05]"= alpha("dodgerblue",1),
                                                                                   "High [fdr<0.05]"= alpha("firebrick",1),
                                                                                   "Low [fdr<0.2]"= alpha("dodgerblue",0.5),
                                                                                   "High [fdr<0.2]"= alpha("firebrick",0.5)))
dev.off()


clin_act <- da_ensemble(ps2, prev = 0.9, agglom = "Genus", variable = "Clinical_activity", groups = c("Low", "High"), normalise = "TSS")

pdf("ensemble_DA/heatmap_clin_act.pdf", width = 4, height = 4)
clin_act$heatmap
dev.off()

pdf("ensemble_DA/effect_clin_act.pdf", width = 5, height = 3.5)
clin_act$effectplot + xlim(c(-1.5,3)) + scale_fill_manual("Legend",values = c("Low [fdr<0.05]"= alpha("dodgerblue",1),
                                                                              "High [fdr<0.05]"= alpha("firebrick",1),
                                                                              "Low [fdr<0.2]"= alpha("dodgerblue",0.5),
                                                                              "High [fdr<0.2]"= alpha("firebrick",0.5)))
dev.off()


psc <- da_ensemble(ps2.ibd, prev = 0.9, agglom = "Genus", variable = "PSC_yes_no", groups = c("No", "Yes"), normalise = "TSS")

pdf("ensemble_DA/heatmap_psc.pdf", width = 4, height = 4)
psc$heatmap
dev.off()

pdf("ensemble_DA/effect_psc.pdf", width = 5, height = 3.5)
psc$effectplot + scale_fill_manual("Legend",values = c("No [fdr<0.05]"= alpha("dodgerblue",1),
                                                                              "Yes [fdr<0.05]"= alpha("firebrick",1),
                                                                              "No [fdr<0.2]"= alpha("dodgerblue",0.5),
                                                                              "Yes [fdr<0.2]"= alpha("firebrick",0.5)))
dev.off()




pdf("ensemble_DA/Agathobacter_UC.pdf", width = 2.3, height = 2.3)
boxplot_individual_taxa(ps2.uc, taxa = "Agathobacter", norm = "TSS", agglom = "Genus", variable = "Neoplasia2", 
                        comparisons = list(c("N0", "Nx"))) + 
  scale_colour_manual("Neoplasia", values = c("blue", "red4")) + 
  scale_fill_manual("Neoplasia",values = c("blue", "red4"))
dev.off()


pdf("ensemble_DA/Escherichia_Shigella_UC.pdf", width = 2.3, height = 2.3)
boxplot_individual_taxa(ps2.uc, taxa = "Escherichia-Shigella", norm = "TSS", agglom = "Genus", variable = "Neoplasia2", 
                        comparisons = list(c("N0", "Nx"))) + 
  scale_colour_manual("Neoplasia", values = c("blue", "red4")) + 
  scale_fill_manual("Neoplasia",values = c("blue", "red4"))
dev.off()

pdf("ensemble_DA/Lachnospira_UC.pdf", width = 2.3, height = 2.3)
boxplot_individual_taxa(ps2.uc, taxa = "Lachnospira", norm = "TSS", agglom = "Genus", variable = "Neoplasia2", 
                        comparisons = list(c("N0", "Nx"))) + 
  scale_colour_manual("Neoplasia", values = c("blue", "red4")) + 
  scale_fill_manual("Neoplasia",values = c("blue", "red4"))
dev.off()
