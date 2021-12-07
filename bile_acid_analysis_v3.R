#Bile acids

#import bile acid data
bile <- read.table("Bile_acids/Bile_acids_import.txt", header = T)

#subset to include only subjects with microbiome and bile acid data
ps2.bile <- subset_samples(ps2, rownames(sample_data(ps2)) %in% rownames(bile))
all(rownames(bile) == rownames(sample_data(ps2.bile)))

#subset to IBD
ps2.bile.ibd <- subset_samples(ps2.bile, sample_data(ps2.bile)$Subtype != "Control")
bile.ibd <- bile[rownames(bile) %in% rownames(sample_data(ps2.bile.ibd)),]






###########################################################

                      ### PCA ###

###########################################################


#replace min
addMin <- function(df){
  cm <- apply(df, 2, function (x) min(x[x>0])/5)
  df2 <- df
  for(i in 1:length(cm)){
    df2[,i][df2[,i] == 0] <- cm[i]
  }
  return(df2)
}
bile.ibd.comp <- bile.ibd[,1:28][,colSums(bile.ibd[,1:28]) > 0]
bile.ibd.comp <- addMin(bile.ibd.comp)

#log
bile.ibd.comp <- log2(bile.ibd.comp)

#center
bile.ibd.comp <- sweep(bile.ibd.comp,2,colMeans(bile.ibd.comp),"-")

library(ade4)
library(factoextra)
#all
sampleTable <- sample_data(ps2.bile.ibd)

gene.pac <- dudi.pca(bile.ibd.comp, nf=3, scale = FALSE, center = FALSE, scannf = FALSE)
scrs <- gene.pac$li

if(!(is.null(rownames(sample_data(ps2.bile.ibd)))) & !(is.null(rownames(scrs)))){
  if(all(rownames(sample_data(ps2.bile.ibd)) == rownames(scrs))){
    scrs <- data.frame(scrs, sample_data(ps2.bile.ibd), bile.ibd[,29:ncol(bile.ibd)])
  }
}



#Clusters PCA

######## Figure 5a ########
#Neoplasia2 PCA

pal <- c("blue", "red4")

pca.df <- scrs

eval.1 <- round(gene.pac$eig[1]/sum(gene.pac$eig),3)*100
eval.2 <- round(gene.pac$eig[2]/sum(gene.pac$eig),3)*100


p2 <- ggplot(pca.df, aes(x=Axis1, y=Axis2, colour=Neoplasia2, fill=Neoplasia2)) + geom_point(size=1.2, alpha=1, shape=21) + scale_fill_manual("Neoplasia2", values = pal) + 
  scale_colour_manual("Neoplasia", values = pal) + scale_fill_manual("Neoplasia", values = pal) + stat_ellipse(level = 0.8) + theme_classic() + 
  xlab(paste("PCA 1 [", round(eval.1, 2), "%]", sep = "")) + ylab(paste("PCA 2 [", round(eval.2, 2), "%]", sep = "")) + 
  theme(legend.position = "right",
        axis.line = element_line(colour = 'black', size = 0.56), 
        axis.ticks = element_line(colour = 'black', size = 0.56))

pdf("Bile_acids/PCA_neoplasia2.pdf", width = 3.3, height = 1.8)
p2
dev.off()

adonis(bile.ibd.comp~sampleTable$Neoplasia2, method = "euclidean")

######## Figure 5b ########
pal <- c("darkseagreen1","purple","indianred1")

pca.df <- scrs

eval.1 <- round(gene.pac$eig[1]/sum(gene.pac$eig),3)*100
eval.2 <- round(gene.pac$eig[2]/sum(gene.pac$eig),3)*100

p2 <- ggplot(pca.df, aes(x=Axis1, y=Axis2, colour=Cluster, fill=Cluster)) + geom_point(size=1.2, alpha=1, shape=21) + scale_fill_manual("Cluster", values = pal) + 
  scale_colour_manual("Cluster", values = pal) + scale_fill_manual("Cluster", values = pal) + stat_ellipse(level = 0.8) + theme_classic() + 
  xlab(paste("PCA 1 [", round(eval.1, 2), "%]", sep = "")) + ylab(paste("PCA 2 [", round(eval.2, 2), "%]", sep = "")) + 
  theme(legend.position = "right",
        axis.line = element_line(colour = 'black', size = 0.56), 
        axis.ticks = element_line(colour = 'black', size = 0.56))

pdf("Bile_acids/PCA_cluster.pdf", width = 3.3, height = 1.8)
p2
dev.off()

adonis(bile.ibd.comp~sampleTable$Cluster, method = "euclidean")



######## Figure S6a ########
#Neoplasia PCA

pal <- c("blue", "purple", "peachpuff", "salmon", "red1")

pca.df <- scrs

eval.1 <- round(gene.pac$eig[1]/sum(gene.pac$eig),3)*100
eval.2 <- round(gene.pac$eig[2]/sum(gene.pac$eig),3)*100

p2 <- ggplot(pca.df, aes(x=Axis1, y=Axis2, colour=Neoplasia, fill=Neoplasia)) + geom_point(size=1.2, alpha=1, shape=21) + scale_fill_manual("Neoplasia", values = pal) + 
  scale_colour_manual("Neoplasia", values = pal) + scale_fill_manual("Neoplasia", values = pal) + stat_ellipse(level = 0.8) + theme_classic() + 
  xlab(paste("PCA 1 [", round(eval.1, 2), "%]", sep = "")) + ylab(paste("PCA 2 [", round(eval.2, 2), "%]", sep = "")) + 
  theme(legend.position = "right",
        axis.line = element_line(colour = 'black', size = 0.56), 
        axis.ticks = element_line(colour = 'black', size = 0.56))

pdf("Bile_acids/PCA_neoplasia.pdf", width = 3.3, height = 1.8)
p2
dev.off()


adonis(bile.ibd.comp[sampleTable$Neoplasia %in% c("N0", "N3"),]~sampleTable$Neoplasia[sampleTable$Neoplasia %in% c("N0", "N3")], method = "euclidean")






####### Figure S6b #######
#Primary:Secondary BAs PCA

#add min for zero value in primary bile acids and re-calculate ratio column
pca.df[,c("Primary_BA", "Secondary_BA")] <- addMin(pca.df[,c("Primary_BA", "Secondary_BA")])
pca.df$Primary_BA.Secondary_BA <- pca.df$Primary_BA/pca.df$Secondary_BA
p.pca <- ggplot(pca.df, aes(x=Axis1, y=Axis2, colour=log10(Primary_BA.Secondary_BA), group=Cluster)) + 
  geom_point(size=1.2, alpha=1) + stat_ellipse(level = 0.8, color="grey25", lty=2) + theme_classic() + 
  xlab(paste("PCA 1 [", round(eval.1, 2), "%]", sep = "")) + 
  ylab(paste("PCA 2 [", round(eval.2, 2), "%]", sep = "")) + 
  theme(legend.position = "right",
        axis.line = element_line(colour = 'black', size = 0.56), 
        axis.ticks = element_line(colour = 'black', size = 0.56))

p3 <- p.pca + scale_fill_gradient2(parse(text=paste("1^o", "2^o", "l", sep = ":")),
                                    low = "blue",
                                    mid = "grey92",
                                    high = "red",
                                    midpoint = 0,
                                    space = "Lab",
                                    na.value = "green",
                                    guide = "colourbar",
                                    aesthetics = "colour"
)

pdf("Bile_acids/PCA_prim_sec_gradient.pdf", width = 3.3, height = 1.8)
p3
dev.off()

####### Figure S6c #######
#Proportion unconjugated

p.pca <- ggplot(pca.df, aes(x=Axis1, y=Axis2, colour=1-BA_nonconj.AB_Tot, group=Cluster)) + 
  geom_point(size=1.2, alpha=1) + stat_ellipse(level = 0.8, color="grey25", lty=2) + theme_classic() + 
  xlab(paste("PCA 1 [", round(eval.1, 2), "%]", sep = "")) + 
  ylab(paste("PCA 2 [", round(eval.2, 2), "%]", sep = "")) + 
  theme(legend.position = "right",
        axis.line = element_line(colour = 'black', size = 0.56), 
        axis.ticks = element_line(colour = 'black', size = 0.56))

p3 <- p.pca + scale_fill_gradient("conj",
                                   low = "grey",
                                   
                                   high = "red2",
                                   #midpoint = 0.971,
                                   space = "Lab",
                                   na.value = "green",
                                   guide = "colourbar",
                                   aesthetics = "colour"
)

pdf("Bile_acids/PCA_unconj.pdf", width = 3.3, height = 1.8)
p3
dev.off()


####### Figure S6d #######
#Proportion sulphated

p.pca <- ggplot(pca.df, aes(x=Axis1, y=Axis2, colour=Sulfo.BA_Total, group=Cluster)) + 
  geom_point(size=1.2, alpha=1) + stat_ellipse(level = 0.8, color="grey25", lty=2) + theme_classic() + 
  xlab(paste("PCA 1 [", round(eval.1, 2), "%]", sep = "")) + 
  ylab(paste("PCA 2 [", round(eval.2, 2), "%]", sep = "")) + 
  theme(legend.position = "right",
        axis.line = element_line(colour = 'black', size = 0.56), 
        axis.ticks = element_line(colour = 'black', size = 0.56))

p3 <- p.pca + scale_fill_gradient("sulf",
                                   low = "grey",
                                   
                                   high = "orange",
                                   #midpoint = 0.971,
                                   space = "Lab",
                                   na.value = "green",
                                   guide = "colourbar",
                                   aesthetics = "colour"
)

pdf("Bile_acids/PCA_sulfated.pdf", width = 3.3, height = 1.8)
p3
dev.off()



###########################################################

              ### Differential abundance ###

###########################################################
######## Not plotted ########

#individual bile acid metabolites - differential abundance between different neoplasia groups in IBD

#IBD N0 vs N3
x.all <- aldex2_wrap_metab(ps2.bile, metabolites=bile[,1:28], prev = 0.9, variable = "Neoplasia3", groups = c("IBD.N0", "IBD.N3"))

x.sig <- x.all[x.all$wi.eBH < 0.2,]
x.sig <- x.sig[order(x.sig$diff.btw),]
x.sig$taxa <- rownames(x.sig)
x.sig$diff <- unlist(lapply(x.sig$diff.btw, function(y) if(y > 0){"Nx"} else {"N0"}))
x.sig$FDR <- unlist(lapply(x.sig$wi.eBH, function(y) if(y < 0.05){"[fdr<0.05]"} else {"[fdr<0.2]"}))

#no DA metabolites
x.sig

#IBD N0 vs Nx
x.all <- aldex2_wrap_metab(ps2.bile, metabolites=bile[,1:28], prev = 0.9, variable = "Neoplasia.st", groups = c("N0.IBD", "Nx.IBD"))

x.sig <- x.all[x.all$wi.eBH < 0.2,]
x.sig <- x.sig[order(x.sig$diff.btw),]
x.sig$taxa <- rownames(x.sig)
x.sig$diff <- unlist(lapply(x.sig$diff.btw, function(y) if(y > 0){"Nx"} else {"N0"}))
x.sig$FDR <- unlist(lapply(x.sig$wi.eBH, function(y) if(y < 0.05){"[fdr<0.05]"} else {"[fdr<0.2]"}))

#no DA metabolites
x.sig

#UC N0 vs Nx
x.all <- aldex2_wrap_metab(ps2.bile, metabolites=bile[,1:28], prev = 0.9, variable = "Subtype_neoplasia", groups = c("UC.N0", "UC.Nx"))

x.sig <- x.all[x.all$wi.eBH < 0.2,]
x.sig <- x.sig[order(x.sig$diff.btw),]
x.sig$taxa <- rownames(x.sig)
x.sig$diff <- unlist(lapply(x.sig$diff.btw, function(y) if(y > 0){"Nx"} else {"N0"}))
x.sig$FDR <- unlist(lapply(x.sig$wi.eBH, function(y) if(y < 0.05){"[fdr<0.05]"} else {"[fdr<0.2]"}))

#no DA metabolites
x.sig

#CD N0 vs Nx
x.all <- aldex2_wrap_metab(ps2.bile, metabolites=bile[,1:28], prev = 0.9, variable = "Subtype_neoplasia", groups = c("CD.N0", "CD.Nx"))

x.sig <- x.all[x.all$wi.eBH < 0.2,]
x.sig <- x.sig[order(x.sig$diff.btw),]
x.sig$taxa <- rownames(x.sig)
x.sig$diff <- unlist(lapply(x.sig$diff.btw, function(y) if(y > 0){"Nx"} else {"N0"}))
x.sig$FDR <- unlist(lapply(x.sig$wi.eBH, function(y) if(y < 0.05){"[fdr<0.05]"} else {"[fdr<0.2]"}))

#no DA metabolites
x.sig



###########################################################

                    ### Barplots ###

###########################################################
######## Figure 5c ########

library(tidyr)
#Bile acid composition
bile.ibd.prop <- bile.ibd[,1:27]
bile.ibd.prop <- sweep(bile.ibd.prop, 1, rowSums(bile.ibd.prop), "/")

if(!(is.null(rownames(bile.ibd))) & !(is.null(rownames(sample_data(ps2.bile.ibd))))){
  if(all(rownames(bile.ibd) == rownames(sample_data(ps2.bile.ibd)))){
    bile.ibd.df <- cbind(SampleID = rownames(sample_data(ps2.bile.ibd)),sample_data(ps2.bile.ibd)[,c("Neoplasia", "Cluster")], bile.ibd.prop)
  }
}

bile.ibd.df.long <-  gather(bile.ibd.df, key = "Bile_acid", value = "Proportion", TUDCA:HCA)
bile.ibd.df.long$Bile_acid <- factor(bile.ibd.df.long$Bile_acid, 
                                     levels = c("CA", "CA.3S", "GCA", "TCA",
                                                "CDCA", "CDCA.3S", "GCDCA", "TCDCA",
                                                "HCA",
                                                "DCA", "DCA.3S", "GDCA", "TDCA",
                                                "LCA", "LCA.3S", "GLCA", "GLCA.3S", "TLCA", "TLCA.3S",
                                                "UDCA", "UDCA.3S", "GUDCA", "GUDCA.3S", "TUDCA", "TUDCA.3S",
                                                "HDCA", "THDCA"))
bile.ibd.df.long$SampleID <- factor(bile.ibd.df.long$SampleID, 
                                    levels = rev(rownames(bile.ibd.df)[order(bile.ibd.df$CA)]))
pal <- c("#FF3300", "#FF9966", "#CC3300", "#CC6666", "#990033", "#CC6600", "#FFCC33", "#FFFF66", "#CC9900",
         "#0000FF", "#99CCFF", "#003366", "#66FFFF", "#0099FF", "#33FF66", "#3399FF", "#006633", "#009966", 
         "#CCFF33", "#003300", "#9933FF", "#6600CC", "#CC33FF", "#9900CC", "#FF00FF", "#996699", "#9999CC")
p.bar <- ggplot(bile.ibd.df.long, aes(x=SampleID, y=Proportion, fill=Bile_acid)) + geom_bar(stat = "identity") + 
  facet_grid(~Cluster, scales = "free") + scale_fill_manual("Bile acids",values = pal) + xlab("") + 
  theme_classic() + 
  theme(text = element_text(size=20), axis.text.x=element_blank(), axis.text = element_text(size=16))

#plot
p.bar


###########################################################

            ### Bile acid boxplots ###

###########################################################
######## Figure 5d & 5e ########

bile.b <- bile[,c("CA.1", "CDCA.1", "DCA.1", "LCA.1", "UDCA.1", "BA_nonconj.AB_Tot", "Glyco.BA_Total", "Tauro.BA_Total", "Sulfo.BA_Total")]
bile.b[,c("CA.1", "CDCA.1", "DCA.1", "LCA.1", "UDCA.1")] <- sweep(bile.b[,c("CA.1", "CDCA.1", "DCA.1", "LCA.1", "UDCA.1")],1,bile$Total_BA,"/")
bile.m <- reshape2::melt(t(as.data.frame(bile.b)))
colnames(bile.m)[1] <- "BA"
colnames(bile.m)[2] <- "Sample"
colnames(bile.m)[3] <- "Abundance"
bile.m$BA <- as.character(bile.m$BA)
bile.m$Sample <- as.character(bile.m$Sample)
sampdf <- data.frame(sample_data(ps2.bile), stringsAsFactors = FALSE)
sampdf$Sample <- rownames(sampdf)
bile.merge <- merge(bile.m, sampdf, by.x = "Sample")

ggplot(bile.merge, aes(x=Cluster, y=Abundance, colour=Cluster)) + geom_boxplot(outlier.colour = NA, colour="black") + theme_classic() +
  geom_jitter(alpha=0.5) + facet_wrap(~BA, scales = "free_y", ncol = 5) + 
  stat_compare_means(comparisons = list(c("C.2", "C.3"),c("C.1", "C.2"),c("C.1", "C.3")), label = "p.signif", tip.length = 0) + 
  scale_colour_manual(values = c("darkseagreen1","purple","indianred1")) + ylim(c(-0.1,1.3))


###########################################################

           ### Genus bile acid correlation ###

###########################################################
######## Figure 5f ########

library(compositions)
library(zCompositions)

otus.b <- assign_taxa_names(ps2.bile, "Genus")
otus.b <- otus.b[,colSums(otus.b == 0)/nrow(otus.b) <= 0.9]
otus.b <- cmultRepl(otus.b, method = "CZM")
otus.b <- clr(otus.b)
write.table(t(otus.b), "otus.b.txt", sep = "\t")

bile.b <- round(bile[,c("CA.1", "CDCA.1", "DCA.1", "LCA.1", "UDCA.1")])
bile.b <- bile.b[,colSums(bile.b == 0)/nrow(bile.b) <= 0.9]
bile.b <- cmultRepl(bile.b, method = "CZM")
bile.b <- clr(bile.b)
write.table(t(bile.b), "bile.b.txt", sep = "\t")

library(Hmisc)
otus.x <- as.matrix(otus.b)
meta_cont <- as.matrix(bile.b)
prac.corr <- rcorr(meta_cont, otus.x, type = "spearman")
prac.corrR <- prac.corr$r[1:5,6:128]
prac.corrP <- prac.corr$p[1:5,6:128]
prac.corrR[is.na(prac.corrP)] <- 0
prac.corrR[p.adjust(prac.corrP, "fdr") >= 0.05] <- 0
prac.corrR[abs(prac.corrR) < 0.4] <- 0
prac.corrR <- prac.corrR[!(rowSums(prac.corrR == 0) == ncol(prac.corrR)),!(colSums(prac.corrR== 0) == nrow(prac.corrR))]
library(pheatmap)
library(RColorBrewer)
df <- as.data.frame(cbind(Mean=p0[o], p3[o,], diff=diff[o], cdiff))

df <- data.frame(df, cluster=apply(df[,c(2:4)], 1, function(x) c("C.1", "C.2", "C.3")[x == max(x)]))

df.ss <- df[rownames(df) %in% colnames(prac.corrR),]
df.ss <- df.ss[match(colnames(prac.corrR), rownames(df.ss)),]
annotation.row <- data.frame(Cluster=df.ss[,7])
rownames(annotation.row) <- rownames(df.ss)
ann_colors = list(
  Cluster = c(C.2="purple", C.3="indianred1")
)

pdf("Bile_acids/heatmap_cor_bile_bact.pdf", width = 6, height = 8)
pheatmap(t(prac.corrR), color = colorRampPalette(rev(brewer.pal(n=7,name = "RdBu")))(100)[7:100], annotation_row = annotation.row, annotation_colors = ann_colors)
dev.off()


