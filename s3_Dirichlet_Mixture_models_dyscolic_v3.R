

### Determining optimal number of clusters using Dirichlet Multinomial Mixtures (Holmes I, Harris K, Quince C. PlosONE, 2012: doi.org/10.1371/journal.pone.0030126) ###

library(ggrepel)
library(ggplot2)
library(DirichletMultinomial)
library(pheatmap)
library(RColorBrewer)


# 1. Agglomeration at genus level and make count OTU table


otus <- assign_taxa_names(ps2, agglom = "Genus")
otus <- as.matrix(as.data.frame(otus))


# 2. Cluster - even with set.seed there appears to be some variability in this so test for convergence



#checking for convergence of DMM
fit.list <- vector("list", 20)


for(i in 1:20){
  # 1. Agglomeration at genus level and make count OTU table
  
  
  otus <- assign_taxa_names(ps2, agglom = "Genus")
  otus <- as.matrix(as.data.frame(otus))
  
  # 2. Cluster - even with set.seed there appears to be some variability in this
  fit <- mclapply(1:7, dmn, count=otus, verbose=TRUE)
  fit.list[[i]] <- mclapply(1:7, dmn, count=otus, verbose=TRUE)
  
  # 3. Laplace
  lplc <- sapply(fit, laplace)

  
  # 4. Find best
  
  best <- fit[[which.min(lplc)]]
  
  
  
  
  if(i == 1){
    run1 <- paste("C", mixture(best, assign = TRUE), sep = ".")
  } else {
    run1 <- cbind(run1, paste("C", mixture(best, assign = TRUE), sep = "."))
  }
}

rownames(run1) <- rownames(sample_data(ps2))
colnames(run1) <- paste("run", seq(1,ncol(run1),1), sep = ".")
run2 <- run1

run2 <- data.frame(gsub("C.", "", run2), stringsAsFactors = F)
run2 <- apply(run2, 2, function(x) as.numeric(x))
dist.runs <- as.matrix(dist(t(run2)))
diag(dist.runs) <- 1
convergent_runs <- dist.runs[!(rowSums(dist.runs == 0) == 0), !(colSums(dist.runs == 0) == 0)]
non_convergent_runs <- dist.runs[rowSums(dist.runs == 0) == 0, colSums(dist.runs == 0) == 0]
ord.vec <- dist.runs[,colnames(convergent_runs)[1]]
ord.vec[names(ord.vec) == colnames(convergent_runs)[1]] <- 0

pdf("Clusters/convergence.pdf", width = 5, height = 5)
pheatmap(run2[,order(ord.vec)], cluster_cols = F)
dev.off()



#extract first converged run
index.con <- seq(1,ncol(convergent_runs),1)[colnames(convergent_runs) == colnames(convergent_runs)[1]]
fit <- fit.list[[index.con]]


# 3. Laplace
lplc <- sapply(fit, laplace)
plot(lplc, type="b", xlab="Number of Dirichlet Components",
     ylab="Model Fit")

# 4. Find best

(best <- fit[[which.min(lplc)]])

mixturewt(best)
pi     theta
1 0.3701018 38.327875
2 0.3441818 42.378170
3 0.2857165  8.699905

#confirm match before assigning
all(names(mixture(best, assign = TRUE)) == rownames(sample_data(ps2)))


sample_data(ps2)$Cluster <- paste("C", mixture(best, assign = TRUE), sep = ".")

#add
ps2.ibd <- subset_samples(ps2, sample_data(ps2)$Subtype != "Control")

# 7. Identify bacteria

p0 <- fitted(fit[[1]], scale=TRUE)     # scale by theta
p3 <- fitted(best, scale=TRUE)
colnames(p3) <- paste("c", 1:3, sep="")
(meandiff <- colSums(abs(p3 - as.vector(p0))))
sum(meandiff)

diff <- rowSums(abs(p3 - as.vector(p0)))
o <- order(diff, decreasing=TRUE)
cdiff <- cumsum(diff[o]) / sum(diff)
df <- as.data.frame(head(cbind(Mean=p0[o], p3[o,], diff=diff[o], cdiff), 40))



#Figure 4a - row-scaled heatmap of top 40 genus contributors to DMM clusters
pdf("Clusters/heatmap_40_3.pdf", width = 4.2, height = 6.2)
pheatmap(df[,c(2:4)], scale = "row", cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(brewer.pal(n=7,name = "RdBu")))(100))
dev.off()





###################################################

          #### Exploratory triplot ####

###################################################

#Figure 4b
library(seqgroup)

otus <- assign_taxa_names(ps2.ibd, "Genus")
otus <- sweep(otus, 1, rowSums(otus), "/")

ibd_metadata2 <- as.data.frame(
  as.matrix(sample_data(ps2.ibd)[,c("PSC_yes_no", "Montreal_combined", "Subtype_neoplasia", "Combined_endoscopic_activity.bin",
                                    "Combined_activity.bin", "anti_TNF", "Immunomodulators", "Oral_ASA", 
                                    "Cluster", "Age", "Duration", "BMI")]))

colnames(ibd_metadata2) <- c("PSC:", "Montreal:", "Neoplasia:", "Activity(endoscopy):", "Activity(clinical):", 
                             "anti_TNF:", "Immunomodulators:", "Oral_ASA", "C.", "Age", "Duration", "BMI")


ibd_metadata2[,"C."] <- as.character(ibd_metadata2[,"C."])
ibd_metadata2[,"C."][ibd_metadata2[,"C."] == "C.1"] <- "1"
ibd_metadata2[,"C."][ibd_metadata2[,"C."] == "C.2"] <- "2"
ibd_metadata2[,"C."][ibd_metadata2[,"C."] == "C.3"] <- "3"

ibd_metadata2 <- assignMetadataTypes(ibd_metadata2, categoric = c("PSC:", "Montreal:", "Neoplasia:", "Activity(endoscopy):", "Activity(clinical):",
                                                                  "anti_TNF:", "Immunomodulators:", "Oral_ASA", "C."))


na.indices=unique(which(is.na(ibd_metadata2),arr.ind=TRUE)[,1])
indices.to.keep=setdiff(1:nrow(ibd_metadata2),na.indices)
ibd_metadata2=ibd_metadata2[indices.to.keep,]
ibd_taxa2=t(otus)[,indices.to.keep]
set.seed(101)

seqPCoA(ibd_taxa2,groups=ibd_metadata2$'C.', metadata = ibd_metadata2, topTaxa=10, 
        xlim = c(-0.31,0.37), ylim = c(-0.25,0.35), drawEllipse = TRUE, centroidFactor = 0.5, arrowFactor = 0.4, colors = unlist(lapply(ibd_metadata2$'C.', function(x) if(x == "1"){"darkseagreen1"} else if(x == "2") {"purple"}else{"indianred1"})), taxonColor = "blue",
        metadataColor = "grey12", qvalThreshold = 0.2)


###################################################

        #### Correspondence analysis ####

###################################################

#Figure 4c
library("FactoMineR")
library("factoextra")
if(!dir.exists("CA")) dir.create("CA")


#make a data frame from clinical data
mod.df <- as.data.frame(sample_data(ps2.ibd))
mod.df$Montreal_combined <- recode(mod.df$Montreal_combined, E1="MontrealUC:E1/2", E2="MontrealUC:E1/2", E3="MontrealUC:E3", L1="MontrealCD:L1/L3", L3="MontrealCD:L1/L3", L2="MontrealCD:L2")
mod.df$Neoplasia2 <- recode(mod.df$Neoplasia2, N0="No Neoplasia", Nx="Neoplasia")
mod.df$Subtype_neoplasia <- paste(mod.df$Subtype, mod.df$Neoplasia2, sep = ":")
mod.df$PSC_yes_no <- recode(mod.df$PSC_yes_no, Yes="PSC:Yes", No="PSC:No")
mod.df$Combined_endoscopic_activity.bin <- recode(mod.df$Combined_endoscopic_activity.bin, "Low"="Endoscopically active:No", "High"="Endoscopically active:Yes")
mod.df$Clinical_activity <- recode(mod.df$Clinical_activity, High="Clinically active:Yes", Low="Clinically active:No")


cont.tab <- rbind(table(mod.df$Subtype_neoplasia, mod.df$Cluster), 
                  table(mod.df$PSC_yes_no, mod.df$Cluster),
                  table(mod.df$Combined_endoscopic_activity.bin, mod.df$Cluster),
                  table(mod.df$Clinical_activity, mod.df$Cluster),
                  table(mod.df$Montreal_combined, mod.df$Cluster))


res.ca <- CA(cont.tab, graph = TRUE)

pdf("CA/Biplot_arrow_ibd.pdf", width = 5.7, height = 4.5)
fviz_ca_biplot(res.ca, map ="colgreen", arrow = c(TRUE, FALSE),
               repel = TRUE, title = "CA Biplot - IBD Samples")
dev.off()

#model
summary(res.ca)

#The chi square of independence between the two variables is equal to 57.57947 (p-value =  0.000351622 ).








#############################################################

                 ### Relative risk ###

#############################################################

#Figure 4d

ps2.ibd <- subset_samples(ps2, sample_data(ps2)$Subtype != "Control")

# Cross tabulate
ps2.cd <- subset_samples(ps2.ibd, sample_data(ps2.ibd)$Subtype == "CD")
ps2.uc <- subset_samples(ps2.ibd, sample_data(ps2.ibd)$Subtype == "UC")

counts.cd <- table(sample_data(ps2.cd)$Neoplasia, sample_data(ps2.cd)$Cluster)
counts.uc <- table(sample_data(ps2.uc)$Neoplasia, sample_data(ps2.uc)$Cluster)

counts.cd <- rbind(counts.cd,Nx=colSums(counts.cd[2:5,]))
counts.uc <- rbind(counts.uc,Nx=colSums(counts.uc[2:5,]))

cd.prop <- sweep(counts.cd, 2, colSums(counts.cd[1:5,]), "/")
uc.prop <- sweep(counts.uc, 2, colSums(counts.uc[1:5,]), "/")



####### Figure 4d #######
# Make a bar graph
library(tidyr)
pal <- c("blue", "purple", "pink2", "salmon", "red", "red4")

cd.prop.long <- gather(as.data.frame(cbind(cd.prop, Neoplasia = rownames(cd.prop))), key = "Cluster", value = "Percentage", C.1:C.3)
uc.prop.long <- gather(as.data.frame(cbind(uc.prop, Neoplasia = rownames(uc.prop))), key = "Cluster", value = "Percentage", C.1:C.3)
cd.prop.long$Percentage <- as.numeric(as.character(cd.prop.long$Percentage))
uc.prop.long$Percentage <- as.numeric(as.character(uc.prop.long$Percentage))

cd.prop.long$Percentage <- cd.prop.long$Percentage*100
uc.prop.long$Percentage <- uc.prop.long$Percentage*100
x <- rbind(cbind(cd.prop.long, Condition = rep("CD", nrow(cd.prop.long))), cbind(uc.prop.long, Condition = rep("UC", nrow(uc.prop.long))))
x$Neoplasia <- factor(x$Neoplasia, levels = c("N0", "A1", "N1", "N2", "N3", "Nx"))

pdf("Clusters/rel_risk_barplots.pdf", width = 6.5, height = 4)
ggplot(x, aes(x=Neoplasia, y=Percentage, fill=Neoplasia)) + geom_bar(stat = "identity") + 
  facet_wrap(Condition~Cluster, scales = "free") + theme_classic() + scale_fill_manual(values = pal) +
  geom_text(aes(label=round(Percentage, 1)), vjust=-0.3, size=2.5) + ylim(0,100)
dev.off()

mod.df <- data.frame(sample_data(ps2.ibd))
calc_relative_risk(x=mod.df[mod.df$Subtype == "UC",], variable1 = "Cluster", variable2 = "Neoplasia2", contr1 = c("C.3", "C.1"), contr2 = c("Nx","N0"))
[1]  4.073333333  1.616591304 10.263598726  0.002894705
calc_relative_risk(x=mod.df[mod.df$Subtype == "UC",], variable1 = "Cluster", variable2 = "Neoplasia2", contr1 = c("C.2", "C.1"), contr2 = c("Nx","N0"))
[1] 3.565517241 1.378314289 9.223522746 0.008749873

#proportion of UC with PSC in each cluster
sweep(table(sample_data(ps2.uc)$PSC_yes_no, sample_data(ps2.uc)$Cluster),2,colSums(table(sample_data(ps2.uc)$PSC_yes_no, sample_data(ps2.uc)$Cluster)),"/")


















