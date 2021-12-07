
#################################################################

                 ### Functions ###

#################################################################

calc_relative_risk <- function(x, variable1=NULL, variable2=NULL, contr1=NULL,contr2=NULL){
  log.vec <- is.na(x[,variable1]) | is.na(x[,variable2])
  x <- x[!log.vec,]
  tab_rr <- table(x[,variable1], x[,variable2])
  exposed=contr1[1]
  notExposed=contr1[2]
  event=contr2[1]
  noEvent=contr2[2]
  a=tab_rr[exposed,event]
  b=tab_rr[exposed,noEvent]
  c=tab_rr[notExposed,event]
  d=tab_rr[notExposed,noEvent]
  risk_exposed <- a/sum(a+b)
  risk_notExposed <- c/sum(c+d)
  
  rel_risk <- risk_exposed/risk_notExposed
  ci_low <- exp(log(rel_risk) - 1.96*sqrt((1/a) + (1/c) - (1/(a+b)) - (1/(c+d))))
  ci_high <- exp(log(rel_risk) + 1.96*sqrt((1/a) + (1/c) - (1/(a+b)) - (1/(c+d))))
  q.val <- log(rel_risk)/sqrt((1/a) + (1/c) - (1/(a+b)) - (1/(c+d)))
  
  p.val <- (1-pnorm(abs(q.val)))*2
  return(c(rel_risk, ci_low, ci_high, p.val))
}

clean_taxa_names <- function(ps.x=NULL){
  #print warning
  warning("This function starts naming at 'Genus' level and assumes that the 6th column of 'tax_table(ps.x)' is 'Genus': in other words, the order is:\n
          \t'Kingdom' 'Phylum' 'Class' 'Order' 'Family' 'Genus' 'Species'\n\n\n", immediate. = TRUE)
  #helper functions
  log.fun <- function(x) {is.na(x) | grepl("uncultur|ambiguous", ignore.case = TRUE, x)}
  log.fun.ucg <- function(x) {grepl("^UCG", ignore.case = TRUE, x)}
  
  otu.vec <- c()
  otu.vec <- tax_table(ps.x)[,6]
  if(sum(unlist(lapply(otu.vec, log.fun))) > 0){
    temp.vec <- unlist(lapply(otu.vec, log.fun)) & !unlist(lapply(tax_table(ps.x)[,5], log.fun))
    otu.vec[temp.vec] <- paste(tax_table(ps.x)[,5][temp.vec],"g.NA", sep = "_")
  }
  if(sum(unlist(lapply(otu.vec, log.fun))) > 0){
    temp.vec <- unlist(lapply(otu.vec, log.fun)) & !unlist(lapply(tax_table(ps.x)[,4], log.fun))
    otu.vec[temp.vec] <- paste(tax_table(ps.x)[,4][temp.vec],"g.NA", sep = "_")
  }
  if(sum(unlist(lapply(otu.vec, log.fun))) > 0){
    temp.vec <- unlist(lapply(otu.vec, log.fun)) & !unlist(lapply(tax_table(ps.x)[,3], log.fun))
    otu.vec[temp.vec] <- paste(tax_table(ps.x)[,3][temp.vec],"g.NA", sep = "_")
  }
  
  if(sum(unlist(lapply(otu.vec, log.fun))) > 0){
    temp.vec <- unlist(lapply(otu.vec, log.fun)) & !unlist(lapply(tax_table(ps.x)[,2], log.fun))
    otu.vec[temp.vec] <- paste(tax_table(ps.x)[,2][temp.vec],"g.NA", sep = "_")
  }
  if(sum(unlist(lapply(otu.vec, log.fun))) > 0){
    temp.vec <- unlist(lapply(otu.vec, log.fun)) & !unlist(lapply(tax_table(ps.x)[,1], log.fun))
    otu.vec[temp.vec] <- paste(tax_table(ps.x)[,1][temp.vec],"g.NA", sep = "_")
  }
  if(sum(unlist(lapply(otu.vec, log.fun.ucg))) > 0){
    temp.vec <- unlist(lapply(otu.vec, log.fun.ucg))
    otu.vec[temp.vec] <- paste(tax_table(ps.x)[,5][temp.vec],tax_table(ps.x)[,6][temp.vec], sep = "_")
  }
  otu.vec <- as.vector(otu.vec)
  return(otu.vec)
}

assign_taxa_names <- function(ps.x=NULL, agglom=NULL){
  #ensure phyloseq object
  if(is.null(ps.x) | class(ps.x) != "phyloseq"){
    stop("'ps.x' must supply a phyloseq object")
  }
  #check that ps in correct orientation
  if(taxa_are_rows(ps.x)){
    stop("Taxa are rows - please re-orientate phyloseq object.")
  }
  if(!is.null(agglom)){
    if(!(agglom %in% c("Genus", "genus"))){
      stop("Only agglomeration at 'Genus' level currently supported by this function")
    }
  }
  #agglomerate
  if(!(is.null(agglom))){
    ps.x <- tax_glom(ps.x, taxrank = agglom, NArm = FALSE)
    otus <- otu_table(ps.x)
    colnames(otus) <- clean_taxa_names(ps.x)
  } else {
    otus <- otu_table(ps.x)
    colnames(otus) <- clean_taxa_names(ps.x)
    colnames(otus) <- paste(colnames(otus), "ASV", seq(1,length(colnames(otus)),1), sep = ".")
  }
  return(otus)
}



aldex2_wrap_metab <- function(ps.x=NULL, metabolites=NULL, prev=NULL, variable=NULL, groups=NULL, mc_samples=1000, rnd.seed=101){
  library(ALDEx2)
  #ensure phyloseq object
  if(is.null(ps.x) | class(ps.x) != "phyloseq"){
    stop("'ps.x' must supply a phyloseq object")
  }
  #check that ps in correct orientation
  if(taxa_are_rows(ps.x)){
    stop("Taxa are rows - please re-orientate phyloseq object.")
  }
  #check that rownames of 'metabolites' and 'sample_data(ps.x)' match
  if(!all(rownames(sample_data(ps.x)) == rownames(metabolites))){
    stop("The names of the samples in the metadata table from 'ps.x' and the metabolites from 'y' do not match.")
  }
  
  if(!is.null(prev)){
    metabolites <- metabolites[,colSums(metabolites == 0)/nrow(metabolites) <= prev]
  }
  metabolites <- round(metabolites)
  metabolites <- metabolites[rowSums(metabolites) > 0,]
  metabolites <- t(metabolites)
  meta_data <- data.frame(sample_data(ps.x))
  metabolites <- metabolites[,meta_data[,variable] %in% groups]
  meta_data_s <- meta_data[meta_data[,variable] %in% groups,]
  if(!all(rownames(meta_data_s) == colnames(metabolites))){
    stop("Error - sample names do not match - some samples may have been removed after prevalence filtering?")
  }
  meta_data_s[,variable] <- factor(as.character(meta_data_s[,variable]), levels = groups)
  set.seed(rnd.seed)
  x <- aldex.clr(metabolites, meta_data_s[,variable], mc.samples=mc_samples, denom="all", verbose=F)
  x.tt <- aldex.ttest(x, paired.test=FALSE, verbose=FALSE)
  x.effect <- aldex.effect(x, CI=T, verbose=FALSE)
  x.combined <- data.frame(x.tt,x.effect)
  return(x.combined)
}


boxplot_individual_taxa <- function(ps.x, taxa=NULL, norm=NULL, agglom=NULL, variable=NULL, comparisons=NULL){
  library(zCompositions)
  library(compositions)
  library(ggplot2)
  library(ggpubr)
  otus.df <- assign_taxa_names(ps.x = ps.x, agglom = agglom)
  if(norm == "TSS"){
    otus.df <- sweep(otus.df, 1, rowSums(otus.df), "/")
    plot_df <- data.frame(sample_data(ps.x), temp = otus.df[,taxa]*100)
    y.lab <- paste0(taxa, " (%)")
  } else if(norm == "CLR"){
    otus.df <- cmultRepl(otus.df, method = "CZM")
    otus.df <- clr(otus.df)
    plot_df <- data.frame(sample_data(ps.x), temp = otus.df[,taxa])
    y.lab <- paste0(taxa, " (CLR)")
  } else {
    stop("Error - 'norm' must be one of 'TSS' or 'CLR'")
  }
  taxa <- gsub("-|[[]|[]]|[ ]", ".", taxa)
  names(plot_df)[names(plot_df) == 'temp'] <- taxa
  
  g.1 <- ggplot(plot_df, aes_string(x=variable,y=taxa, colour=variable, fill=variable)) + geom_boxplot(outlier.colour = NA, alpha=0.4) + 
    geom_jitter(aes(alpha=0.7)) + theme_classic() + 
    stat_compare_means(comparisons = comparisons, label = "p.signif") + xlab("") + ylab(y.lab) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())
  return(g.1)
}




#Step 1 - load phyloseq and set working directory
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(dplyr)
#import phyloseq object
ps <- readRDS("dada2_230_210/ps_5_5_228_228.rds")


#################################################################

                  ### Setting factors ###

#################################################################

#set factors
sample_data(ps)$Subtype <- factor(sample_data(ps)$Subtype, levels = c("Control", "CD", "UC"))

sample_data(ps)$Neoplasia <- factor(sample_data(ps)$Neoplasia, levels = c("N0", "A1", "N1", "N2", "N3"))

sample_data(ps)$Neoplasia.st <- factor(sample_data(ps)$Neoplasia.st, 
                                       levels = c("N0.Control", "Nx.Control",
                                                  "N0.IBD", "Nx.IBD"))

sample_data(ps)$Subtype_neoplasia <- factor(sample_data(ps)$Subtype_neoplasia, levels=c("Control.N0", "Control.Nx", "CD.N0", "CD.Nx", "UC.N0", "UC.Nx"))

sample_data(ps)$Clinical_activity <- factor(sample_data(ps)$Clinical_activity, levels = c("Low", "High"))

sample_data(ps)$Combined_endoscopic_activity.bin <- factor(sample_data(ps)$Combined_endoscopic_activity.bin, levels = c("Low", "High"))

########################################################################

                  ###### Filtering #######

########################################################################


# Create table, number of features for each phyla and remove NA phyla
table(tax_table(ps)[, "Phylum"], exclude = NULL)
ps_phy_noNA <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))


# Define prevalence of each taxa
# (in how many samples did each taxa appear at least once)
prev0 = apply(X = otu_table(ps_phy_noNA),
              MARGIN = ifelse(taxa_are_rows(ps_phy_noNA), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})

# Define prevalence threshold as 1 - will remove ASVs present only in 1 sample
prevalenceThreshold = 1
# Execute prevalence filter, using `prune_taxa()` function
ps2 = prune_taxa((prev0 > prevalenceThreshold), ps_phy_noNA)



# check retained seqs post filtering
sum(colSums(otu_table(ps2)))/sum(colSums(otu_table(ps)))

# check range of sequences per sample
min(rowSums(otu_table(ps2))[order(rowSums(otu_table(ps2)))])
max(rowSums(otu_table(ps2))[order(rowSums(otu_table(ps2)))])






######### subset

ps2.ibd <- subset_samples(ps2, sample_data(ps2)$Subtype != "Control")


#################################################################

                    ### Alpha diversity ###

#################################################################



######alpha diversity#####
library(RColorBrewer)
library(ggpubr)
dir.create(paste(getwd(), "alpha_diversity", sep = "/"))


####### Alpha diversity comparisons between different cohorts (control, CD, UC) #######

####### Results #########
#use unfiltered data
diversity_df <- data.frame(sample_data(ps_phy_noNA), estimate_richness(ps_phy_noNA))

#overall Kruskal-Wallis
kruskal.test(diversity_df$Shannon~diversity_df$Subtype)
Kruskal-Wallis chi-squared = 8.7368, df = 2, p-value = 0.01267

#get medians and IQRs
aggregate(diversity_df$Shannon, by=list(diversity_df$Subtype), median)

Group.1        x
1 Control 3.742534
2      CD 3.442696
3      UC 3.412482

aggregate(diversity_df$Shannon, by=list(diversity_df$Subtype), IQR)

Group.1         x
1 Control 0.6316464
2      CD 0.8639390
3      UC 0.6112766

#check pairwise Wilcoxon test
wilcox.test(diversity_df$Shannon[diversity_df$Subtype != "UC"]~diversity_df$Subtype[diversity_df$Subtype != "UC"])
W = 3487, p-value = 0.009631
wilcox.test(diversity_df$Shannon[diversity_df$Subtype != "CD"]~diversity_df$Subtype[diversity_df$Subtype != "CD"])
W = 3455, p-value = 0.004851
wilcox.test(diversity_df$Shannon[diversity_df$Subtype != "Control"]~diversity_df$Subtype[diversity_df$Subtype != "Control"])
W = 5702, p-value = 0.8702

#confirm findings with rarefied data to control for library size (use sample size of 10000 - loosing 3 samples)
ps.rare <- rarefy_even_depth(ps_phy_noNA, sample.size = 10000, rngseed = 101)
diversity_df <- data.frame(sample_data(ps.rare), estimate_richness(ps.rare))

#K-W test
kruskal.test(diversity_df$Shannon~diversity_df$Subtype)
Kruskal-Wallis chi-squared = 7.6522, df = 2, p-value = 0.02179

#pairwise Wilcoxon
wilcox.test(diversity_df$Shannon[diversity_df$Subtype != "UC"]~diversity_df$Subtype[diversity_df$Subtype != "UC"])
W = 3382, p-value = 0.01488
wilcox.test(diversity_df$Shannon[diversity_df$Subtype != "CD"]~diversity_df$Subtype[diversity_df$Subtype != "CD"])
W = 3281, p-value = 0.008646
wilcox.test(diversity_df$Shannon[diversity_df$Subtype != "Control"]~diversity_df$Subtype[diversity_df$Subtype != "Control"])
W = 5581, p-value = 0.8474


####### Figure 1a #######
my_comparisons <- list(c("Control", "UC"), c("Control", "CD"), c("UC", "CD"))
pdf("alpha_diversity/alpha_div_subtype.pdf", width = 2, height = 2.5)
ggplot(diversity_df, aes(x=Subtype,y=Shannon)) + 
  geom_boxplot(outlier.size = 0.7, fill=c("skyblue1", "goldenrod1", "chocolate3")) + xlab("") + ylab("Shannon Index") + theme_classic2() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + stat_compare_means(comparisons = rev(my_comparisons),label = "p.signif")
dev.off()



####### Alpha diversity comparisons between no neoplasia and neoplasia for different cohorts (control, CD, UC) #######

####### Results #########
#use unfiltered data
diversity_df <- data.frame(sample_data(ps_phy_noNA), estimate_richness(ps_phy_noNA))

#overall Kruskal-Wallis
kruskal.test(diversity_df$Shannon~diversity_df$Subtype_neoplasia)
Kruskal-Wallis chi-squared = 10.298, df = 5, p-value = 0.06723

#get medians and IQRs
aggregate(diversity_df$Shannon, by=list(diversity_df$Subtype_neoplasia), median)

Group.1        x
1 Control.N0 3.576770
2 Control.Nx 3.760782
3      CD.N0 3.439092
4      CD.Nx 3.678170
5      UC.N0 3.419582
6      UC.Nx 3.386761

aggregate(diversity_df$Shannon, by=list(diversity_df$Subtype_neoplasia), IQR)

Group.1         x
1 Control.N0 0.6344798
2 Control.Nx 0.5894159
3      CD.N0 0.8856138
4      CD.Nx 0.7297078
5      UC.N0 0.5197393
6      UC.Nx 0.6493047

#check pairwise Wilcoxon test
wilcox.test(diversity_df$Shannon[diversity_df$Subtype == "Control"]~diversity_df$Subtype_neoplasia[diversity_df$Subtype == "Control"])
W = 253, p-value = 0.2149
wilcox.test(diversity_df$Shannon[diversity_df$Subtype == "CD"]~diversity_df$Subtype_neoplasia[diversity_df$Subtype == "CD"])
W = 417, p-value = 0.7205
wilcox.test(diversity_df$Shannon[diversity_df$Subtype == "UC"]~diversity_df$Subtype_neoplasia[diversity_df$Subtype == "UC"])
W = 1128, p-value = 0.9379

#confirm findings with rarefied data to control for library size (use sample size of 10000 - loosing 3 samples)
ps.rare <- rarefy_even_depth(ps_phy_noNA, sample.size = 10000, rngseed = 101)
diversity_df <- data.frame(sample_data(ps.rare), estimate_richness(ps.rare))

#K-W test
kruskal.test(diversity_df$Shannon~diversity_df$Subtype)
Kruskal-Wallis chi-squared = 7.6522, df = 2, p-value = 0.02179

#pairwise Wilcoxon
wilcox.test(diversity_df$Shannon[diversity_df$Subtype == "Control"]~diversity_df$Subtype_neoplasia[diversity_df$Subtype == "Control"])
W = 237, p-value = 0.1697
wilcox.test(diversity_df$Shannon[diversity_df$Subtype == "CD"]~diversity_df$Subtype_neoplasia[diversity_df$Subtype == "CD"])
W = 418, p-value = 0.7287
wilcox.test(diversity_df$Shannon[diversity_df$Subtype == "UC"]~diversity_df$Subtype_neoplasia[diversity_df$Subtype == "UC"])
W = 1114, p-value = 0.8505

####### Figure 2a #######
diversity_df$Subtype_neoplasia <- recode(diversity_df$Subtype_neoplasia, "Control.N0"="Control:No neoplasia","Control.Nx"="Control:Neoplasia","CD.N0"="CD:No neoplasia","CD.Nx"="CD:Neoplasia", "UC.N0"="UC:No neoplasia","UC.Nx"="UC:Neoplasia")
my_comparisons <- list(c("Control:No neoplasia", "Control:Neoplasia"), c("CD:No neoplasia", "CD:Neoplasia"),c("UC:No neoplasia", "UC:Neoplasia"))
pdf("alpha_diversity/alpha_div_subtype_neoplasia.pdf", width = 2.5, height = 3)
ggplot(diversity_df, aes(x=Subtype_neoplasia,y=Shannon)) + 
  geom_boxplot(outlier.size = 0.7, fill=c("skyblue1", "skyblue1", "goldenrod1", "goldenrod1", "chocolate3", "chocolate3")) + xlab("") + ylab("Shannon Index") + theme_classic2() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + stat_compare_means(comparisons = rev(my_comparisons),label = "p.signif", step.increase = 0)
dev.off()



#################################################################

                ### Beta diversity ###

#################################################################


#TSS normalise data
ps.transf <- transform_sample_counts(ps2, function(x) {x/sum(x)})

#ordination Bray-Curtis/MDS
out.donor.bray <- ordinate(ps.transf, method = "MDS", distance = "bray")

####### Figure 2C #######

#set variable, legend name and palate

bray.df <- cbind(out.donor.bray$vectors[,1:2], sample_data(ps.transf))
eval.1 <- (out.donor.bray$values$Eigenvalues[1]/sum(out.donor.bray$values$Eigenvalues))*100
eval.2 <- (out.donor.bray$values$Eigenvalues[2]/sum(out.donor.bray$values$Eigenvalues))*100


pdf("alpha_diversity/beta_div.pdf", width = 3.3, height = 1.8)
ggplot(bray.df, aes(x=Axis.1, y=Axis.2, colour=Subtype, fill=Subtype)) + 
  geom_point(size=1.2, alpha=1, shape=21) + scale_colour_manual(values = c("skyblue1", "goldenrod1", "chocolate3")) + scale_fill_manual(values = c("skyblue1", "pink", "firebrick1")) + 
  stat_ellipse(level = 0.8) + theme_classic() + 
  xlab(paste("PCoA 1 [", round(eval.1, 2), "%]", sep = "")) + 
  ylab(paste("PCoA 2 [", round(eval.2, 2), "%]", sep = "")) + 
  theme(legend.position = "right",
        axis.line = element_line(colour = 'black', size = 0.56), 
        axis.ticks = element_line(colour = 'black', size = 0.56))
dev.off()

#now check PERMANOVA and add resulting values into plot
metadata <- as(sample_data(ps.transf), "data.frame")
datamat <- as(otu_table(ps.transf), "matrix")
set.seed(101)
adonis(datamat ~ metadata$Subtype, permutations = 999, method = "bray")
#R^2 0.01213 p-value 0.002

metadata.ss <- metadata[metadata$Subtype %in% c("Control", "CD"),]
metadata.ss$Subtype <- as.factor(as.character(metadata.ss$Subtype))
datamat.ss <- datamat[metadata$Subtype %in% c("Control", "CD"),]
all(rownames(datamat.ss) == rownames(metadata.ss))
set.seed(101)
adonis(datamat.ss ~ metadata.ss$Subtype, permutations = 999, method = "bray")
#R^2 0.00956 p-value 0.021

metadata.ss <- metadata[metadata$Subtype %in% c("Control", "UC"),]
metadata.ss$Subtype <- as.factor(as.character(metadata.ss$Subtype))
datamat.ss <- datamat[metadata$Subtype %in% c("Control", "UC"),]
all(rownames(datamat.ss) == rownames(metadata.ss))
set.seed(101)
adonis(datamat.ss ~ metadata.ss$Subtype, permutations = 999, method = "bray")
#R^2 0.01135 p-value 0.007

metadata.ss <- metadata[metadata$Subtype %in% c("CD", "UC"),]
metadata.ss$Subtype <- as.factor(as.character(metadata.ss$Subtype))
datamat.ss <- datamat[metadata$Subtype %in% c("CD", "UC"),]
all(rownames(datamat.ss) == rownames(metadata.ss))
set.seed(101)
adonis(datamat.ss ~ metadata.ss$Subtype, permutations = 999, method = "bray")
#R^2 0.0073 p-value 0.026

> p.adjust(c(0.026,0.007,0.021, 0.002), "fdr")
[1] 0.026 0.014 0.026 0.008







#TSS normalise data
ps.transf <- transform_sample_counts(ps2.ibd, function(x) {x/sum(x)})

#ordination Bray-Curtis/MDS
out.donor.bray <- ordinate(ps.transf, method = "MDS", distance = "bray")

####### Figure 2C #######

#set variable, legend name and palate


bray.df <- cbind(out.donor.bray$vectors[,1:2], sample_data(ps.transf))
eval.1 <- (out.donor.bray$values$Eigenvalues[1]/sum(out.donor.bray$values$Eigenvalues))*100
eval.2 <- (out.donor.bray$values$Eigenvalues[2]/sum(out.donor.bray$values$Eigenvalues))*100


pdf("alpha_diversity/beta_div_dysplasia.pdf", width = 3.3, height = 1.8)
ggplot(bray.df, aes(x=Axis.1, y=Axis.2, colour=Subtype_neoplasia)) + 
  geom_point(size=1.2, alpha=1) + scale_colour_manual("Subtype",values = c("grey", "pink", "black", "firebrick1")) + 
  stat_ellipse(level = 0.8) + theme_classic() + 
  xlab(paste("PCoA 1 [", round(eval.1, 2), "%]", sep = "")) + 
  ylab(paste("PCoA 2 [", round(eval.2, 2), "%]", sep = "")) + 
  theme(legend.position = "right",
        axis.line = element_line(colour = 'black', size = 0.56), 
        axis.ticks = element_line(colour = 'black', size = 0.56))
dev.off()


#now check PERMANOVA and add resulting values into plot
metadata <- as(sample_data(ps.transf), "data.frame")
datamat <- as(otu_table(ps.transf), "matrix")
set.seed(101)
adonis(datamat ~ metadata$Subtype_neoplasia, permutations = 999, method = "bray")
#R^2 0.02147 p-value 0.002

metadata.ss <- metadata[metadata$Subtype_neoplasia %in% c("CD:N0", "CD:Nx"),]
metadata.ss$Subtype_neoplasia <- as.factor(as.character(metadata.ss$Subtype_neoplasia))
datamat.ss <- datamat[metadata$Subtype_neoplasia %in% c("CD:N0", "CD:Nx"),]
all(rownames(datamat.ss) == rownames(metadata.ss))
set.seed(101)
adonis(datamat.ss ~ metadata.ss$Subtype_neoplasia, permutations = 999, method = "bray")
#R^2 0.00847 p-value 0.629

metadata.ss <- metadata[metadata$Subtype_neoplasia %in% c("UC:N0", "UC:Nx"),]
metadata.ss$Subtype_neoplasia <- as.factor(as.character(metadata.ss$Subtype_neoplasia))
datamat.ss <- datamat[metadata$Subtype_neoplasia %in% c("UC:N0", "UC:Nx"),]
all(rownames(datamat.ss) == rownames(metadata.ss))
set.seed(101)
adonis(datamat.ss ~ metadata.ss$Subtype_neoplasia, permutations = 999, method = "bray")
#R^2 0.02066 p-value 0.003

> p.adjust(c(0.002,0.629,0.003), "fdr")
[1] 0.0045 0.6290 0.0045









pdf("alpha_diversity/beta_div_neoplasia.pdf", width = 3.3, height = 1.8)
ggplot(bray.df, aes(x=Axis.1, y=Axis.2, colour=Neoplasia)) + 
  geom_point(size=1.2, alpha=1) + scale_colour_manual("Neoplasia",values = c("blue", "purple", "peachpuff", "salmon", "red1")) + 
  stat_ellipse(level = 0.8) + theme_classic() + 
  xlab(paste("PCoA 1 [", round(eval.1, 2), "%]", sep = "")) + 
  ylab(paste("PCoA 2 [", round(eval.2, 2), "%]", sep = "")) + 
  theme(legend.position = "right",
        axis.line = element_line(colour = 'black', size = 0.56), 
        axis.ticks = element_line(colour = 'black', size = 0.56))
dev.off()


#now check PERMANOVA and add resulting values into plot
metadata <- as(sample_data(ps.transf), "data.frame")
datamat <- as(otu_table(ps.transf), "matrix")
set.seed(101)
adonis(datamat ~ metadata$Neoplasia, permutations = 999, method = "bray")
#R^2 0.02366 p-value 0.015

metadata.ss <- metadata[metadata$Neoplasia %in% c("N0", "N3"),]
metadata.ss$Neoplasia <- as.factor(as.character(metadata.ss$Neoplasia))
datamat.ss <- datamat[metadata$Neoplasia %in% c("N0", "N3"),]
all(rownames(datamat.ss) == rownames(metadata.ss))
set.seed(101)
adonis(datamat.ss ~ metadata.ss$Neoplasia, permutations = 999, method = "bray")
#R^2 0.00782 p-value 0.056



################# Go to DMM script ################





