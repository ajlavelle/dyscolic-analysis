#################################################################

                    ### Figure S1 ###

#################################################################


######alpha diversity#####
library(RColorBrewer)
library(ggpubr)
dir.create(paste(getwd(), "alpha_diversity", sep = "/"))



####### Figure S1b #######
demo_df <- data.frame(sample_data(ps2))
my_comparisons <- list(c("Control", "UC"), c("Control", "CD"), c("UC", "CD"))

pdf("Demographics/age_all.pdf", width = 2, height = 2.7)
ggplot(demo_df, aes(x=Subtype,y=Age,fill=Subtype)) + 
  geom_boxplot(outlier.size = 0.7, fill=c("skyblue1", "goldenrod1", "chocolate3")) + xlab("") + ylab("Age (years)") + theme_classic2() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + stat_compare_means(comparisons = rev(my_comparisons),label = "p.signif")
dev.off()

####### Figure S1c #######
pdf("Demographics/sample_gap.pdf", width = 2, height = 2.7)
ggplot(demo_df, aes(x=Subtype,y=Gap_between_stool_sampling_and_colonoscopy,fill=Subtype)) + 
  geom_boxplot(outlier.size = 0.7, fill=c("skyblue1", "goldenrod1", "chocolate3")) + xlab("") + ylab("Gap between sample &\ncolonoscopy (days/365)") + theme_classic2() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + stat_compare_means(comparisons = rev(my_comparisons),label = "p.signif")
dev.off()


####### Figure S1d #######
demo_df <- data.frame(sample_data(ps2.ibd))
my_comparisons <- list(c("N0", "N3"), c("N0", "N2"), c("N0", "N1"), c("N0", "A1"))

pdf("Demographics/duration_ibd.pdf", width = 2.5, height = 2.7)
ggplot(demo_df, aes(x=Neoplasia,y=Duration,fill=Neoplasia)) + 
  geom_boxplot(outlier.size = 0.7, fill=c("blue", "purple", "peachpuff", "salmon", "red1")) + xlab("") + ylab("Duration (years)") + theme_classic2() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + stat_compare_means(comparisons = rev(my_comparisons),label = "p.signif")
dev.off()

####### Figure S1e #######
demo_df <- data.frame(sample_data(ps2.ibd))
my_comparisons <- list(c("N0", "N3"), c("N0", "N2"), c("N0", "N1"), c("N0", "A1"))

pdf("Demographics/age_ibd.pdf", width = 2.5, height = 2.7)
ggplot(demo_df, aes(x=Neoplasia,y=Age,fill=Neoplasia)) + 
  geom_boxplot(outlier.size = 0.7, fill=c("blue", "purple", "peachpuff", "salmon", "red1")) + xlab("") + ylab("Age (years)") + theme_classic2() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + stat_compare_means(comparisons = rev(my_comparisons),label = "p.signif")
dev.off()


























####### Figure S1f #######

#read in neoplasia location table
neoplasia.df <- read.table("neoplasia_location.txt")

neoplasia.df$Location <- factor(neoplasia.df$Location, levels = c("Right colon", "Transverse colon", "Left colon", "Rectosigmoid"))
neoplasia.df$Neoplasia <- factor(neoplasia.df$Neoplasia, levels = c("Adenoma", "Dysplasia", "Cancer"))
neoplasia.df$Subtype <- factor(neoplasia.df$Subtype, levels = c("Control", "CD", "UC"))


ggplot(neoplasia.df, aes(x=Location, y=Count, fill=Neoplasia)) + geom_bar(stat = "identity") + 
  theme_classic() + facet_grid(Neoplasia~Subtype) + scale_fill_manual(values = c("purple", "salmon2", "red")) + 
  ylab("Count")


####### Figure S1g #######

#get percentages
table(mapping_file$Montreal_combined[mapping_file$Subtype == "CD"])/length(mapping_file$Subtype[mapping_file$Subtype == "CD"])
table(mapping_file$Montreal_combined[mapping_file$Subtype == "UC"])/length(mapping_file$Subtype[mapping_file$Subtype == "UC"])


####### Figure S1h #######
mapping_file <- data.frame(sample_data(ps2))

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

pie.df <- data.frame(matrix(nrow = 2, ncol = 4))
colnames(pie.df) <- c("Control", "CD", "UC", "PSC")
pie.df[1,1] <- length(mapping_file$PSC_yes_no[mapping_file$Subtype == "Control" & mapping_file$PSC_yes_no == "Yes"])
pie.df[2,1] <- length(mapping_file$PSC_yes_no[mapping_file$Subtype == "Control" & mapping_file$PSC_yes_no == "No"])
pie.df[1,2] <- length(mapping_file$PSC_yes_no[mapping_file$Subtype == "CD" & mapping_file$PSC_yes_no == "Yes"])
pie.df[2,2] <- length(mapping_file$PSC_yes_no[mapping_file$Subtype == "CD" & mapping_file$PSC_yes_no == "No"])
pie.df[1,3] <- length(mapping_file$PSC_yes_no[mapping_file$Subtype == "UC" & mapping_file$PSC_yes_no == "Yes"])
pie.df[2,3] <- length(mapping_file$PSC_yes_no[mapping_file$Subtype == "UC" & mapping_file$PSC_yes_no == "No"])
pie.df[1,4] <- "Yes"
pie.df[2,4] <- "No"
pie.df$PSC <- factor(pie.df$PSC, levels = c("Yes", "No"))

psc <- ggplot(pie.df, aes(x="",y=Control, fill=PSC)) + geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("brown3", "green3")) + blank_theme
psc + coord_polar("y", start=0) + theme(axis.text.x=element_blank()) +
  geom_text(aes(y = c(50.5, 25), 
                label = Control), size=5, colour="white")

psc.cd <- ggplot(pie.df, aes(x="",y=CD, fill=PSC)) + geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("brown3", "green3")) + blank_theme
psc.cd + coord_polar("y", start=0) + theme(axis.text.x=element_blank()) +
  geom_text(aes(y = c(102, 50), 
                label = CD), size=5, colour="white")


psc.UC <- ggplot(pie.df, aes(x="",y=UC, fill=PSC)) + geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("brown3", "green3")) + blank_theme
psc.UC + coord_polar("y", start=0) + theme(axis.text.x=element_blank()) +
  geom_text(aes(y = c(90, 35), 
                label = UC), size=5, colour="white")







