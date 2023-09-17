library(ggvenn)
library(ggrepel)
library(ggplot2)
library(ggpubr)
library(plotly)
library(plyr)
library(stringr)
library(dplyr)
library(tidyr)
library(vegan)
library(phyloseq)
library(openxlsx)

amr <- read.xlsx("./Sum_AMR_analytic_matrix_normalised_updated.xlsx", sheet = "Normalized_AMR_analytic_matrix")

amrdt <- amr[,8:ncol(amr)]
#Shannon diversity index (Phyloseq)
head(amrdt)
shan <- c()
for(i in 1:ncol(amrdt)){
  shan <- c(shan, diversity(amrdt[,i], index = "shannon")  )
}

#Simpson diversity index (Phyloseq)
simp <- c()
for(i in 1:ncol(amrdt)){
  simp <- c(simp, diversity(amrdt[,i], index = "simpson")  )
}

#Inverse Simpson diversity index (Phyloseq)
invsimp <- c()
for(i in 1:ncol(amrdt)){
  invsimp <- c(invsimp, 1/diversity(amrdt[,i], index = "simpson")  )
}
#diversity(amrdt[,i], index = "invsimpson")

#Evenness (Shannon equitability index, Simpson evenness)
evenn <- c()
for(i in 1:ncol(amrdt)){
  eachev <- diversity(amrdt[,i], index = "shannon")/(log(specnumber(amrdt)[i]))
  evenn <- c(evenn, eachev  )
}



#Chao1 (Richness) (Phyloseq)

chao <- c()
for(i in 1:ncol(amrdt)){
  echao <- specpool(t(amrdt[,i]))
  chao <- c(chao, echao$chao)
}

amrdiv <- cbind.data.frame(colnames(amrdt), shan, simp, invsimp, evenn, chao)

#update evenness
amrdiv$evenn <- amrdiv$shan/log(amrdiv$chao)

colnames(amrdiv)[1] <- "ID"
amrdiv$ID <- gsub("X", "", amrdiv$ID)
amrdiv$treatment <- ifelse(grepl("N", amrdiv$ID), "none", "enzyme")
amrdiv$treatment <- factor(amrdiv$treatment, levels = c("none", "enzyme"))

lamrdiv <- amrdiv %>% pivot_longer(cols = c(shan, simp, invsimp, evenn, chao),
                        names_to="metric")


#box plot
ggboxplot(amrdiv, x = "treatment", y = "shan",
          color = "treatment", add = "dotplot") + 
  stat_compare_means(method = "t.test") +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "none") +
  ggtitle("Shanon") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Relative abundance of AMR") + xlab("Treatments")
ggsave("./AMR_shan.jpg", 
       width = 20, height = 20, units = "cm")

ggboxplot(amrdiv, x = "treatment", y = "simp",
          color = "treatment", add = "dotplot") + 
  stat_compare_means(method = "t.test") +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "none") +
  ggtitle("Simpson diversity index") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Relative abundance of AMR") + xlab("Treatments")
ggsave("./AMR_simpson.jpg", 
       width = 20, height = 20, units = "cm")

ggboxplot(amrdiv, x = "treatment", y = "invsimp",
          color = "treatment", add = "dotplot", palette =c("#00AFBB", "#E7B800", "#FC4E07")) + 
  stat_compare_means(method = "t.test") +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "none") +
  ggtitle("Inverse Simpson diversity index") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Relative abundance of AMR") + xlab("Treatments")
ggsave("./AMR_invsimpson.jpg", 
       width = 20, height = 20, units = "cm")

ggboxplot(amrdiv, x = "treatment", y = "evenn",
          color = "treatment", add = "dotplot", palette =c("#00AFBB", "#E7B800", "#FC4E07")) + 
  stat_compare_means(method = "t.test") +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "none") +
  ggtitle("Evenness diversity index") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Relative abundance of AMR") + xlab("Treatments")
ggsave("./AMR_evenness.jpg", 
       width = 20, height = 20, units = "cm")

ggboxplot(amrdiv, x = "treatment", y = "chao",
          color = "treatment", add = "dotplot", palette =c("#00AFBB", "#E7B800", "#FC4E07")) + 
  stat_compare_means(method = "t.test") +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "none") +
  ggtitle("Richness") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Relative abundance of AMR") + xlab("Treatments")
ggsave("./AMR_richness.jpg", 
       width = 20, height = 20, units = "cm")

#NMDR
#amrdt <- amrdt[!is.na(amrdt$ID),]

rownames(amrdt) <- make.names(amr$ID, unique = T)
colnames(amrdt) <- gsub("X", "", colnames(amrdt))

varespec <- t(amrdt)
vare.mds <- metaMDS(varespec)

#permova
nmstress <- vare.mds$stress
lamrdt <- amrdt
tlamrdt <- t(lamrdt)
amrcondition <- as.data.frame(row.names(tlamrdt))
colnames(amrcondition) <- "format"
amrcondition <- amrcondition %>% mutate(sample = gsub("N", "", format), 
                                        treatment = ifelse(grepl("N", format), "none", "enzyme"))
amrcondition$sample <- factor(amrcondition$sample)
amrcondition$treatment <- factor(amrcondition$treatment)
#amrcondition$treatment_num <- ifelse(grepl("non", amrcondition$treatment), 0, 1)
perm <- adonis2(tlamrdt ~ treatment, data = amrcondition, permutations = 999, method="bray", by = NULL)
amr.dist<-vegdist(tlamrdt, method='bray')
#perm <- adonis2(amr.dist ~ treatment, data = amrcondition, permutations = 9, method="bray")
pval <- perm$`Pr(>F)`[1]
r2 <- perm$R2[1]
stat_tab <- cbind.data.frame("AMR", pval, r2, nmstress)
write.csv(stat_tab, paste0("/Users/ksongsom/Library/CloudStorage/OneDrive-Personal/SideHustle/Ball/update_result_Aug2023/amr_stat.csv"), row.names = F,  na = "", quote = F)

#assign group
data.scores <- scores(vare.mds)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
grp <- ifelse(grepl("N", rownames(varespec)), "none", "enzyme")

amrtoplot <- as.data.frame(data.scores$sites)
amrtoplot$treatment <- grp
amrtoplot$ID <- row.names(amrtoplot)
amrtoplot$ID <- ifelse(grepl("N", amrtoplot$ID), paste0(gsub("N", "", amrtoplot$ID), "_none"), paste0(amrtoplot$ID, "_enz"))

find_hull12_id <- function(amrtoplot)amrtoplot[chull(amrtoplot[,1], amrtoplot[,2]), ]
hulls12_id <- ddply(amrtoplot, "treatment", find_hull12_id)

amrtoplot$samples <- paste(gsub("N", "", amrtoplot$ID), amrtoplot$treatment)
ggplot(data = amrtoplot, aes(x = NMDS1, y = NMDS2, col=treatment, 
                                         fill=treatment, group=treatment, label = ID)) +
  geom_point(size = 3) + 
  labs(x = "NMDS1", y = "NMDS2") +
  geom_polygon(data = hulls12_id, alpha = 0.3) + theme_bw() +
  theme(text = element_text(size=20)) +
  scale_color_manual("Treatments", breaks = c("none", "enzyme"),
                     values=c("#00AFBB", "#E7B800", "#FC4E07"))  +
  scale_fill_manual("Treatments", breaks = c("none", "enzyme"),
                    values=c("#00AFBB", "#E7B800", "#FC4E07"))  +
  #scale_shape_manual(values=1:nlevels(as.factor(amrtoplot$treatment))) +
  geom_text_repel(max.overlaps = 10000) + labs(title = paste0("Bray–Curtis dissimilarity metric of AMR"))

ggsave("./AMR_Bray_curtis_group.jpg", 
       width = 20, height = 20, units = "cm")
#no frame

ggplot(data = amrtoplot, aes(x = NMDS1, y = NMDS2, col=treatment, 
                             fill=treatment, group=treatment, label = ID)) +
  geom_point(size = 3) + 
  labs(x = "NMDS1", y = "NMDS2") +
  #geom_polygon(data = hulls12_id, alpha = 0.3) + 
  theme_bw() +
  theme(text = element_text(size=20)) +
  scale_color_manual("Treatments", breaks = c("none", "enzyme"),
                     values=c("#00AFBB", "#E7B800", "#FC4E07"))  +
  scale_fill_manual("Treatments", breaks = c("none", "enzyme"),
                    values=c("#00AFBB", "#E7B800", "#FC4E07"))  +
  #scale_shape_manual(values=1:nlevels(as.factor(amrtoplot$treatment))) +
  geom_text_repel(max.overlaps = 10000) + labs(title = paste0("Bray–Curtis dissimilarity metric of AMR"))

ggsave("./AMR_Bray_curtis.jpg", 
       width = 20, height = 20, units = "cm")
#circle
ggplot(data = amrtoplot, aes(x = NMDS1, y = NMDS2, col=treatment, 
                             fill=treatment, group=treatment, label = ID)) +
  geom_point(size = 3) + 
  labs(x = "NMDS1", y = "NMDS2") +
  #geom_polygon(data = hulls12_id, alpha = 0.3) + 
  theme_bw() +
  theme(text = element_text(size=20)) +
  stat_ellipse() +
  scale_color_manual("Treatments", breaks = c("none", "enzyme"),
                     values=c("#00AFBB", "#E7B800", "#FC4E07"))  +
  scale_fill_manual("Treatments", breaks = c("none", "enzyme"),
                    values=c("#00AFBB", "#E7B800", "#FC4E07"))  +
  #scale_shape_manual(values=1:nlevels(as.factor(amrtoplot$treatment))) +
  geom_text_repel(max.overlaps = 10000) + labs(title = paste0("Bray–Curtis dissimilarity metric of AMR"))
ggsave("./AMR_Bray_curtis_circle.jpg", 
       width = 20, height = 20, units = "cm")

write.xlsx(amrtoplot, "./AMR_Bray_curtis_table.xlsx")

#venn
head(amr)
sumamrdt <- cbind.data.frame(amrdt %>% dplyr::select(contains("N")) %>% rowSums(), amrdt %>% dplyr::select(!contains("N")) %>% rowSums())
head(sumamrdt)
colnames(sumamrdt) <- c("none", "enzyme")
#rownames(sumamrdt) <- amr$Genes

sumamrdt$none_name <- ifelse(sumamrdt$none > 0.00001, rownames(sumamrdt), NA)
sumamrdt$enzyme_name <- ifelse(sumamrdt$enzyme > 0.00001, rownames(sumamrdt), NA)
amrdt_list <- list(
  none_enzyme = na.omit(sumamrdt$none_name),
  enzyme = na.omit(sumamrdt$enzyme_name)
)
ggvenn(
  amrdt_list, 
  fill_color = c("#00AFBB", "#E7B800", "#FC4E07"),
  stroke_size = 0.5, set_name_size = 4
) + labs(title = paste0("Venn diagram amr "))
ggsave(paste0("./amr_venn", ".jpg"), 
       width = 20, height = 20, units = "cm")

both <- intersect(amrdt_list$none_enzyme, amrdt_list$enzyme)
only_none <- setdiff(amrdt_list$none_enzyme, amrdt_list$enzyme)
only_enzyme <- setdiff(amrdt_list$enzyme, amrdt_list$none_enzyme)
n <- max(max(length(both), length(only_none)), length(only_enzyme))
length(both) <- n                      
length(only_none) <- n
length(only_enzyme) <- n


sum_name <- cbind.data.frame(cbind(only_none, both), only_enzyme)
head(sum_name)
write.csv(sum_name, paste0("./amr_summary_name", ".csv"), row.names = F,  na = "")


#group box plot
lamrdiv <- as.data.frame(lamrdiv)
lamrdiv$metric <- factor(lamrdiv$metric, levels = c("shan", "simp", "invsimp", "evenn", "chao"))

#only shan and evenness
ggboxplot(lamrdiv[lamrdiv$metric %in% c("shan", "evenn"),], x = "treatment", y = "value",
          color = "treatment", add = "dotplot",
          facet.by = "metric", short.panel.labs = T, 
          panel.labs = list(metric = c("Shannon diversity index",
                                       "Simpson diversity index",
                                       "Inverse Simpson diversity index",
                                       "Evenness", "Richness"))) + 
  stat_compare_means(method = "t.test", label.y = max(amrdiv$shan), label.x = 1) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "none") +
  ggtitle("Relative abundance of AMR") +
  theme(plot.title = element_text(hjust = 0.5)) + ylab("Relative abundance of AMR")

ggboxplot(lamrdiv[lamrdiv$metric %in% c("shan", "chao"),], x = "treatment", y = "value",
          color = "treatment", add = "dotplot",
          facet.by = "metric", short.panel.labs = T, 
          panel.labs = list(metric = c("Shannon diversity index",
                                       "Simpson diversity index",
                                       "Inverse Simpson diversity index",
                                       "Evenness", "Richness"))) + 
  stat_compare_means(method = "t.test", label.y = max(amrdiv$shan), label.x = 1) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "none") +
  ggtitle("Relative abundance of AMR") +
  theme(plot.title = element_text(hjust = 0.5)) + ylab("Relative abundance of AMR")


#all five
ggboxplot(lamrdiv, x = "treatment", y = "value",
          color = "treatment", add = "dotplot",
          facet.by = "metric", short.panel.labs = T, 
          panel.labs = list(metric = c("Shannon diversity index",
                                       "Simpson diversity index",
                                       "Inverse Simpson diversity index",
                                       "Evenness", "Richness"))) + 
  stat_compare_means(method = "t.test", label.y = max(amrdiv$shan), label.x = 1) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "none") +
  ggtitle("Relative abundance of AMR") +
  theme(plot.title = element_text(hjust = 0.5)) + ylab("Relative abundance of AMR") + 
  facet_grid(metric ~ treatment, scales='free')


