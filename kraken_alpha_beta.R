library(stringr)
library(ggvenn)
library(reshape)
library(ggplot2)
library(ggpubr)
library(plotly)
library(plyr)
library(stringr)
library(dplyr)
library(tidyr)
library(plyr)
library(vegan)
library(phyloseq)
library(metagenomeSeq)
library(vroom)
library(xlsx)
library(ggrepel)


files <- list.files(path="./", pattern="*.csv", full.names=TRUE, recursive=FALSE)
files <- files[1:8]

each <- gsub("\\..*", "", gsub(".*\\/", "", files))
each
each[c(1,3,4,6,7)]

head(P_dt)
stat_tab <- c()
#use N_taxon "Number of fragments assigned directly to this taxon"
for(j in 2){#6:(length(files)-1)){c(7,2,6)
  #j = 6
      print(paste(each[j], j))
      X_dt <- read.csv(files[j])
      X_dt$newname <- paste(X_dt$sampleID, X_dt$enzyme, sep = "_")
      OTU_read_count <- reshape(X_dt[, c("newname", "name", "N_taxon")], idvar = "name", timevar = "newname", direction = "wide")
      #head(OTU_read_count)
      #check
      head(X_dt)
      
      row.names(OTU_read_count) <- OTU_read_count$name
      OTU_read_count <- OTU_read_count[, -1]
      metaSeqObject      = newMRexperiment(OTU_read_count)
      metaSeqObject_CSS  = cumNorm( metaSeqObject , p=cumNormStatFast(metaSeqObject) )
      OTU_read_count_CSS = data.frame(MRcounts(metaSeqObject_CSS, norm=TRUE, log=TRUE))
      row.names(OTU_read_count_CSS) <- make.unique(trimws(row.names(OTU_read_count_CSS)))
      
      to_OTU <- OTU_read_count
      for(k in 1:ncol(to_OTU)){
        to_OTU[,k] <- to_OTU[,k]/sum(to_OTU[,k])
      }
      
      krakendt <- OTU_read_count_CSS
      head(krakendt)
      krakendt <- krakendt
      for(k in 1:ncol(krakendt)){
        krakendt[,k] <- krakendt[,k]/sum(krakendt[,k])
      }
      
      #
      write.csv(krakendt, paste0("./", each[j], ".csv"),
                row.names = T)
      
      #plot abundance
      krakendt <- to_OTU
      rownames(krakendt) <- make.unique(gsub(" ", "", rownames(krakendt)))#str_replace_all(rownames(krakendt), fixed(" "), "")
    
      kdf <- krakendt
      kdf$Taxon <- row.names(kdf)
      kdf <- melt(kdf, id.vars = "Taxon", variable.name = "Treatments", value.name = "Abundance")
      colnames(kdf) <- c("Taxon", "Treatments", "Abundance")
      kdf <- kdf[kdf$Abundance != 0,]
      
      kdf$Taxon <- ifelse(kdf$Abundance > quantile(kdf$Abundance, 0.95), kdf$Taxon, "Other")
      ggplot(data = kdf) + 
        geom_bar(aes(x = Treatments, y = Abundance, fill = Taxon), stat = "identity")
      ggsave(paste0("./kraken_abundance_percent_", each[j], ".jpg"), 
             width = 20, height = 20, units = "cm")
      
      krakendt <- OTU_read_count_CSS
      head(krakendt)
      krakendt <- krakendt
      for(k in 1:ncol(krakendt)){
        krakendt[,k] <- krakendt[,k]/sum(krakendt[,k])
      }
      
      kdf <- krakendt
      kdf$Taxon <- row.names(kdf)
      kdf <- melt(kdf, id.vars = "Taxon", variable.name = "Treatments", value.name = "Abundance")
      colnames(kdf) <- c("Taxon", "Treatments", "Abundance")
      kdf <- kdf[kdf$Abundance != 0,]
      
      kdf$Taxon <- ifelse(kdf$Abundance > quantile(kdf$Abundance, 0.95), kdf$Taxon, "Other")
      ggplot(data = kdf) + 
        geom_bar(aes(x = Treatments, y = Abundance, fill = Taxon), stat = "identity")
      ggsave(paste0("./kraken_abundance_CSS_", each[j], ".jpg"), 
             width = 20, height = 20, units = "cm")
      
      #Shannon diversity index (Phyloseq)
      head(krakendt)
      shan <- c()
      for(i in 1:ncol(krakendt)){
        shan <- c(shan, diversity(krakendt[,i], index = "shannon")  )
      }
      
      #Simpson diversity index (Phyloseq)
      simp <- c()
      for(i in 1:ncol(krakendt)){
        simp <- c(simp, diversity(krakendt[,i], index = "simpson")  )
      }
      
      #Inverse Simpson diversity index (Phyloseq)
      invsimp <- c()
      for(i in 1:ncol(krakendt)){
        invsimp <- c(invsimp, diversity(krakendt[,i], index = "invsimpson")  )
      }
      
      #Evenness (Shannon equitability index, Simpson evenness)
      evenn <- c()
      for(i in 1:ncol(krakendt)){
        eachev <- diversity(krakendt[,i], index = "shannon")/(log(specnumber(krakendt)[i]))
        evenn <- c(evenn, eachev  )
      }
      
      #Chao1 (Richness) (Phyloseq)
      
      chao <- c()
      for(i in 1:ncol(krakendt)){
        echao <- specpool(t(krakendt[,i]))
        chao <- c(chao, echao$chao)
      }
      
      krakendiv <- cbind.data.frame(colnames(krakendt), shan, simp, invsimp, evenn, chao)
      colnames(krakendiv)[1] <- "ID"
      row.names(krakendiv) <- NULL
      krakendiv$treatment <- ifelse(!grepl("enzyme", krakendiv$ID), "none", "enzyme")
      krakendiv$treatment <- factor(krakendiv$treatment, levels = c("none", "enzyme"))
      
      #update evenness
      krakendiv$evenn <- krakendiv$shan/log(krakendiv$chao)
      
      lkrakendiv <- krakendiv %>% pivot_longer(cols = c(shan, simp, invsimp, evenn, chao),
                                         names_to="metric")
      
      write.csv(krakendiv, paste0("./kraken_alpha_", each[j], ".csv"))
      write.xlsx(krakendiv, paste0("./kraken_alpha_", each[j], ".xlsx"))
      
      #box plot
      ggboxplot(krakendiv, x = "treatment", y = "simp",
                color = "treatment", add = "dotplot", palette =c("#00AFBB", "#E7B800", "#FC4E07")) + 
        stat_compare_means(method = "t.test") +
        stat_compare_means(label = "p.signif", method = "t.test",
                           ref.group = "none") +
        ggtitle(paste0("Simpson diversity index ", each[j])) +
        theme(plot.title = element_text(hjust = 0.5)) +
        ylab("Relative abundance of Kraken") + xlab("Treatments")
      ggsave(paste0("./kraken_simpson_", each[j], ".jpg"), 
             width = 20, height = 20, units = "cm")
      
      ggboxplot(krakendiv, x = "treatment", y = "invsimp",
                color = "treatment", add = "dotplot", palette =c("#00AFBB", "#E7B800", "#FC4E07")) + 
        stat_compare_means(method = "t.test") +
        stat_compare_means(label = "p.signif", method = "t.test",
                           ref.group = "none") +
        ggtitle(paste0("Inverse Simpson diversity index ", each[j])) +
        theme(plot.title = element_text(hjust = 0.5)) +
        ylab("Relative abundance of Kraken") + xlab("Treatments")
      ggsave(paste0("./kraken_invsimpson_", each[j], ".jpg"), 
             width = 20, height = 20, units = "cm")
      
      ggboxplot(krakendiv, x = "treatment", y = "evenn",
                color = "treatment", add = "dotplot", palette =c("#00AFBB", "#E7B800", "#FC4E07")) + 
        stat_compare_means(method = "t.test") +
        stat_compare_means(label = "p.signif", method = "t.test",
                           ref.group = "none") +
        ggtitle(paste0("Evenness diversity index ", each[j])) +
        theme(plot.title = element_text(hjust = 0.5)) +
        ylab("Relative abundance of Kraken") + xlab("Treatments")
      ggsave(paste0("./kraken_evenness_", each[j], ".jpg"), 
             width = 20, height = 20, units = "cm")
      
      ggboxplot(krakendiv, x = "treatment", y = "chao",
                color = "treatment", add = "dotplot", palette =c("#00AFBB", "#E7B800", "#FC4E07")) + 
        stat_compare_means(method = "t.test") +
        stat_compare_means(label = "p.signif", method = "t.test",
                           ref.group = "none") +
        ggtitle(paste0("Richness ", each[j])) +
        theme(plot.title = element_text(hjust = 0.5)) +
        ylab("Relative abundance of Kraken") + xlab("Treatments")
      ggsave(paste0("./kraken_richness_", each[j], ".jpg"), 
             width = 20, height = 20, units = "cm")
      
      #NMDR in R following https://chrischizinski.github.io/rstats/vegan-ggplot2/ 
      #https://rpubs.com/collnell/manova
      #rownames(krakendt) <- kraken$ID
      #colnames(krakendt) <- gsub("X", "", colnames(krakendt))
      
      varespec <- t(krakendt)
      vare.mds <- metaMDS(varespec)
      
      #vare.mds$stress
      #vare.mds$grstress
      
      #assign group
      data.scores <- scores(vare.mds)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
      data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
      grp <- ifelse(!grepl("enzyme", rownames(varespec)), "none", "enzyme")
      
      krakentoplot <- as.data.frame(data.scores$sites)
      krakentoplot$treatment <- grp
      krakentoplot$ID <- row.names(krakentoplot)
      krakentoplot$ID <- gsub("yme", "", gsub("N_clade.", "", krakentoplot$ID))
      
      find_hull12_id <- function(krakentoplot)krakentoplot[chull(krakentoplot[,1], krakentoplot[,2]), ]
      hulls12_id <- ddply(krakentoplot, "treatment", find_hull12_id)
      
      #permanova
      nmstress <- vare.mds$stress
      lkrakendt <- krakendt
      tlkrakendt <- t(lkrakendt)
      krakencondition <- as.data.frame(row.names(tlkrakendt))
      colnames(krakencondition) <- "format"
      krakencondition <- krakencondition %>% separate(format, sep = "\\.", c("type", "condition")) %>%
        separate(condition, sep = "_", c("sample", "treatment"))
      krakencondition$sample <- factor(krakencondition$sample)
      krakencondition$treatment <- factor(krakencondition$treatment)
      #krakencondition$treatment_num <- ifelse(grepl("non", krakencondition$treatment), 0, 1)
      perm <- adonis2(tlkrakendt ~ treatment, data = krakencondition, permutations = 999, method="bray", by = NULL)
      kraken.dist<-vegdist(tlkrakendt, method='bray')
      #perm <- adonis2(kraken.dist ~ treatment, data = krakencondition, permutations = 9, method="bray")
      pval <- perm$`Pr(>F)`[1]
      r2 <- perm$R2[1]
      
      stat_tab <- rbind.data.frame(stat_tab, cbind.data.frame(each[j], pval, r2, nmstress))
      
      #plot scatter
      ########################
      krakentoplot$samples <- gsub("yme*", "", gsub(".*\\.", "", krakentoplot$ID))
      ggplot(data = krakentoplot, aes(x = NMDS1, y = NMDS2, col=treatment, 
                                   fill=treatment, group=treatment, label = ID)) +
        geom_point(size = 3) + 
        labs(x = "NMDS1", y = "NMDS2") +
        geom_polygon(data = hulls12_id, alpha = 0.3) + theme_bw() +
        theme(text = element_text(size=20)) +
        scale_color_manual("Treatments", breaks = c("none", "enzyme"),
                           values=c("#00AFBB", "#E7B800", "#FC4E07"))  +
        scale_fill_manual("Treatments", breaks = c("none", "enzyme"),
                          values=c("#00AFBB", "#E7B800", "#FC4E07"))  +
        #scale_shape_manual(values=1:nlevels(as.factor(krakentoplot$treatment))) +
        geom_text_repel(max.overlaps = 10000) + labs(title = paste0("Bray–Curtis dissimilarity metric of kraken"))
      
      ggsave(paste0("./kraken_Bray_curtis_group", each[j], ".jpg"), 
             width = 20, height = 20, units = "cm")
      #no frame
      
      ggplot(data = krakentoplot, aes(x = NMDS1, y = NMDS2, col=treatment, 
                                   fill=treatment, group=treatment, label = samples)) +
        geom_point(size = 3) + 
        labs(x = "NMDS1", y = "NMDS2") +
        #geom_polygon(data = hulls12_id, alpha = 0.3) + 
        theme_bw() +
        theme(text = element_text(size=20)) +
        scale_color_manual("Treatments", breaks = c("none", "enzyme"),
                           values=c("#00AFBB", "#E7B800", "#FC4E07"))  +
        scale_fill_manual("Treatments", breaks = c("none", "enzyme"),
                          values=c("#00AFBB", "#E7B800", "#FC4E07"))  +
        #scale_shape_manual(values=1:nlevels(as.factor(krakentoplot$treatment))) +
        geom_text_repel(max.overlaps = 10000) + labs(title = paste0("Bray–Curtis dissimilarity metric of kraken"))
      
      ggsave(paste0("./kraken_Bray_curtis", each[j], ".jpg"), 
             width = 20, height = 20, units = "cm")
      #circle
      ggplot(data = krakentoplot, aes(x = NMDS1, y = NMDS2, col=treatment, 
                                   fill=treatment, group=treatment, label = samples)) +
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
        #scale_shape_manual(values=1:nlevels(as.factor(krakentoplot$treatment))) +
        geom_text_repel(max.overlaps = 10000) + labs(title = paste0("Bray–Curtis dissimilarity metric of kraken"))
      ggsave(paste0("./kraken_Bray_curtis_circle", each[j], ".jpg"), 
             width = 20, height = 20, units = "cm")
      
      write.xlsx(krakentoplot, paste0("./kraken_Bray_curtis_table_", each[j], ".xlsx"))
      
      #venn diagram
      colnames(krakendt)
      sumkrakendt <- cbind.data.frame(krakendt %>% dplyr::select(contains("non")) %>% rowSums(), krakendt %>% dplyr::select(!contains("non")) %>% rowSums())
      head(sumkrakendt)
      colnames(sumkrakendt) <- c("none", "enzyme")
      sumkrakendt$none_name <- ifelse(sumkrakendt$none > 0.00001, rownames(sumkrakendt), NA)
      sumkrakendt$enzyme_name <- ifelse(sumkrakendt$enzyme > 0.00001, rownames(sumkrakendt), NA)
      krakendt_list <- list(
        none_enzyme = na.omit(sumkrakendt$none_name),
        enzyme = na.omit(sumkrakendt$enzyme_name)
      )
      ggvenn(
        krakendt_list, 
        fill_color = c("#00AFBB", "#E7B800", "#FC4E07"),
        stroke_size = 0.5, set_name_size = 4
      ) + labs(title = paste0("Venn diagram kraken ", each[j]))
      ggsave(paste0("./kraken_venn_", each[j], ".jpg"), 
             width = 20, height = 20, units = "cm")
      
      both <- intersect(krakendt_list$none_enzyme, krakendt_list$enzyme)
      only_none <- setdiff(krakendt_list$none_enzyme, krakendt_list$enzyme)
      only_enzyme <- setdiff(krakendt_list$enzyme, krakendt_list$none_enzyme)
      n <- max(max(length(both), length(only_none)), length(only_enzyme))
      length(both) <- n                      
      length(only_none) <- n
      length(only_enzyme) <- n
      
      
      sum_name <- cbind.data.frame(cbind(only_none, both), only_enzyme)
      head(sum_name)
      write.csv(sum_name, paste0("./kraken_summary_name_", each[j], ".csv"), row.names = F,  na = "")
      
      
      
      ###################################################################################
      ###################################################################################
      ###################################################################################
      ###################################################################################
      ###################################################################################
      ###################################################################################
      #positive negative gram CSS
      if(j ==3){
      head(background)
      
      #Shannon diversity index (Phyloseq)
      head(krakendt)
      
      all(background$Family %in% row.names(krakendt))
      updatekrakendt <- krakendt
      updatekrakendt$Family <- row.names(updatekrakendt)
      updatekrakendt$Family <- gsub("\\'", "", updatekrakendt$Family)
      newkrakendt <- merge(updatekrakendt, background[, c("Family", "Gram")], by = "Family", all.x = T, all.y = F)
      newkrakendt$Gram[is.na(newkrakendt$Gram)] <- "unlabelled"
      newkrakendt <- newkrakendt[, !(colnames(newkrakendt) %in% "Family")]
      newkrakendt <- newkrakendt %>% group_by(Gram) %>% summarise_each(list(sum = sum))
      colnames(newkrakendt) <- gsub("_sum", "", colnames(newkrakendt))
      newkrakendt <- as.data.frame(newkrakendt)
      rownames(newkrakendt) <- newkrakendt$Gram
      newkrakendt <- newkrakendt[!(newkrakendt$Gram %in% "unlabelled"),]
      newkrakendt <- newkrakendt[,-1]
      
      
      rownames(krakendt) <- make.unique(gsub(" ", "", rownames(krakendt)))#str_replace_all(rownames(krakendt), fixed(" "), "")
      
      kdf <- newkrakendt
      kdf$Taxon <- row.names(kdf)
      kdf <- melt(kdf, id.vars = "Taxon", variable.name = "Treatments", value.name = "Abundance")
      colnames(kdf) <- c("Taxon", "Treatments", "Abundance")
      kdf <- kdf[kdf$Abundance != 0,]
      kdf$Taxon <- factor(kdf$Taxon, levels = c("Positive", "Negative", "Others"))
      
      #kdf$Taxon <- ifelse(kdf$Abundance > quantile(kdf$Abundance, 0.95), kdf$Taxon, "Other")
      ggplot(data = kdf) + 
        geom_bar(aes(x = Treatments, y = Abundance, fill = Taxon), stat = "identity") +
        scale_fill_manual("Taxon", breaks = c("Positive", "Negative", "Others"),
                          values=c("#00AFBB", "#E7B800", "#FC4E07")) 
      ggsave(paste0("./kraken_abundance_percent_pna", each[j], ".jpg"), 
             width = 20, height = 20, units = "cm")
      
      shan <- c()
      for(i in 1:ncol(newkrakendt)){
        shan <- c(shan, diversity(newkrakendt[,i], index = "shannon")  )
      }
      
      #Simpson diversity index (Phyloseq)
      simp <- c()
      for(i in 1:ncol(newkrakendt)){
        simp <- c(simp, diversity(newkrakendt[,i], index = "simpson")  )
      }
      
      #Inverse Simpson diversity index (Phyloseq)
      invsimp <- c()
      for(i in 1:ncol(newkrakendt)){
        invsimp <- c(invsimp, diversity(newkrakendt[,i], index = "invsimpson")  )
      }
      
      #Evenness (Shannon equitability index, Simpson evenness)
      evenn <- c()
      for(i in 1:ncol(newkrakendt)){
        eachev <- diversity(newkrakendt[,i], index = "shannon")/(log(specnumber(newkrakendt)[i]))
        evenn <- c(evenn, eachev  )
      }
      
      #Chao1 (Richness) (Phyloseq)
      
      chao <- c()
      for(i in 1:ncol(newkrakendt)){
        echao <- specpool(t(newkrakendt[,i]))
        chao <- c(chao, echao$chao)
      }
      
      krakendiv <- cbind.data.frame(colnames(newkrakendt), shan, simp, invsimp, evenn, chao)
      colnames(krakendiv)[1] <- "ID"
      row.names(krakendiv) <- NULL
      krakendiv$treatment <- ifelse(!grepl("enzyme", krakendiv$ID), "none", "enzyme")
      krakendiv$treatment <- factor(krakendiv$treatment, levels = c("none", "enzyme"))
      
      lkrakendiv <- krakendiv %>% pivot_longer(cols = c(shan, simp, invsimp, evenn, chao),
                                               names_to="metric")
      
      write.csv(krakendiv, paste0("./kraken_alpha_pna_", each[j], ".csv"))
      write.xlsx(krakendiv, paste0("./kraken_alpha_pna_", each[j], ".xlsx"))
      
      #box plot
      ggboxplot(krakendiv, x = "treatment", y = "simp",
                color = "treatment", add = "dotplot", palette =c("#00AFBB", "#E7B800", "#FC4E07")) + 
        stat_compare_means(method = "t.test") +
        stat_compare_means(label = "p.signif", method = "t.test",
                           ref.group = "none") +
        ggtitle(paste0("Simpson diversity index ", each[j])) +
        theme(plot.title = element_text(hjust = 0.5)) +
        ylab("Relative abundance of Kraken") + xlab("Treatments")
      ggsave(paste0("./kraken_simpson_pna_", each[j], ".jpg"), 
             width = 20, height = 20, units = "cm")
      
      ggboxplot(krakendiv, x = "treatment", y = "invsimp",
                color = "treatment", add = "dotplot", palette =c("#00AFBB", "#E7B800", "#FC4E07")) + 
        stat_compare_means(method = "t.test") +
        stat_compare_means(label = "p.signif", method = "t.test",
                           ref.group = "none") +
        ggtitle(paste0("Inverse Simpson diversity index ", each[j])) +
        theme(plot.title = element_text(hjust = 0.5)) +
        ylab("Relative abundance of Kraken") + xlab("Treatments")
      ggsave(paste0("./kraken_invsimpson_pna_", each[j], ".jpg"), 
             width = 20, height = 20, units = "cm")
      
      ggboxplot(krakendiv, x = "treatment", y = "evenn",
                color = "treatment", add = "dotplot", palette =c("#00AFBB", "#E7B800", "#FC4E07")) + 
        stat_compare_means(method = "t.test") +
        stat_compare_means(label = "p.signif", method = "t.test",
                           ref.group = "none") +
        ggtitle(paste0("Evenness diversity index ", each[j])) +
        theme(plot.title = element_text(hjust = 0.5)) +
        ylab("Relative abundance of Kraken") + xlab("Treatments")
      ggsave(paste0("./kraken_evenness_pna_", each[j], ".jpg"), 
             width = 20, height = 20, units = "cm")
      
      ggboxplot(krakendiv, x = "treatment", y = "chao",
                color = "treatment", add = "dotplot", palette =c("#00AFBB", "#E7B800", "#FC4E07")) + 
        stat_compare_means(method = "t.test") +
        stat_compare_means(label = "p.signif", method = "t.test",
                           ref.group = "none") +
        ggtitle(paste0("Richness ", each[j])) +
        theme(plot.title = element_text(hjust = 0.5)) +
        ylab("Relative abundance of Kraken") + xlab("Treatments")
      ggsave(paste0("./kraken_richness_pna_", each[j], ".jpg"), 
             width = 20, height = 20, units = "cm")
      
      #NMDR in R following https://chrischizinski.github.io/rstats/vegan-ggplot2/ 
      #https://rpubs.com/collnell/manova
      #rownames(newkrakendt) <- kraken$ID
      #colnames(newkrakendt) <- gsub("X", "", colnames(newkrakendt))
      
      varespec <- t(newkrakendt)
      vare.mds <- metaMDS(varespec)
      
      #vare.mds$stress
      #vare.mds$grstress
      
      #assign group
      data.scores <- scores(vare.mds)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
      data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
      grp <- ifelse(!grepl("enzyme", rownames(varespec)), "none", "enzyme")
      
      krakentoplot <- as.data.frame(data.scores$sites)
      krakentoplot$treatment <- grp
      krakentoplot$ID <- row.names(krakentoplot)
      krakentoplot$ID <- gsub("yme", "", gsub("N_clade.", "", krakentoplot$ID))
      
      find_hull12_id <- function(krakentoplot)krakentoplot[chull(krakentoplot[,1], krakentoplot[,2]), ]
      hulls12_id <- ddply(krakentoplot, "treatment", find_hull12_id)
      
      #permova
      nmstress <- vare.mds$stress
      lnewkrakendt <- newkrakendt
      tlnewkrakendt <- t(lnewkrakendt)
      krakencondition <- as.data.frame(row.names(tlnewkrakendt))
      colnames(krakencondition) <- "format"
      krakencondition <- krakencondition %>% separate(format, sep = "\\.", c("type", "condition")) %>%
        separate(condition, sep = "_", c("sample", "treatment"))
      krakencondition$sample <- factor(krakencondition$sample)
      krakencondition$treatment <- factor(krakencondition$treatment)
      #krakencondition$treatment_num <- ifelse(grepl("non", krakencondition$treatment), 0, 1)
      perm <- adonis2(tlnewkrakendt ~ treatment, data = krakencondition, permutations = 999, method="bray", by = NULL)
      kraken.dist<-vegdist(tlnewkrakendt, method='bray')
      #perm <- adonis2(kraken.dist ~ treatment, data = krakencondition, permutations = 9, method="bray")
      pval <- perm$`Pr(>F)`[1]
      r2 <- perm$R2[1]
      
      stat_add <- cbind.data.frame(paste0(each[j], "_pna"), pval, r2, nmstress)
      colnames(stat_add) <- colnames(stat_tab)
      stat_tab <- rbind.data.frame(stat_tab, stat_add)
      
      #plot scatter
      ########################
      krakentoplot$samples <- gsub("yme*", "", gsub(".*\\.", "", krakentoplot$ID))
      ggplot(data = krakentoplot, aes(x = NMDS1, y = NMDS2, col=treatment, 
                                      fill=treatment, group=treatment, label = ID)) +
        geom_point(size = 3) + 
        labs(x = "NMDS1", y = "NMDS2") +
        geom_polygon(data = hulls12_id, alpha = 0.3) + theme_bw() +
        theme(text = element_text(size=20)) +
        scale_color_manual("Treatments", breaks = c("none", "enzyme"),
                           values=c("#00AFBB", "#E7B800", "#FC4E07"))  +
        scale_fill_manual("Treatments", breaks = c("none", "enzyme"),
                          values=c("#00AFBB", "#E7B800", "#FC4E07"))  +
        #scale_shape_manual(values=1:nlevels(as.factor(krakentoplot$treatment))) +
        geom_text_repel(max.overlaps = 10000) + labs(title = paste0("Bray–Curtis dissimilarity metric of kraken"))
      
      ggsave(paste0("./kraken_Bray_curtis_group_pna_", each[j], ".jpg"), 
             width = 20, height = 20, units = "cm")
      #no frame
      
      ggplot(data = krakentoplot, aes(x = NMDS1, y = NMDS2, col=treatment, 
                                      fill=treatment, group=treatment, label = samples)) +
        geom_point(size = 3) + 
        labs(x = "NMDS1", y = "NMDS2") +
        #geom_polygon(data = hulls12_id, alpha = 0.3) + 
        theme_bw() +
        theme(text = element_text(size=20)) +
        scale_color_manual("Treatments", breaks = c("none", "enzyme"),
                           values=c("#00AFBB", "#E7B800", "#FC4E07"))  +
        scale_fill_manual("Treatments", breaks = c("none", "enzyme"),
                          values=c("#00AFBB", "#E7B800", "#FC4E07"))  +
        #scale_shape_manual(values=1:nlevels(as.factor(krakentoplot$treatment))) +
        geom_text_repel(max.overlaps = 10000) + labs(title = paste0("Bray–Curtis dissimilarity metric of kraken"))
      
      ggsave(paste0("./kraken_Bray_curtis_pna_", each[j], ".jpg"), 
             width = 20, height = 20, units = "cm")
      #circle
      ggplot(data = krakentoplot, aes(x = NMDS1, y = NMDS2, col=treatment, 
                                      fill=treatment, group=treatment, label = samples)) +
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
        #scale_shape_manual(values=1:nlevels(as.factor(krakentoplot$treatment))) +
        geom_text_repel(max.overlaps = 10000) + labs(title = paste0("Bray–Curtis dissimilarity metric of kraken"))
      ggsave(paste0("./kraken_Bray_curtis_circle_pna_", each[j], ".jpg"), 
             width = 20, height = 20, units = "cm")
      
      write.xlsx(krakentoplot, paste0("./kraken_Bray_curtis_table_pna_", each[j], ".xlsx"))
      
      #venn diagram
      colnames(newkrakendt)
      sumnewkrakendt <- cbind.data.frame(newkrakendt %>% dplyr::select(contains("non")) %>% rowSums(), newkrakendt %>% dplyr::select(!contains("non")) %>% rowSums())
      head(sumnewkrakendt)
      colnames(sumnewkrakendt) <- c("none", "enzyme")
      sumnewkrakendt$none_name <- ifelse(sumnewkrakendt$none > 0.00001, rownames(sumnewkrakendt), NA)
      sumnewkrakendt$enzyme_name <- ifelse(sumnewkrakendt$enzyme > 0.00001, rownames(sumnewkrakendt), NA)
      newkrakendt_list <- list(
        none_enzyme = na.omit(sumnewkrakendt$none_name),
        enzyme = na.omit(sumnewkrakendt$enzyme_name)
      )
      ggvenn(
        newkrakendt_list, 
        fill_color = c("#00AFBB", "#E7B800", "#FC4E07"),
        stroke_size = 0.5, set_name_size = 4
      ) + labs(title = paste0("Venn diagram kraken ", each[j]))
      ggsave(paste0("./kraken_venn_pna_", each[j], ".jpg"), 
             width = 20, height = 20, units = "cm")
      
      both <- intersect(newkrakendt_list$none_enzyme, newkrakendt_list$enzyme)
      only_none <- setdiff(newkrakendt_list$none_enzyme, newkrakendt_list$enzyme)
      only_enzyme <- setdiff(newkrakendt_list$enzyme, newkrakendt_list$none_enzyme)
      n <- max(max(length(both), length(only_none)), length(only_enzyme))
      length(both) <- n                      
      length(only_none) <- n
      length(only_enzyme) <- n
      
      
      sum_name <- cbind.data.frame(cbind(only_none, both), only_enzyme)
      head(sum_name)
      write.csv(sum_name, paste0("./kraken_summary_name_pna_", each[j], ".csv"), row.names = F,  na = "")
      }
      
}
colnames(stat_tab)[1] <- "case"
write.csv(stat_tab, paste0("./kraken_stat.csv"), row.names = F,  na = "")      
