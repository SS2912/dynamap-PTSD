################################################################################
#################   PTSD - Functional Connectivity Stats   #####################
#############   SS -  started August 2023, last edit: 17/10/23   ###############
################################################################################

# 1. upload packages and functions
rm(list=ls())
library(lme4)
library(ggplot2)
library(grid)
library(lmerTest) # show p values in anova of lmer
library(reshape2)
library(wesanderson) # palette
library(paletteer) # 
library(RColorBrewer)
library(emmeans) # posthoc test avec ajustement des pvalues
library(tidyverse)
library(gapminder)
library(MuMIn)
library(ggpubr)
library(readxl)
library("stringr")
library(rstatix)
#library(gghalves)

source(paste("C:/Users/simula/OneDrive - Aix-Marseille Université/PhD/MyProjects/R", '/summary_data.R', sep = ""), 
       encoding = 'UTF-8', echo=TRUE) # charge cette fonction pour les stats descriptives


############################# 2. read the datasets  #############################    
setwd("//dynaserv/meg/nicolas/PTSD/analysis/FC/FC_tables/With_repetitions")
path_data <- getwd()
path_RES <- "//dynaserv/meg/nicolas/PTSD/analysis/FC/Rstats/WIth_repetitions/Inv-inv" 

bands <- c('delta', 'theta', 'alpha', 'beta', 'lowgamma') #"broad"
band <- "broad"
raw <- read_excel(paste(path_data, "FCtable-inv_ni-alpha-AMYHPCmean-29-Mar-2024.xlsx", sep = "/"), sheet = 1, col_names = TRUE)

for (ii in 1:length(bands)) {
  band <- bands[ii]
  rm(raw)
  rm(fc)
  raw <- read_excel(paste(path_data, paste("FCtable-inv_ni", band, "with_repetitions-22-Mar-2024.xlsx", sep = "-"), sep = "/"), sheet = 1, col_names = TRUE)
  

#############    define groups and pre-process data    #########################
  thr_GAD <- 7    # anxiety
  thr_NDDI <- 15  # depression in epilepsy
  thr_PCL <- 31   # PTSD score
  
  raw <- cbind(raw, PTSDgroup=NA)
  raw$PTSDgroup[raw$PCL_5 >= thr_PCL & is.na(raw$PCL_5) == FALSE] <- "PTSD+"
  raw$PTSDgroup[raw$PCL_5 < thr_PCL & is.na(raw$PCL_5) == FALSE] <- "PTSD-"
  
  raw <- cbind(raw, GADgroup=NA)
  raw$GADgroup[raw$GAD_7 >= thr_GAD & is.na(raw$GAD_7) == FALSE] <- "anxiety"
  raw$GADgroup[raw$GAD_7 < thr_GAD & is.na(raw$GAD_7) == FALSE] <- "no_anxiety"
  
  raw <- cbind(raw, DEPRgroup=NA)
  raw$DEPRgroup[raw$NDDI_E >= thr_NDDI & is.na(raw$NDDI_E) == FALSE] <- "depression"
  raw$DEPRgroup[raw$NDDI_E < thr_NDDI & is.na(raw$NDDI_E) == FALSE] <- "no_depression"

  ##### melt table
  psychosubj <- "9af5a2bd4d70"
  fc <- raw[raw$subject != psychosubj,]
  fc$involved[str_detect(fc$brain_area, "Anterior-cingulate-cortex")] <- "inv"
  my_colnames <- colnames(fc)    
  my_colnames
  
  fc <- melt(as.data.frame(fc), id=my_colnames[c(1:17,28:31)], measure =my_colnames[18:27], value.name = "node_strength", variable.name = "chan_selection")
  # only TOT: fc <- melt(as.data.frame(fc), id=my_colnames[c(1:17,28:31)], measure =my_colnames[seq(19, 27, by = 2)], value.name = "node_strength", variable.name = "chan_selection")
  colnames(fc)  
  fc[, c(9:17,23)] <- lapply(fc[, c(9:17,23)], as.numeric)
  fc[, c(1:8,18:22)] <- lapply(fc[, c(1:8,18:22)], as.factor)
  
  # Filter for "inv" channels
  inv_channels <- fc %>% filter(involved == "inv")
  inv_channels <- inv_channels[inv_channels$chan_selection != "TOTrest_ni_ni", ]
  
  # filter for ipsilateral channels
  ipsiL_channels <- fc %>% filter(EpilepsySide == "L" & str_detect(brain_area, "Left"))
  ipsiR_channels <- fc %>% filter(EpilepsySide == "R" & str_detect(brain_area, "Right"))
  bilat_channels <- fc %>% filter(EpilepsySide == "Bilateral")
  ipsi_all <- rbind(ipsiL_channels, ipsiR_channels, bilat_channels)
  
  # Filter for contralateral channels
  contraL_channels <- fc %>% filter(EpilepsySide == "L" & str_detect(brain_area, "Right"))
  contraR_channels <- fc %>% filter(EpilepsySide == "R" & str_detect(brain_area, "Left"))
  contra_all <- rbind(contraL_channels, contraR_channels)
  
  
  #############  correlation with indexes   #############
  ipsi_all <- cbind(ipsi_all, ipsi="ipsi")
  contra_all <- cbind(contra_all, ipsi="contra")
  all <- rbind(ipsi_all, contra_all)
  allregions <- c("Orbito-frontal-cortex", "Insula", "Amygdala", "Hippocampus", "Rhinal-cortex", "Anterior-cingulate-cortex")
  
  
  for (i in 1:length(allregions)) {
    cur_region <- allregions[i]
    
    sel <- all[str_detect(fc$brain_area, cur_region),]
   
    
    tmp <- sel[sel$chan_selection == "TOTrest_all",]
    my_colnames <- colnames(sel)
    scores <- melt(as.data.frame(tmp), id=my_colnames[c(1:8,18:24)], measure =my_colnames[c(9:17)], value.name = "score_value", variable.name = "score_name")
    
    ymin <- -0.22
    ymax <- 0.22
    
    titre = paste("Correlation of node strength with different scores - ", cur_region)
    AEcorr <- ggplot(data = scores, aes(x = score_value, y = node_strength)) + 
      xlab("score value") + ylab("Node strength h2") + ggtitle(titre)+
      geom_smooth(method = "lm", se = T, colour="#71AA8E", linewidth = 0.5) +
      facet_wrap(score_name~ipsi, ncol =6) +
      theme_light() +
      theme(plot.title = element_text(size=15, hjust = 0.5),
            axis.title = element_text(size = 12, face = "bold"), 
            axis.text = element_text(colour = "black", size = 12),
            strip.text = element_text(size = 12), strip.background = element_rect(fill = "skyblue2")) + 
      geom_point(size = 0.4) +
      stat_cor(method = "spearman", label.x = 0, label.y = -0.08, p.accuracy = 0.001, r.accuracy = 0.01, size = 4)+
      scale_y_continuous(breaks = seq(ymin, ymax, 0.1), limits = c(ymin, ymax))
    
    print(AEcorr)
    # png(paste(path_RES, "Correlation_with_scores", paste(cur_region, ".png", sep = ""), sep = "/"), width = 30, height = 15, units = "cm",
    #     res = 300)
    # print(AEcorr)
    # dev.off()
    
  }
  
  ############################ ipsi and contra
  colmag <- 9
  
  allregions <- c("Hippocampus", "Amygdala")
  
  for (i in 1:length(allregions)) {
    cur_region <- allregions[i]
    rm(pwc)
    rm(cur_data)
    rm(sel)
    
    #### IPSILATERAL ######
    sel <- all
    
    if (cur_region == "all channels") {
      cur_data <- sel
    } else{
      cur_data <- sel %>%
        filter(str_detect(brain_area, cur_region))
      
    }
    
    cur_data <- cur_data[is.na(cur_data$node_strength) == FALSE,]
    ymax = round(max(cur_data$node_strength), 2) +0.03
    ymin = 0
    
    pwc <- cur_data %>% 
      group_by(chan_selection, ipsi) %>%
      pairwise_wilcox_test(node_strength ~ PTSDgroup, p.adjust.method = "fdr", paired = FALSE)
    pwc <- pwc %>% add_xy_position(scales = "free_y")
    
    pwc$y.position <- rep(ymax-0.02, nrow(pwc))
    coco <- cur_data %>%  #|d| < 0.2 "negligible", |d| < 0.5 "small", |d| < 0.8 "medium", otherwise "large"
      group_by(chan_selection, ipsi) %>%
      cohens_d(node_strength ~ PTSDgroup, paired = F)  
    pwc <- cbind(pwc, coco[,4], coco[,colmag])
    pwc$p <- round(pwc$p, 3)
    pwc$effsize <- round(pwc$effsize, 1)
    
    
    titre = paste("Ipsi- and contra-lateral", cur_region)
    allgraph <- ggboxplot(cur_data, x = "PTSDgroup", y = "node_strength",
                           fill = "PTSDgroup", facet.by = c("chan_selection", "ipsi"), outlier.shape = NA) +
      #geom_line(aes(group = channel), alpha = 0.5, colour = "grey") +
      geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 0.7, alpha = 0.6, colour = 'gray20')+
     stat_pvalue_manual(pwc, hide.ns = FALSE, label =  "{p.adj.signif}: p={p.adj}, d={effsize}", bracket.nudge.y = 0, tip.length = 0.005) +
      xlab("group") + ylab("Node strength") + ggtitle(titre) +
      theme_classic() + 
      theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
      theme(axis.title = element_text(size = 10, face = "bold"), 
            axis.text = element_text(colour = "black", size = 10),
            axis.text.x = element_blank()) +
      theme(legend.text = element_text(size = 10),
            legend.title=element_blank(),
            legend.margin = margin(1,1,1,1, unit = "mm"),
            legend.direction = "horizontal",
            legend.position = "bottom")+
      scale_y_continuous(breaks = seq(ymin, ymax, 0.05), limits = c(ymin, ymax)) 
    # scale_fill_manual(values = c( "#cce3de", "#f6d5f7"))
    print(allgraph)
    
  }
  
  #################### REPETAED REGIONS GRAPHS ##########################
  ####################  All inv pooled together  ########################
  colmag <- 8
  
  allregions <- c("all channels", "Orbito-frontal-cortex", "Amygdala", "Hippocampus", "Anterior-cingulate-cortex", "Rhinal-cortex", "Insula")
  #allregions <- c("all channels", "Orbito-frontal-cortex", "Amygdala", "Hippocampus", "Anterior-cingulate-cortex", "Rhinal-cortex")
  # allregions <- "Insula"
  
  for (i in 1:length(allregions)) {
    cur_region <- allregions[i]
    rm(pwc)
    rm(cur_data)
    rm(sel)
  
  #### IPSILATERAL ######
    sel <- ipsiR_channels
    
    if (cur_region == "all channels") {
      cur_data <- sel
    } else{
      cur_data <- sel %>%
        filter(str_detect(brain_area, cur_region))
      
    }
    
    cur_data <- cur_data[is.na(cur_data$node_strength) == FALSE,]
    ymax = round(max(cur_data$node_strength), 2) +0.03
    ymin = 0
    
    pwc <- cur_data %>% 
      group_by(chan_selection) %>%
      pairwise_wilcox_test(node_strength ~ PTSDgroup, p.adjust.method = "fdr", paired = FALSE)
    pwc <- pwc %>% add_xy_position(scales = "free_y")
    
    pwc$y.position <- rep(ymax-0.02, nrow(pwc))
    coco <- cur_data %>%  #|d| < 0.2 "negligible", |d| < 0.5 "small", |d| < 0.8 "medium", otherwise "large"
      group_by(chan_selection) %>%
      cohens_d(node_strength ~ PTSDgroup, paired = F)  
    pwc <- cbind(pwc, coco[,4], coco[,colmag])
    pwc$p <- round(pwc$p, 3)
    pwc$effsize <- round(pwc$effsize, 1)
    
    
    titre = paste("Ipsilateral", cur_region)
    ipsigraph <- ggboxplot(cur_data, x = "PTSDgroup", y = "node_strength",
                          fill = "PTSDgroup", facet.by = c("chan_selection"), outlier.shape = NA) +
      #geom_line(aes(group = channel), alpha = 0.5, colour = "grey") +
      geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 0.7, alpha = 0.6, colour = 'gray20')+
      stat_pvalue_manual(pwc, hide.ns = FALSE, label =  "{p.adj.signif}: p={p.adj}, d={effsize}", bracket.nudge.y = 0, tip.length = 0.005) +
      xlab("group") + ylab("Node strength") + ggtitle(titre) +
      theme_classic() + 
      theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
      theme(axis.title = element_text(size = 10, face = "bold"), 
            axis.text = element_text(colour = "black", size = 10),
            axis.text.x = element_blank()) +
      theme(legend.text = element_text(size = 10),
            legend.title=element_blank(),
            legend.margin = margin(1,1,1,1, unit = "mm"),
            legend.direction = "horizontal",
            legend.position = "bottom")+
      scale_y_continuous(breaks = seq(ymin, ymax, 0.05), limits = c(ymin, ymax)) 
    # scale_fill_manual(values = c( "#cce3de", "#f6d5f7"))
    print(ipsigraph)
    
    ##### CONTRALATERAL ######
    rm(pwc)
    rm(cur_data)
    rm(sel)

    sel <- contraR_channels

    if (cur_region == "all channels") {
      cur_data <- sel
    } else{
      cur_data <- sel %>%
        filter(str_detect(brain_area, cur_region))
    }

    if ((length(unique(cur_data$PTSDgroup)) > 1) & (nrow(cur_data)>6)) {

        cur_data <- cur_data[is.na(cur_data$node_strength) == FALSE,]
        ymax = round(max(cur_data$node_strength), 2) +0.03
        ymin = 0

        pwc <- cur_data %>%
          group_by(chan_selection) %>%
          pairwise_wilcox_test(node_strength ~ PTSDgroup, p.adjust.method = "fdr", paired = FALSE)
        pwc <- pwc %>% add_xy_position(scales = "free_y")

        pwc$y.position <- rep(ymax-0.02, nrow(pwc))
        coco <- cur_data %>%  #|d| < 0.2 "negligible", |d| < 0.5 "small", |d| < 0.8 "medium", otherwise "large"
          group_by(chan_selection) %>%
          cohens_d(node_strength ~ PTSDgroup, paired = F)
        pwc <- cbind(pwc, coco[,4], coco[,colmag])
        pwc$p <- round(pwc$p, 3)
        pwc$effsize <- round(pwc$effsize, 1)


        titre = paste("Contralateral", cur_region)
        contragraph <- ggboxplot(cur_data, x = "PTSDgroup", y = "node_strength",
                              fill = "PTSDgroup", facet.by = c("chan_selection"), outlier.shape = NA) +
          #geom_line(aes(group = channel), alpha = 0.5, colour = "grey") +
          geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 0.7, alpha = 0.6, colour = 'gray20')+
          stat_pvalue_manual(pwc, hide.ns = FALSE, label =  "{p.adj.signif}: p={p.adj}, d={effsize}", bracket.nudge.y = 0, tip.length = 0.005) +
          xlab("group") + ylab("Node strength") + ggtitle(titre) +
          theme_classic() +
          theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
          theme(axis.title = element_text(size = 10, face = "bold"),
                axis.text = element_text(colour = "black", size = 10),
                axis.text.x = element_blank()) +
          theme(legend.text = element_text(size = 10),
                legend.title=element_blank(),
                legend.margin = margin(1,1,1,1, unit = "mm"),
                legend.direction = "horizontal",
                legend.position = "bottom")+
          scale_y_continuous(breaks = seq(ymin, ymax, 0.05), limits = c(ymin, ymax))
        # scale_fill_manual(values = c( "#cce3de", "#f6d5f7"))
        print(contragraph)


      #### PLOT ALL TOGETHER
      figure <- ggarrange(ipsigraph, contragraph,  nrow = 2, labels = c("A", "B"),
                          common.legend = TRUE, legend = "bottom")
    
    }

    else {
        figure <- ipsigraph
        
    }

    title_fig <- paste("FC by PTSD in Right lobe epi-", cur_region, "-", band, "band")
    figure <- annotate_figure(figure, top = text_grob(title_fig, face = "bold", size = 14))
    
    print(figure)
    svg(paste(path_RES, "Right_epi_only", paste(title_fig, ".svg", sep = ""), sep = "/"), width = 15, height = 15)
    print(figure)
    dev.off()
    print(figure)
    png(paste(path_RES, "Right_epi_only", paste(title_fig, ".png", sep = ""), sep = "/"), width = 25, height = 25, units = "cm", res = 700)
    print(figure)
    dev.off()
  
  }
  
  
}


################################################################################
#####           graphs (all single channels for each subj)                ######
################################################################################
allregions <- c("all channels", "Orbito-frontal-cortex", "Amygdala", "Hippocampus", "Anterior-cingulate-cortex", "Rhinal-cortex", "Insula")
#allregions <- "all channels"

for (i in 1:length(allregions)) {
  cur_region <- allregions[i]

if (cur_region == "all channels") {
  sel <- fc
} else{
  sel <- fc[str_detect(fc$brain_area, cur_region),]
}

sel <- sel[!is.na(sel$subject), ]
sel$chan_selection <- factor(sel$chan_selection, levels = c("TOTrest_all", "TOTrest_inv_inv", "TOTrest_inv_ni"))    

ipsisel <- sel[sel$involved == "ipsi" & sel$chan_selection != "TOTrest_ni_ni" & sel$chan_selection != "TOTrest_ni_inv" & sel$node_strength != 0,]
ipsisel$chan_selection <- factor(ipsisel$chan_selection, levels = c("TOTrest_all", "TOTrest_inv_inv", "TOTrest_inv_ni"))
ipsisel <- ipsisel[!is.na(ipsisel$subject), ]
contrasel <- sel[sel$involved == "contra" & sel$chan_selection != "TOTrest_inv_inv" & sel$chan_selection != "TOTrest_inv_ni" & sel$node_strength != 0,]
contrasel$chan_selection <- factor(contrasel$chan_selection, levels = c("TOTrest_all", "TOTrest_ni_ni", "TOTrest_ni_inv"))
contrasel <- contrasel[!is.na(contrasel$subject), ]


## boxplots PTSD+ vs PTSD-
rm(tmp)
rm(pwc)
rm(cur_data)

ymax = round(max(ipsisel$node_strength), 2) +0.03
ymin = 0

tmp <-ipsisel

# shapiro.test(tmp$node_strength)
pwc <- tmp %>% 
  group_by(chan_selection, involved) %>%
  pairwise_wilcox_test(node_strength ~ PTSDgroup, p.adjust.method = "fdr", paired = FALSE)
pwc <- pwc %>% add_xy_position(scales = "free_y")

pwc$y.position <- rep(ymax-0.02, nrow(pwc))
coco <- tmp %>%  #|d| < 0.2 "negligible", |d| < 0.5 "small", |d| < 0.8 "medium", otherwise "large"
  group_by(chan_selection, involved) %>%
  cohens_d(node_strength ~ PTSDgroup, paired = F)  
pwc <- cbind(pwc, coco[,4], coco[,9])
pwc$p <- round(pwc$p, 3)
pwc$effsize <- round(pwc$effsize, 1)


titre = paste("Ipsi-lateral channels")
ipsigraph <- ggboxplot(tmp, x = "PTSDgroup", y = "node_strength",
                      fill = "PTSDgroup", facet.by = "chan_selection", outlier.shape = NA) +
  #geom_line(aes(group = channel), alpha = 0.5, colour = "grey") +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 0.7, alpha = 0.6, colour = 'gray20')+
  stat_pvalue_manual(pwc, hide.ns = FALSE, label =  "{p.adj.signif}: p={p.adj}, d={effsize}", bracket.nudge.y = 0, tip.length = 0.005) +
  xlab("group") + ylab("Node strength") + ggtitle(titre) +
  theme_classic() + 
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  theme(axis.title = element_text(size = 10, face = "bold"), 
        axis.text = element_text(colour = "black", size = 10),
        axis.text.x = element_blank()) +
  theme(legend.text = element_text(size = 10),
        legend.title=element_blank(),
        legend.margin = margin(1,1,1,1, unit = "mm"),
        legend.direction = "horizontal",
        legend.position = "bottom")+
  scale_y_continuous(breaks = seq(ymin, ymax, 0.05), limits = c(ymin, ymax)) 
  # scale_fill_manual(values = c( "#cce3de", "#f6d5f7"))

## contra
rm(tmp)
rm(pwc)
ymax = round(max(contrasel$node_strength), 2) +0.03

tmp <-contrasel
pwc <- tmp %>% 
  group_by(chan_selection, involved) %>%
  pairwise_wilcox_test(node_strength ~ PTSDgroup, p.adjust.method = "fdr", paired = FALSE)
pwc <- pwc %>% add_xy_position(scales = "free_y")
pwc$y.position <- rep(ymax-0.02, nrow(pwc))
coco <- tmp %>%  #|d| < 0.2 "negligible", |d| < 0.5 "small", |d| < 0.8 "medium", otherwise "large"
  group_by(chan_selection, involved) %>%
  cohens_d(node_strength ~ PTSDgroup, paired = F)  
pwc <- cbind(pwc, coco[,4], coco[,9])
pwc$p <- round(pwc$p, 3)
pwc$effsize <- round(pwc$effsize, 1)

titre = paste("Contra-lateral channels")
contragraph <- ggboxplot(tmp, x = "PTSDgroup", y = "node_strength",
                       fill = "PTSDgroup", facet.by = "chan_selection", outlier.shape = NA) +
  #geom_line(aes(group = channel), alpha = 0.5, colour = "grey") +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 0.7, alpha = 0.6, colour = 'gray20')+
  stat_pvalue_manual(pwc, hide.ns = FALSE, label =  "{p.adj.signif}: p={p.adj}, d={effsize}", bracket.nudge.y = 0, tip.length = 0.005) +
  xlab("group") + ylab("Node strength") + ggtitle(titre) +
  theme_classic() + 
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  theme(axis.title = element_text(size = 10, face = "bold"), 
        axis.text = element_text(colour = "black", size = 10),
        axis.text.x = element_blank()) +
  theme(legend.text = element_text(size = 10),
        legend.title=element_blank(),
        legend.margin = margin(1,1,1,1, unit = "mm"),
        legend.direction = "horizontal",
        legend.position = "bottom")+
  scale_y_continuous(breaks = seq(ymin, ymax, 0.05), limits = c(ymin, ymax)) 
# scale_fill_manual(values = c( "#cce3de", "#f6d5f7"))

print(contragraph)
#### PLOT ALL TOGETHER
title_fig <- paste("Node strength by PTSD group in", cur_region, "-", band, "band")

figure <- ggarrange(ipsigraph, contragraph,  nrow = 2, labels = c("A", "B"), 
                    common.legend = TRUE, legend = "bottom")
figure <- annotate_figure(figure, 
                          top = text_grob(title_fig, 
                                          face = "bold", size = 14))

print(figure)
svg(paste(path_RES, band, paste(title_fig, ".svg", sep = ""), sep = "/"), width = 15, height = 15)
print(figure)
dev.off()
print(figure)
png(paste(path_RES, band, paste(title_fig, ".png", sep = ""), sep = "/"), width = 25, height = 25, units = "cm", res = 700)
print(figure)
dev.off()

################################ LEFT vs RIGHT  

sel <- sel[!is.na(sel$subject), ]
sel$chan_selection <- factor(sel$chan_selection, levels = c("TOTrest_all", "TOTrest_inv_inv", "TOTrest_inv_ni"))   
  
leftsel  <- sel[str_detect(sel$brain_area, "Left") & str_detect(sel$EpilepsySide, "L")
                & sel$chan_selection != "TOTrest_ni_ni" & sel$chan_selection != "TOTrest_ni_inv",]
leftsel  <- leftsel[!is.na(leftsel$subject), ]
rightsel <- sel[str_detect(sel$brain_area, "Right") & str_detect(sel$EpilepsySide, "R")
                & sel$chan_selection != "TOTrest_ni_ni" & sel$chan_selection != "TOTrest_ni_inv",]
rightsel <- rightsel[!is.na(rightsel$subject), ]


## boxplots PTSD+ vs PTSD-
rm(tmp)
rm(pwc)

tmp <-rightsel
ymax = round(max(tmp$node_strength), 2) +0.03
ymin = 0
# shapiro.test(tmp$node_strength)
pwc <- tmp %>% 
  group_by(chan_selection, involved) %>%
  pairwise_wilcox_test(node_strength ~ PTSDgroup, p.adjust.method = "fdr", paired = FALSE)
pwc <- pwc %>% add_xy_position(scales = "free_y")
pwc$y.position <- rep(ymax-0.02, nrow(pwc))
coco <- tmp %>%  #|d| < 0.2 "negligible", |d| < 0.5 "small", |d| < 0.8 "medium", otherwise "large"
  group_by(chan_selection, involved) %>%
  cohens_d(node_strength ~ PTSDgroup, paired = F)  
pwc <- cbind(pwc, coco[,4], coco[,9])
pwc$p <- round(pwc$p, 3)
pwc$effsize <- round(pwc$effsize, 1)


titre = paste("Ipsi channels in RIGHT-lobe epilepsy")
ipsigraph <- ggboxplot(tmp, x = "PTSDgroup", y = "node_strength",
                       fill = "PTSDgroup", facet.by = "chan_selection", outlier.shape = NA) +
  #geom_line(aes(group = channel), alpha = 0.5, colour = "grey") +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 0.7, 
              alpha = 0.6, colour = 'gray20')+
  #stat_pvalue_manual(pwc, hide.ns = FALSE, label =  "{p.adj.signif}: p={p.adj}, d={effsize}", 
                    # bracket.nudge.y = 0, tip.length = 0.005) +
  xlab("group") + ylab("Node strength") + ggtitle(titre) +
  theme_classic() + 
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  theme(axis.title = element_text(size = 10, face = "bold"), 
        axis.text = element_text(colour = "black", size = 10),
        axis.text.x = element_blank()) +
  theme(legend.text = element_text(size = 10),
        legend.title=element_blank(),
        legend.margin = margin(1,1,1,1, unit = "mm"),
        legend.direction = "horizontal",
        legend.position = "bottom")+
  scale_y_continuous(breaks = seq(ymin, ymax, 0.05), limits = c(ymin, ymax)) 
# scale_fill_manual(values = c( "#cce3de", "#f6d5f7"))

## contra
rm(tmp)
rm(pwc)

tmp <-leftsel
ymax = round(max(tmp$node_strength), 2) +0.03
pwc <- tmp %>% 
  group_by(chan_selection, involved) %>%
  pairwise_wilcox_test(node_strength ~ PTSDgroup, p.adjust.method = "fdr", paired = FALSE)
pwc <- pwc %>% add_xy_position(scales = "free_y")
pwc$y.position <- rep(ymax-0.02, nrow(pwc))
coco <- tmp %>%  #|d| < 0.2 "negligible", |d| < 0.5 "small", |d| < 0.8 "medium", otherwise "large"
  group_by(chan_selection, involved) %>%
  cohens_d(node_strength ~ PTSDgroup, paired = F)  
pwc <- cbind(pwc, coco[,4], coco[,9])
pwc$p <- round(pwc$p, 3)
pwc$effsize <- round(pwc$effsize, 1)

titre = paste("Ipsi channels in LEFT-lobe epilepsy")
contragraph <- ggboxplot(tmp, x = "PTSDgroup", y = "node_strength",
                         fill = "PTSDgroup", facet.by = "chan_selection", outlier.shape = NA) +
  #geom_line(aes(group = channel), alpha = 0.5, colour = "grey") +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 0.7, 
              alpha = 0.6, colour = 'gray20')+
  stat_pvalue_manual(pwc, hide.ns = FALSE, label =  "{p.adj.signif}: p={p.adj}, d={effsize}",
                     bracket.nudge.y = 0, tip.length = 0.005) +
  xlab("group") + ylab("Node strength") + ggtitle(titre) +
  theme_classic() + 
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  theme(axis.title = element_text(size = 10, face = "bold"), 
        axis.text = element_text(colour = "black", size = 10),
        axis.text.x = element_blank()) +
  theme(legend.text = element_text(size = 10),
        legend.title=element_blank(),
        legend.margin = margin(1,1,1,1, unit = "mm"),
        legend.direction = "horizontal",
        legend.position = "bottom")+
  scale_y_continuous(breaks = seq(ymin, ymax, 0.05), limits = c(ymin, ymax)) 
# scale_fill_manual(values = c( "#cce3de", "#f6d5f7"))
print(contragraph)

#### PLOT ALL TOGETHER
title_fig <- paste("Left vs Right epilepsy - node strength by PTSD group in", cur_region, "-", band, "band")

figure <- ggarrange(ipsigraph, contragraph,  nrow = 2, labels = c("A", "B"), 
                    common.legend = TRUE, legend = "bottom")
figure <- annotate_figure(figure, 
                          top = text_grob(title_fig, 
                                          face = "bold", size = 14))

print(figure)
svg(paste(path_RES, band, paste(title_fig, ".svg", sep = ""), sep = "/"), width = 15, height = 15)
print(figure)
dev.off()
print(figure)
png(paste(path_RES, band, paste(title_fig, ".png", sep = ""), sep = "/"), width = 25, height = 25, units = "cm", res = 700)
print(figure)
dev.off()

}



############ MIXED MODELS  #######################
# first with PTSD group only
fc_all <- fc[fc$chan_selection == "TOTrest_all",]
model1 <- lmer(node_strength ~ PTSDgroup * EpilepsySide * involved + (1|subject) + (1|Sex), data = fc_all) #only one band at the time
summary(model1)
anova(model1)
step(model1)
model2 <- lmer(node_strength ~ PTSDgroup + EpilepsySide + involved + (1 | subject) + PTSDgroup:EpilepsySide + 
                 PTSDgroup:involved + EpilepsySide:involved + PTSDgroup:EpilepsySide:involved, data = fc_all)
emmeans(model2, pairwise ~ PTSDgroup:EpilepsySide:involved, lmer.df = "satterthwaite", adjust = "tukey")

### only in AMY 
model1 <- lmer(node_strength ~ PTSDgroup * EpilepsySide * brain_area + (1|subject) + (1|Sex), data = amy[amy$chan_selection=="TOTrest_all",]) #only one band at the time
summary(model1)
anova(model1)
step(model1)
model2 <- lmer(node_strength ~ PTSDgroup + EpilepsySide + brain_area + (1 | subject) + PTSDgroup:EpilepsySide + PTSDgroup:brain_area + 
                 EpilepsySide:brain_area + PTSDgroup:EpilepsySide:brain_area, data = amy[amy$chan_selection=="TOTrest_all",])
emmeans(model2, pairwise ~ PTSDgroup:EpilepsySide:brain_area, lmer.df = "satterthwaite", adjust = "tukey")


################################################################################
##############################     OLD       ############################## 
################################################################################

### division contra/ipsi
contra <- fc[fc$involved == "contra" & fc$chan_selection != "TOTrest_inv_inv" & fc$chan_selection != "TOTrest_inv_ni" & fc$node_strength != 0,,]
contra$chan_selection <- factor(contra$chan_selection, levels = c("TOTrest_all", "TOTrest_ni_ni", "TOTrest_ni_inv"))

ipsi <- fc[fc$involved == "ipsi" & fc$chan_selection != "TOTrest_ni_ni" & fc$chan_selection != "TOTrest_ni_inv" & fc$node_strength != 0,]
ipsi$chan_selection <- factor(ipsi$chan_selection, levels = c("TOTrest_all", "TOTrest_inv_inv", "TOTrest_inv_ni"))

#### division left/right
leftfc <- fc[str_detect(fc$brain_area, "Left") & str_detect(fc$EpilepsySide, "L") & fc$chan_selection != "TOTrest_ni_ni" & fc$chan_selection != "TOTrest_ni_inv",]
rightfc <- fc[str_detect(fc$brain_area, "Right") & str_detect(fc$EpilepsySide, "R") & fc$chan_selection != "TOTrest_ni_ni" & fc$chan_selection != "TOTrest_ni_inv",]
leftfc$chan_selection <- factor(leftfc$chan_selection, levels = c("TOTrest_all", "TOTrest_inv_inv", "TOTrest_inv_ni"))    
rightfc$chan_selection <- factor(rightfc$chan_selection, levels = c("TOTrest_all", "TOTrest_inv_inv", "TOTrest_inv_ni"))    

#### select only contacts in hpc and amygdala
amyhpc <- fc[(str_detect(fc$brain_area, "Amygdala") | str_detect(fc$brain_area, "Hippocampus")) & fc$chan_selection != "TOTrest_inv_ni" & fc$chan_selection != "TOTrest_ni_inv",]
amyhpc <- amyhpc[!is.na(amyhpc$subject), ]
amyhpc <- cbind(amyhpc, struc=NA)
amyhpc$struc[str_detect(amyhpc$brain_area, "Amygdala") & is.na(amyhpc$brain_area) == FALSE] <- "Amygdala"
amyhpc$struc[str_detect(amyhpc$brain_area, "Hippocampus") & is.na(amyhpc$brain_area) == FALSE] <- "Hippocampus"

#### select only contacts in amygdala
amy <- fc[(str_detect(fc$brain_area, "Amygdala")) & fc$chan_selection != "TOTrest_inv_ni" & fc$chan_selection != "TOTrest_ni_inv",]
amy <- amy[!is.na(amy$subject), ]
amy <- cbind(amy, struc=NA)
amy$struc[str_detect(amy$brain_area, "Amygdala") & is.na(amy$brain_area) == FALSE] <- "Amygdala"

###############################################################################
colmag <- 8
allregions <- c("all channels", "Orbito-frontal-cortex", "Amygdala", "Hippocampus", "Anterior-cingulate-cortex", "Rhinal-cortex", "Insula")

for (i in 1:length(allregions)) {
  cur_region <- allregions[i]
  
  if (cur_region == "all channels") {
    sel <- fc
  } else{
    sel <- fc[str_detect(fc$brain_area, cur_region),]
  }
  
  sel <- sel[!is.na(sel$subject), ]
  sel$chan_selection <- factor(sel$chan_selection, levels = c("TOTrest_all", "TOTrest_inv_inv", "TOTrest_inv_ni"))    
  
  ipsisel <- sel[sel$involved == "ipsi" | sel$involved == "inv" & sel$chan_selection != "TOTrest_ni_ni" & sel$chan_selection != "TOTrest_ni_inv" & sel$node_strength != 0,]
  ipsisel$chan_selection <- factor(ipsisel$chan_selection, levels = c("TOTrest_all", "TOTrest_inv_inv", "TOTrest_inv_ni"))
  ipsisel <- ipsisel[!is.na(ipsisel$subject), ]
  contrasel <- sel[sel$involved == "contra" | sel$involved == "non-inv"  & sel$chan_selection != "TOTrest_inv_inv" & sel$chan_selection != "TOTrest_inv_ni" & sel$node_strength != 0,]
  contrasel$chan_selection <- factor(contrasel$chan_selection, levels = c("TOTrest_all", "TOTrest_ni_ni", "TOTrest_ni_inv"))
  contrasel <- contrasel[!is.na(contrasel$subject), ]
  
  rm(tmp)
  rm(pwc)
  rm(cur_data)
  ipsisel <- ipsisel[!is.na(ipsisel$node_strength), ]
  ymax = round(max(ipsisel$node_strength), 2) +0.03
  ymin = 0
  
  tmp <-ipsisel
  cur_data <- aggregate(node_strength ~ chan_selection + subject + PTSDgroup, 
                        filter_at(tmp, "involved", ~ .x == "ipsi"),
                        function(x) c(mean = mean(x), sd = sd(x))) 
  
  cur_data <- data.frame(PTSDgroup = cur_data$PTSDgroup, 
                         chan_selection = cur_data$chan_selection, 
                         subject = cur_data$subject, 
                         node_strength = cur_data$node_strength[, "mean"], 
                         node_strength.sd = cur_data$node_strength[, "sd"])
  pwc <- cur_data %>% 
    group_by(chan_selection) %>%
    pairwise_wilcox_test(node_strength ~ PTSDgroup, p.adjust.method = "fdr", paired = FALSE)
  pwc <- pwc %>% add_xy_position(scales = "free_y")
  
  pwc$y.position <- rep(ymax-0.02, nrow(pwc))
  coco <- cur_data %>%  #|d| < 0.2 "negligible", |d| < 0.5 "small", |d| < 0.8 "medium", otherwise "large"
    group_by(chan_selection) %>%
    cohens_d(node_strength ~ PTSDgroup, paired = F)  
  pwc <- cbind(pwc, coco[,4], coco[,colmag])
  pwc$p <- round(pwc$p, 3)
  pwc$effsize <- round(pwc$effsize, 1)
  
  
  titre = paste("Ipsi-lateral channels")
  ipsigraph <- ggboxplot(cur_data, x = "PTSDgroup", y = "node_strength",
                         fill = "PTSDgroup", facet.by = "chan_selection", outlier.shape = NA) +
    #geom_line(aes(group = channel), alpha = 0.5, colour = "grey") +
    geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 0.7, alpha = 0.6, colour = 'gray20')+
    stat_pvalue_manual(pwc, hide.ns = FALSE, label =  "{p.adj.signif}: p={p.adj}, d={effsize}", bracket.nudge.y = 0, tip.length = 0.005) +
    xlab("group") + ylab("Node strength") + ggtitle(titre) +
    theme_classic() + 
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(axis.title = element_text(size = 10, face = "bold"), 
          axis.text = element_text(colour = "black", size = 10),
          axis.text.x = element_blank()) +
    theme(legend.text = element_text(size = 10),
          legend.title=element_blank(),
          legend.margin = margin(1,1,1,1, unit = "mm"),
          legend.direction = "horizontal",
          legend.position = "bottom")+
    scale_y_continuous(breaks = seq(ymin, ymax, 0.05), limits = c(ymin, ymax)) 
  # scale_fill_manual(values = c( "#cce3de", "#f6d5f7"))
  
}

################################################################################
#############         1 mean value per subject             #####################
################################################################################
colmag <- 8
allregions <- c("all channels", "Orbito-frontal-cortex", "Amygdala", "Hippocampus", "Anterior-cingulate-cortex", "Rhinal-cortex", "Insula")

for (i in 1:length(allregions)) {
  cur_region <- allregions[i]
  
  if (cur_region == "all channels") {
    sel <- fc
  } else{
    sel <- fc[str_detect(fc$brain_area, cur_region),]
  }
  
  sel <- sel[!is.na(sel$subject), ]
  sel$chan_selection <- factor(sel$chan_selection, levels = c("TOTrest_all", "TOTrest_inv_inv", "TOTrest_inv_ni"))    
  
  ipsisel <- sel[sel$involved == "ipsi" | sel$involved == "inv" & sel$chan_selection != "TOTrest_ni_ni" & sel$chan_selection != "TOTrest_ni_inv" & sel$node_strength != 0,]
  ipsisel$chan_selection <- factor(ipsisel$chan_selection, levels = c("TOTrest_all", "TOTrest_inv_inv", "TOTrest_inv_ni"))
  ipsisel <- ipsisel[!is.na(ipsisel$subject), ]
  contrasel <- sel[sel$involved == "contra" | sel$involved == "non-inv"  & sel$chan_selection != "TOTrest_inv_inv" & sel$chan_selection != "TOTrest_inv_ni" & sel$node_strength != 0,]
  contrasel$chan_selection <- factor(contrasel$chan_selection, levels = c("TOTrest_all", "TOTrest_ni_ni", "TOTrest_ni_inv"))
  contrasel <- contrasel[!is.na(contrasel$subject), ]
  
  rm(tmp)
  rm(pwc)
  rm(cur_data)
  ipsisel <- ipsisel[!is.na(ipsisel$node_strength), ]
  ymax = round(max(ipsisel$node_strength), 2) +0.03
  ymin = 0
  
  tmp <-ipsisel
  cur_data <- aggregate(node_strength ~ chan_selection + subject + PTSDgroup, 
                        filter_at(tmp, "involved", ~ .x == "ipsi"),
                        function(x) c(mean = mean(x), sd = sd(x))) 
  
  cur_data <- data.frame(PTSDgroup = cur_data$PTSDgroup, 
                         chan_selection = cur_data$chan_selection, 
                         subject = cur_data$subject, 
                         node_strength = cur_data$node_strength[, "mean"], 
                         node_strength.sd = cur_data$node_strength[, "sd"])
  pwc <- cur_data %>% 
    group_by(chan_selection) %>%
    pairwise_wilcox_test(node_strength ~ PTSDgroup, p.adjust.method = "fdr", paired = FALSE)
  pwc <- pwc %>% add_xy_position(scales = "free_y")
  
  pwc$y.position <- rep(ymax-0.02, nrow(pwc))
  coco <- cur_data %>%  #|d| < 0.2 "negligible", |d| < 0.5 "small", |d| < 0.8 "medium", otherwise "large"
    group_by(chan_selection) %>%
    cohens_d(node_strength ~ PTSDgroup, paired = F)  
  pwc <- cbind(pwc, coco[,4], coco[,colmag])
  pwc$p <- round(pwc$p, 3)
  pwc$effsize <- round(pwc$effsize, 1)
  
  
  titre = paste("Ipsi-lateral channels")
  ipsigraph <- ggboxplot(cur_data, x = "PTSDgroup", y = "node_strength",
                         fill = "PTSDgroup", facet.by = "chan_selection", outlier.shape = NA) +
    #geom_line(aes(group = channel), alpha = 0.5, colour = "grey") +
    geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 0.7, alpha = 0.6, colour = 'gray20')+
    stat_pvalue_manual(pwc, hide.ns = FALSE, label =  "{p.adj.signif}: p={p.adj}, d={effsize}", bracket.nudge.y = 0, tip.length = 0.005) +
    xlab("group") + ylab("Node strength") + ggtitle(titre) +
    theme_classic() + 
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(axis.title = element_text(size = 10, face = "bold"), 
          axis.text = element_text(colour = "black", size = 10),
          axis.text.x = element_blank()) +
    theme(legend.text = element_text(size = 10),
          legend.title=element_blank(),
          legend.margin = margin(1,1,1,1, unit = "mm"),
          legend.direction = "horizontal",
          legend.position = "bottom")+
    scale_y_continuous(breaks = seq(ymin, ymax, 0.05), limits = c(ymin, ymax)) 
  # scale_fill_manual(values = c( "#cce3de", "#f6d5f7"))
  
  ################## contra
  rm(tmp)
  rm(pwc)
  
  contrasel <- contrasel[!is.na(contrasel$node_strength), ]
  ymax = round(max(contrasel$node_strength), 2) +0.03
  
  tmp <-contrasel
  cur_data <- aggregate(node_strength ~ chan_selection + subject + PTSDgroup, 
                        filter_at(tmp, "involved", ~ .x == "contra"),
                        function(x) c(mean = mean(x), sd = sd(x))) 
  
  cur_data <- data.frame(PTSDgroup = cur_data$PTSDgroup, 
                         chan_selection = cur_data$chan_selection, 
                         subject = cur_data$subject, 
                         node_strength = cur_data$node_strength[, "mean"], 
                         node_strength.sd = cur_data$node_strength[, "sd"])
  
  pwc <- cur_data %>% 
    group_by(chan_selection) %>%
    pairwise_wilcox_test(node_strength ~ PTSDgroup, p.adjust.method = "fdr", paired = FALSE)
  pwc <- pwc %>% add_xy_position(scales = "free_y")
  pwc$y.position <- rep(ymax-0.02, nrow(pwc))
  coco <- cur_data %>%  #|d| < 0.2 "negligible", |d| < 0.5 "small", |d| < 0.8 "medium", otherwise "large"
    group_by(chan_selection) %>%
    cohens_d(node_strength ~ PTSDgroup, paired = F)  
  pwc <- cbind(pwc, coco[,4], coco[,colmag])
  pwc$p <- round(pwc$p, 3)
  pwc$effsize <- round(pwc$effsize, 1)
  
  titre = paste("Contra-lateral channels")
  contragraph <- ggboxplot(cur_data, x = "PTSDgroup", y = "node_strength",
                           fill = "PTSDgroup", facet.by = "chan_selection", outlier.shape = NA) +
    #geom_line(aes(group = channel), alpha = 0.5, colour = "grey") +
    geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 0.7, alpha = 0.6, colour = 'gray20')+
    stat_pvalue_manual(pwc, hide.ns = FALSE, label =  "{p.adj.signif}: p={p.adj}, d={effsize}", bracket.nudge.y = 0, tip.length = 0.005) +
    xlab("group") + ylab("Node strength") + ggtitle(titre) +
    theme_classic() + 
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(axis.title = element_text(size = 10, face = "bold"), 
          axis.text = element_text(colour = "black", size = 10),
          axis.text.x = element_blank()) +
    theme(legend.text = element_text(size = 10),
          legend.title=element_blank(),
          legend.margin = margin(1,1,1,1, unit = "mm"),
          legend.direction = "horizontal",
          legend.position = "bottom")+
    scale_y_continuous(breaks = seq(ymin, ymax, 0.05), limits = c(ymin, ymax)) 
  # scale_fill_manual(values = c( "#cce3de", "#f6d5f7"))
  
  
  #### PLOT ALL TOGETHER
  title_fig <- paste("Mean node str by subj in", cur_region, "-", band, "band")
  
  figure <- ggarrange(ipsigraph, contragraph,  nrow = 2, labels = c("A", "B"), 
                      common.legend = TRUE, legend = "bottom")
  figure <- annotate_figure(figure, 
                            top = text_grob(title_fig, 
                                            face = "bold", size = 14))
  
  print(figure)
  # svg(paste(path_RES, band, paste(title_fig, ".svg", sep = ""), sep = "/"), width = 15, height = 15)
  # print(figure)
  # dev.off()
  # print(figure)
  png(paste(path_RES, band, paste(title_fig, ".png", sep = ""), sep = "/"), width = 25, height = 25, units = "cm", res = 700)
  print(figure)
  dev.off()
  
  ################################ LEFT vs RIGHT  
  
  sel <- sel[!is.na(sel$subject), ]
  sel$chan_selection <- factor(sel$chan_selection, levels = c("TOTrest_all", "TOTrest_inv_inv", "TOTrest_inv_ni"))   
  
  leftsel  <- sel[str_detect(sel$brain_area, "Left") & str_detect(sel$EpilepsySide, "L") 
                  & sel$chan_selection != "TOTrest_ni_ni" & sel$chan_selection != "TOTrest_ni_inv",]
  leftsel  <- leftsel[!is.na(leftsel$subject), ]
  rightsel <- sel[str_detect(sel$brain_area, "Right") & str_detect(sel$EpilepsySide, "R") 
                  & sel$chan_selection != "TOTrest_ni_ni" & sel$chan_selection != "TOTrest_ni_inv",]
  rightsel <- rightsel[!is.na(rightsel$subject), ]
  
  
  ## boxplots PTSD+ vs PTSD-
  rm(tmp)
  rm(pwc)
  
  tmp <-rightsel
  cur_data <- aggregate(node_strength ~ chan_selection + subject + PTSDgroup, 
                        filter_at(tmp, "involved", ~ .x == "ipsi"),
                        function(x) c(mean = mean(x), sd = sd(x))) 
  
  cur_data <- data.frame(PTSDgroup = cur_data$PTSDgroup, 
                         chan_selection = cur_data$chan_selection, 
                         subject = cur_data$subject, 
                         node_strength = cur_data$node_strength[, "mean"], 
                         node_strength.sd = cur_data$node_strength[, "sd"])
  
  cur_data <- cur_data[!is.na(cur_data$node_strength), ]
  ymax = round(max(cur_data$node_strength), 2) +0.03
  ymin = 0
  # shapiro.test(cur_data$node_strength)
  pwc <- cur_data %>% 
    group_by(chan_selection) %>%
    pairwise_wilcox_test(node_strength ~ PTSDgroup, p.adjust.method = "fdr", paired = FALSE)
  pwc <- pwc %>% add_xy_position(scales = "free_y")
  pwc$y.position <- rep(ymax-0.02, nrow(pwc))
  coco <- cur_data %>%  #|d| < 0.2 "negligible", |d| < 0.5 "small", |d| < 0.8 "medium", otherwise "large"
    group_by(chan_selection) %>%
    cohens_d(node_strength ~ PTSDgroup, paired = F)  
  pwc <- cbind(pwc, coco[,4], coco[,colmag])
  pwc$p <- round(pwc$p, 3)
  pwc$effsize <- round(pwc$effsize, 1)
  
  
  titre = paste("Ipsi channels in RIGHT-lobe epilepsy")
  ipsigraph <- ggboxplot(cur_data, x = "PTSDgroup", y = "node_strength",
                         fill = "PTSDgroup", facet.by = "chan_selection", outlier.shape = NA) +
    #geom_line(aes(group = channel), alpha = 0.5, colour = "grey") +
    geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 0.7, 
                alpha = 0.6, colour = 'gray20')+
    stat_pvalue_manual(pwc, hide.ns = FALSE, label =  "{p.adj.signif}: p={p.adj}, d={effsize}", 
                       bracket.nudge.y = 0, tip.length = 0.005) +
    xlab("group") + ylab("Node strength") + ggtitle(titre) +
    theme_classic() + 
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(axis.title = element_text(size = 10, face = "bold"), 
          axis.text = element_text(colour = "black", size = 10),
          axis.text.x = element_blank()) +
    theme(legend.text = element_text(size = 10),
          legend.title=element_blank(),
          legend.margin = margin(1,1,1,1, unit = "mm"),
          legend.direction = "horizontal",
          legend.position = "bottom")+
    scale_y_continuous(breaks = seq(ymin, ymax, 0.05), limits = c(ymin, ymax)) 
  # scale_fill_manual(values = c( "#cce3de", "#f6d5f7"))
  
  ## left
  rm(tmp)
  rm(pwc)
  
  tmp <-leftsel
  cur_data <- aggregate(node_strength ~ chan_selection + subject + PTSDgroup, 
                        filter_at(tmp, "involved", ~ .x == "ipsi"),
                        function(x) c(mean = mean(x), sd = sd(x))) 
  
  cur_data <- data.frame(PTSDgroup = cur_data$PTSDgroup, 
                         chan_selection = cur_data$chan_selection, 
                         subject = cur_data$subject, 
                         node_strength = cur_data$node_strength[, "mean"], 
                         node_strength.sd = cur_data$node_strength[, "sd"])
  
  cur_data <- cur_data[!is.na(cur_data$node_strength), ]
  ymax = round(max(cur_data$node_strength), 2) +0.03
  pwc <- cur_data %>% 
    group_by(chan_selection) %>%
    pairwise_wilcox_test(node_strength ~ PTSDgroup, p.adjust.method = "fdr", paired = FALSE)
  pwc <- pwc %>% add_xy_position(scales = "free_y")
  pwc$y.position <- rep(ymax-0.02, nrow(pwc))
  coco <- cur_data %>%  #|d| < 0.2 "negligible", |d| < 0.5 "small", |d| < 0.8 "medium", otherwise "large"
    group_by(chan_selection) %>%
    cohens_d(node_strength ~ PTSDgroup, paired = F)  
  pwc <- cbind(pwc, coco[,4], coco[,colmag])
  pwc$p <- round(pwc$p, 3)
  pwc$effsize <- round(pwc$effsize, 1)
  
  titre = paste("Ipsi channels in LEFT-lobe epilepsy")
  contragraph <- ggboxplot(cur_data, x = "PTSDgroup", y = "node_strength",
                           fill = "PTSDgroup", facet.by = "chan_selection", outlier.shape = NA) +
    #geom_line(aes(group = channel), alpha = 0.5, colour = "grey") +
    geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 0.7, 
                alpha = 0.6, colour = 'gray20')+
    stat_pvalue_manual(pwc, hide.ns = FALSE, label =  "{p.adj.signif}: p={p.adj}, d={effsize}",
                       bracket.nudge.y = 0, tip.length = 0.005) +
    xlab("group") + ylab("Node strength") + ggtitle(titre) +
    theme_classic() + 
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(axis.title = element_text(size = 10, face = "bold"), 
          axis.text = element_text(colour = "black", size = 10),
          axis.text.x = element_blank()) +
    theme(legend.text = element_text(size = 10),
          legend.title=element_blank(),
          legend.margin = margin(1,1,1,1, unit = "mm"),
          legend.direction = "horizontal",
          legend.position = "bottom")+
    scale_y_continuous(breaks = seq(ymin, ymax, 0.05), limits = c(ymin, ymax)) 
  # scale_fill_manual(values = c( "#cce3de", "#f6d5f7"))
  
  
  #### PLOT ALL TOGETHER
  title_fig <- paste("Left vs Right-mean node str by subj in", cur_region, "-", band, "band")
  
  figure <- ggarrange(ipsigraph, contragraph,  nrow = 2, labels = c("A", "B"), 
                      common.legend = TRUE, legend = "bottom")
  figure <- annotate_figure(figure, 
                            top = text_grob(title_fig, 
                                            face = "bold", size = 14))
  
  print(figure)
  # svg(paste(path_RES, band, paste(title_fig, ".svg", sep = ""), sep = "/"), width = 15, height = 15)
  # print(figure)
  # dev.off()
  # print(figure)
  png(paste(path_RES, band, paste(title_fig, ".png", sep = ""), sep = "/"), width = 25, height = 25, units = "cm", res = 700)
  print(figure)
  dev.off()
  
}


